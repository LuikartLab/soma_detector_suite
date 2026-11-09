"""
cellpose_engine.py — Cellpose-SAM detection engine for the Soma Detector (v3.0)

Drop-in alternative to the classical threshold/watershed detector.  It produces
the same soma dictionaries (cy, cx, area_um2, cell_mean, peak_intensity, mask2d,
z, z_min, z_max, source) so curation, measurement/export, NPZ round-trips, the
Cre classifier and Group Analysis are untouched.

Pipeline
  1. Cellpose-SAM run in 2-D on every z-plane of a channel (GPU if available).
  2. Labels linked through z by mutual-best mask overlap -> one track per soma.
  3. Each slice outline cleaned: morphological opening removes dendrite stubs,
     the largest component is kept, optional breadth erode/dilate.
  4. Representative slice = brightest slice of the track (same rule as the
     classical deduplicate_3d).  Tracks closer than ~half a soma are merged.
  5. Per-CELL selection: contrast (mean / plane background) and z-extent.
     Cellpose finds essentially every labelled cell; these sliders decide which
     are shown for curation, replacing the classical pixel-threshold sigma.

Model inference is cached per channel so the selection sliders re-run in
seconds; only the model settings (cell-probability threshold, cell diameter)
trigger a new inference.

Cellpose 4.x + PyTorch are imported lazily: the GUI still runs without them
(engine reports unavailable) — launch from the .venv-cellpose environment.
"""
import os
import time
from collections import defaultdict

import numpy as np
from scipy import ndimage
from skimage import measure, morphology

_HERE = os.path.dirname(os.path.abspath(__file__))

# ── model files ────────────────────────────────────────────────────────────────
# Cellpose-SAM base weights (cpsam_v2) shipped in models\cellpose: point Cellpose
# there so a freshly copied install needs no download (and no internet).  Must be
# set before cellpose is first imported anywhere.
_LOCAL_CP_MODELS = os.path.join(_HERE, 'models', 'cellpose')
if os.path.isdir(_LOCAL_CP_MODELS):
    os.environ.setdefault('CELLPOSE_LOCAL_MODELS_PATH', _LOCAL_CP_MODELS)

# Fine-tuned weights: Cellpose-SAM trained 2026-09-23 on 349 crops around the
# lab's hand-circled Fig 2B somas (human polygon = label).  On 6 held-out stacks
# (96 cells) vs base cpsam: area CCC 0.90 vs 0.71, 91 % vs 72 % of cells within
# ±25 % of the human area, plane within ±1 slice 93 % vs 83 %, recall equal.
# None = base cpsam_v2.  The fine-tuned model draws the human boundary itself,
# so with it the outline-breadth default is 0 µm (base cpsam needs −0.5 µm).
CUSTOM_MODEL = 'cpsam_luikart_v2'
CUSTOM_MODEL_PATH = os.path.join(_HERE, 'models', CUSTOM_MODEL) if CUSTOM_MODEL else None
if CUSTOM_MODEL_PATH and not os.path.isfile(CUSTOM_MODEL_PATH):
    CUSTOM_MODEL_PATH = None          # file missing → base model (model_name() says so)

# ── optional human-pick ranking model ("refiner") ──────────────────────────────
# Trained on the lab's hand-circled ROI sets: P(a human would have circled this
# track).  Cellpose finds every labelled cell; this ranks them.  Absent file →
# the contrast / z-extent rules alone apply.
PICK_MODEL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models', 'soma_pick_model.joblib')
_PICK = None
_PICK_TRIED = False


def pick_model():
    global _PICK, _PICK_TRIED
    if not _PICK_TRIED:
        _PICK_TRIED = True
        try:
            import joblib
            if os.path.isfile(PICK_MODEL_PATH):
                _PICK = joblib.load(PICK_MODEL_PATH)
        except Exception:
            _PICK = None
    return _PICK


def track_features(t, nz, px_um, vol, dmax, med_contrast):
    """Per-track features for the pick model (must match build_pick_dataset)."""
    rep = t['rep']; R0, C0, R1, C1 = rep['bbox']; crop = rep['crop']
    rp = measure.regionprops(crop.astype(np.uint8))[0]
    per = rp.perimeter if rp.perimeter > 0 else 1e-6
    vals = vol[rep['z'], R0:R1, C0:C1][crop]
    zs = [s['z'] for s in t['slices']]; areas = [s['area_um2'] for s in t['slices']]
    return {'contrast': t['contrast'], 'rep_contrast': rep['contrast'], 'contrast_rel': t['contrast'] / med_contrast,
            'mean': rep['mean'], 'peak': rep['peak'], 'sat_frac': float((vals >= dmax).mean()),
            'area': rep['area_um2'], 'area_max': max(areas), 'area_vs_max': rep['area_um2'] / max(areas),
            'n_slices': t['n_slices'], 'z_rel': rep['z'] / max(nz - 1, 1),
            'rep_pos': zs.index(rep['z']) / max(len(zs) - 1, 1),
            'touch_surface': int(t['z_min'] == 0 or t['z_max'] == nz - 1), 'cut_by_surface': int(_cut_by_surface(t, nz)),
            'solidity': float(rp.solidity), 'compactness': 4 * np.pi * rp.area / per ** 2, 'eccentricity': float(rp.eccentricity),
            'minor_um': rp.minor_axis_length * px_um, 'major_um': rp.major_axis_length * px_um,
            'border_um': min(rep['cy'], rep['cx'], vol.shape[1] - 1 - rep['cy'], vol.shape[2] - 1 - rep['cx']) * px_um}


def score_pool(pool, vol, px_um, ch):
    """Attach t['pick'] = P(human would circle it) to every track of a pool."""
    pm = pick_model()
    if pm is None or not pool:
        return False
    nz = vol.shape[0]
    dmax = np.iinfo(vol.dtype).max if np.issubdtype(vol.dtype, np.integer) else float(vol.max())
    med_c = float(np.median([t['contrast'] for t in pool]))
    reps = np.array([[t['rep']['cy'], t['rep']['cx'], t['rep']['z']] for t in pool])
    feats = []
    for t in pool:
        f = track_features(t, nz, px_um, vol, dmax, med_c)
        d = np.hypot(reps[:, 0] - t['rep']['cy'], reps[:, 1] - t['rep']['cx']) * px_um
        near_z = np.abs(reps[:, 2] - t['rep']['z']) <= 3
        f['n_nb_10um'] = int(((d < 10) & near_z).sum()) - 1
        f['n_nb_20um'] = int(((d < 20) & near_z).sum()) - 1
        f['ch'] = ch
        feats.append([float(f[k]) for k in pm['features']])
    p = pm['model'].predict_proba(np.array(feats))[:, 1]
    for t, v in zip(pool, p):
        t['pick'] = float(v)
    return True

_MODEL      = None
_AVAILABLE  = None
_IMPORT_ERR = ''
_DEVICE     = ''

MODEL_DEFAULTS = {'cellprob': 0.0, 'flow': 0.4, 'diam': 0.0}   # diam 0 = model default (no rescale)
SELECT_DEFAULTS = {'min_contrast': 3.0, 'min_slices': 2, 'open_um': 1.8,
                   'breadth_um': -0.5,        # 1-px erosion: Cellpose outlines include the halo/boundary
                                              # ring; -0.5 µm puts mean area within ~5 % of hand-circled
                   'min_area_um2': 25.0, 'max_area_um2': 600.0,
                   'min_pick': 0.3,           # ranking model: P(human would circle) — out-of-fold on 76 stacks,
                                              # 0.3 keeps 92 % of hand-circled cells and shows ~150/stack of ~620
                   'exclude_z_edges': True,   # soma touching first/last plane = cut by the section surface
                   'surface_keep_um2': 100.0, # ... unless its largest visible cross-section is already this
                                              # big (a soma-sized cell cut at the surface is still offered;
                                              # only small surface fragments are dropped).  Fig 1C iv + 2B,
                                              # 718 hand-circled cells (2026-09-28): the plain rule hid 33 of
                                              # them (median 159 um2); 100 um2 recovers 27 for ~7 extra
                                              # candidates per image; None = plain rule.
                   'rep_rule': 'adaptive'}    # which plane of a track represents the soma (see pick_rep)
# rep_rule: 713 hand-circled Fig 1C iv + Fig 2B somas with the fine-tuned model (rules_eval, 2026-09-28):
#   rule         |dz|<=1   area ratio   CCC    200-250 um2   >=250 um2
#   brightest      44 %      0.89       0.69       -            -        (peak intensity sits ~1.6 planes
#   largest        81 %      1.05       0.87                             above the human's plane)
#   mean_int       87 %      0.98       0.87      89 %         90 %     (old default; picks a small bright
#   adaptive       89 %      1.00       0.90      98 %         98 %  <- end plane of tall somas)
# 'adaptive' = highest mean intensity among the planes whose cross-section is within
# a fraction of the track's largest (fraction 0 below 120 um2, 0.6 to 180 um2, 0.8 above).
if CUSTOM_MODEL_PATH:
    SELECT_DEFAULTS['breadth_um'] = 0.0        # fine-tuned model: no erosion (see CUSTOM_MODEL)
REP_RULES = {                                  # key -> label shown in the GUI
    'brightest':   'brightest plane (classical rule)',
    'largest':     'largest cross-section',
    'largest_mid': 'largest cross-section, middle of track',
    'mid':         'middle plane of track',
    'mean_int':    'highest mean intensity',
    'adaptive':    'brightest plane near the largest cross-section (size-adaptive)',
}
MERGE_XY_UM = 4.0      # tracks whose representative centroids are closer than
MERGE_Z_GAP = 2        # this (and within this many slices) are the same soma


# ── availability ───────────────────────────────────────────────────────────────

def available():
    """True when cellpose + torch import.  Checked once, lazily."""
    global _AVAILABLE, _IMPORT_ERR, _DEVICE
    if _AVAILABLE is None:
        try:
            import torch
            from cellpose import models  # noqa: F401  # pylint: disable=unused-import
            _AVAILABLE = True
            _DEVICE = (torch.cuda.get_device_name(0) if torch.cuda.is_available()
                       else 'CPU only — slow')
        except Exception as e:           # ImportError, DLL errors, ...
            _AVAILABLE  = False
            _IMPORT_ERR = str(e)
    return _AVAILABLE


def device_name():
    available()
    return _DEVICE


def import_error():
    available()
    return _IMPORT_ERR


def get_model():
    """Load Cellpose-SAM once: the fine-tuned weights when CUSTOM_MODEL is set
    and present, else the base cpsam_v2 (from models\\cellpose or ~/.cellpose)."""
    global _MODEL
    if _MODEL is None:
        from cellpose import models
        if CUSTOM_MODEL_PATH:
            _MODEL = models.CellposeModel(gpu=True, pretrained_model=CUSTOM_MODEL_PATH)
        else:
            _MODEL = models.CellposeModel(gpu=True)
    return _MODEL


def model_name():
    """Which weights detect(): recorded in every NPZ and shown in the GUI."""
    return CUSTOM_MODEL if CUSTOM_MODEL_PATH else 'cpsam_v2 (base)'


# ── step 1: inference ──────────────────────────────────────────────────────────

def segment_planes(vol, model_params, batch_size=8):
    """Cellpose-SAM independently on every plane of a (Z, Y, X) volume.
    Returns an int32 (Z, Y, X) label stack; labels are per-plane (not linked)."""
    model = get_model()
    vol   = np.asarray(vol).astype(np.float32)
    if vol.ndim == 2:
        vol = vol[np.newaxis]
    # Whole-stack percentile normalisation (Cellpose would normalise each plane
    # separately for a list input, which amplifies noise in empty top/bottom
    # planes).  Same 1–99 % rule Cellpose uses internally.
    lo, hi = np.percentile(vol, [1.0, 99.0])
    if hi <= lo:
        hi = lo + 1.0
    vol   = np.clip((vol - lo) / (hi - lo), 0.0, 1.0)
    diam  = float(model_params.get('diam') or 0) or None
    planes = [vol[z] for z in range(vol.shape[0])]          # independent 2-D images
    masks = model.eval(planes, batch_size=batch_size, diameter=diam,
                       flow_threshold=float(model_params.get('flow', 0.4)),
                       cellprob_threshold=float(model_params.get('cellprob', 0.0)),
                       do_3D=False, normalize=False,
                       channel_axis=None, z_axis=None)[0]
    masks = np.stack([np.asarray(m) for m in masks]).astype(np.int32)
    return masks


# ── step 2: link labels through z ──────────────────────────────────────────────

def link_planes(masks, min_overlap=0.5):
    """Union-find over (z, label) nodes.  Consecutive planes are linked only
    when the two labels are each other's best overlap partner AND the overlap
    is at least `min_overlap` of the smaller label — one-to-one links, so
    vertically stacked neighbours are not chained together.
    Returns a list of tracks, each a z-sorted list of (z, label)."""
    nz     = masks.shape[0]
    areas  = [np.bincount(masks[z].ravel()) for z in range(nz)]
    parent = {}

    def find(a):
        while parent.get(a, a) != a:
            parent[a] = parent.get(parent[a], parent[a])
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for z in range(nz - 1):
        a, b = masks[z].ravel(), masks[z + 1].ravel()
        sel  = (a > 0) & (b > 0)
        if not sel.any():
            continue
        nb   = int(b.max()) + 1
        pair = a[sel].astype(np.int64) * nb + b[sel].astype(np.int64)
        uniq, cnt = np.unique(pair, return_counts=True)
        la, lb = uniq // nb, uniq % nb
        best_f, best_b = {}, {}
        for x, y, c in zip(la.tolist(), lb.tolist(), cnt.tolist()):
            if c > best_f.get(x, (0, 0))[0]:
                best_f[x] = (c, y)
            if c > best_b.get(y, (0, 0))[0]:
                best_b[y] = (c, x)
        for x, (c, y) in best_f.items():
            if best_b.get(y, (0, -1))[1] != x:
                continue
            if c < min_overlap * min(areas[z][x], areas[z + 1][y]):
                continue
            union((z, x), (z + 1, y))

    tracks = defaultdict(list)
    for z in range(nz):
        for lab in np.nonzero(areas[z])[0]:
            if lab == 0:
                continue
            tracks[find((z, int(lab)))].append((z, int(lab)))
    return [sorted(t) for t in tracks.values()]


# ── step 3/4: clean outlines, build per-track records ─────────────────────────

def _clean_crop(crop, r_open, breadth_px):
    m = crop
    if r_open > 0:
        m = morphology.binary_opening(m, morphology.disk(r_open))
        if not m.any():
            return m
        lab, n = ndimage.label(m)
        if n > 1:                                   # opening can split: keep the largest piece
            sizes = ndimage.sum(m, lab, range(1, n + 1))
            m = lab == (int(np.argmax(sizes)) + 1)
    if breadth_px > 0:
        m = ndimage.binary_fill_holes(morphology.binary_dilation(m, morphology.disk(breadth_px)))
    elif breadth_px < 0:
        e = morphology.binary_erosion(m, morphology.disk(-breadth_px))
        if e.any():
            m = e
    return m


def pick_rep(slices, rule='brightest'):
    """Representative plane of a z-sorted list of cleaned slice records.
    Brightness falls with depth, so 'brightest' is biased toward the top of
    each soma; the human convention is the largest in-focus cross-section."""
    if rule == 'largest':
        return max(slices, key=lambda s: s['area_um2'])
    if rule == 'largest_mid':                      # ignore the blurred ends of the track
        n = len(slices)
        lo, hi = int(np.floor(n * 0.2)), int(np.ceil(n * 0.8))
        return max(slices[lo:hi] or slices, key=lambda s: s['area_um2'])
    if rule == 'mid':
        return slices[len(slices) // 2]
    if rule == 'mean_int':
        return max(slices, key=lambda s: s['mean'])
    if rule == 'adaptive':
        # Mean intensity alone favours the small bright end plane of a tall soma
        # (26 um2 chosen where the human circled 186); restrict the choice to
        # planes near the largest cross-section, more strictly for bigger cells.
        a_max = max(s['area_um2'] for s in slices)
        frac  = 0.0 if a_max < 120 else (0.6 if a_max < 180 else 0.8)
        band  = [s for s in slices if s['area_um2'] >= frac * a_max]
        return max(band, key=lambda s: s['mean'])
    return max(slices, key=lambda s: s['peak'])    # 'brightest'


def build_pool(masks, vol, px_um, tracks, sel):
    """Per-track records with cleaned per-slice crops and statistics.
    Nothing is filtered here except per-slice size / border rules, so the
    pool can be re-selected with different contrast / z-extent settings and
    searched by force-detect."""
    nz, ny, nx = masks.shape
    r_open     = int(round(float(sel.get('open_um', 0.0)) / px_um))
    b_um       = float(sel.get('breadth_um', 0.0))
    breadth_px = int(round(abs(b_um) / px_um)) * (1 if b_um > 0 else -1) if b_um else 0
    margin     = r_open + abs(breadth_px) + 1
    min_a, max_a = float(sel.get('min_area_um2', 25.0)), float(sel.get('max_area_um2', 600.0))
    bg   = np.array([np.percentile(vol[z], 20) for z in range(nz)], float) + 1.0
    bbox = [{rp.label: rp.bbox for rp in measure.regionprops(masks[z])} for z in range(nz)]

    pool = []
    for tr in tracks:
        slices = []
        for z, lab in tr:
            r0, c0, r1, c1 = bbox[z][lab]
            R0, C0 = max(0, r0 - margin), max(0, c0 - margin)
            R1, C1 = min(ny, r1 + margin), min(nx, c1 + margin)
            crop = masks[z, R0:R1, C0:C1] == lab
            crop = _clean_crop(crop, r_open, breadth_px)
            npx  = int(crop.sum())
            if npx == 0:
                continue
            area = npx * px_um ** 2
            if area < min_a or area > max_a:
                continue
            rr, cc = np.nonzero(crop)
            if (R0 + rr.min() == 0 or C0 + cc.min() == 0 or
                    R0 + rr.max() == ny - 1 or C0 + cc.max() == nx - 1):
                continue                                    # touches the image border
            vals = vol[z, R0:R1, C0:C1][crop]
            slices.append({'z': z, 'bbox': (R0, C0, R1, C1), 'crop': crop,
                           'area_um2': area, 'mean': float(vals.mean()),
                           'peak': float(vals.max()),
                           'cy': float(R0 + rr.mean()), 'cx': float(C0 + cc.mean()),
                           'contrast': float(vals.mean()) / bg[z]})
        if not slices:
            continue
        # The first and last planes (section surfaces) never represent a soma:
        # the representative is chosen among interior planes, and a track that
        # exists only on a boundary plane is dropped by select_tracks.
        inner = [s for s in slices if 0 < s['z'] < nz - 1]
        rep = pick_rep(inner or slices, sel.get('rep_rule', 'brightest'))
        pool.append({'slices': slices, 'rep': rep, 'n_slices': len(slices),
                     'edge_only': not inner,
                     'z_min': slices[0]['z'], 'z_max': slices[-1]['z'],
                     'contrast': max(s['contrast'] for s in slices)})
    return pool


def materialize(track, ny, nx, source='auto', rep=None):
    """Track record -> soma dict in the classical format (full-size mask2d).
    `rep` selects which slice record represents the soma (default: the
    track's brightest slice)."""
    rep = rep or track['rep']
    R0, C0, R1, C1 = rep['bbox']
    mask2d = np.zeros((ny, nx), bool)
    mask2d[R0:R1, C0:C1] = rep['crop']
    return {'cy': rep['cy'], 'cx': rep['cx'], 'area_um2': rep['area_um2'],
            'cell_mean': rep['mean'], 'peak_intensity': rep['peak'],
            'mask2d': mask2d, 'z': rep['z'],
            'z_min': track['z_min'], 'z_max': track['z_max'],
            'n_slices': track['n_slices'], 'contrast': track['contrast'],
            'pick': track.get('pick'), 'engine': 'cellpose', 'source': source}


# 1.0 = drop only when the boundary plane IS the largest cross-section (soma still
# growing when the section surface cut it).  Measured on 444 hand-circled Fig 2B
# somas: 'touches a surface' would drop 37 % of them; this rule drops 7 % of them
# and 6 % of all labels; 0.9 would drop 11 % / 8 %.
Z_EDGE_FRACTION = 1.0


def _cut_by_surface(track, nz, keep_um2=None):
    """True if the track reaches a stack surface and its cross-section on that
    boundary plane is still >= Z_EDGE_FRACTION of the track maximum — unless
    that maximum is already >= keep_um2 (a soma-sized cell cut at the surface
    is still worth offering; only small surface fragments are hidden)."""
    if track['z_min'] > 0 and track['z_max'] < nz - 1:
        return False
    areas = {s['z']: s['area_um2'] for s in track['slices']}
    a_max = max(areas.values())
    if keep_um2 is not None and a_max >= float(keep_um2):
        return False
    edge  = max(areas.get(0, 0.0), areas.get(nz - 1, 0.0))
    return edge >= Z_EDGE_FRACTION * a_max


def select_tracks(pool, sel, px_um, roi_mask=None, nz=None):
    """Apply the per-cell selection rules, then merge near-duplicate tracks
    (a soma whose z-linking broke into two pieces)."""
    keep = [t for t in pool
            if not t.get('edge_only')
            and t['n_slices'] >= int(sel.get('min_slices', 1))
            and t['contrast'] >= float(sel.get('min_contrast', 0.0))
            and t.get('pick', 1.0) >= float(sel.get('min_pick', 0.0))]
    if nz is not None and sel.get('exclude_z_edges', True):
        # Surface-cut somas: a track that reaches the first or last plane AND
        # is still growing there (boundary-plane area >= Z_EDGE_FRACTION of the
        # track's largest cross-section) was cut before its maximum diameter.
        # Tracks that merely touch a surface but peak well inside the stack are
        # kept — the hand-circled sets contain many of those (36 % of the Fig 2B
        # somas touch a surface; only the still-growing ones are dropped).
        keep = [t for t in keep
                if not _cut_by_surface(t, nz, sel.get('surface_keep_um2'))]
    if roi_mask is not None:
        keep = [t for t in keep
                if roi_mask[int(round(t['rep']['cy'])), int(round(t['rep']['cx']))]]
    if len(keep) > 1:
        cy = np.array([t['rep']['cy'] for t in keep]); cx = np.array([t['rep']['cx'] for t in keep])
        zlo = np.array([t['z_min'] for t in keep]); zhi = np.array([t['z_max'] for t in keep])
        pk = np.array([t['rep']['peak'] for t in keep])
        d  = np.hypot(cy[:, None] - cy[None, :], cx[:, None] - cx[None, :]) * px_um
        zc = (zlo[:, None] <= zhi[None, :] + MERGE_Z_GAP) & (zlo[None, :] <= zhi[:, None] + MERGE_Z_GAP)
        close = (d < MERGE_XY_UM) & zc
        np.fill_diagonal(close, False)
        drop = np.zeros(len(keep), bool)
        order = np.argsort(-pk)                                  # brightest wins
        for i in order:
            if drop[i]:
                continue
            for j in np.nonzero(close[i])[0]:
                if not drop[j] and pk[j] <= pk[i]:
                    drop[j] = True
        keep = [t for t, dr in zip(keep, drop) if not dr]
    return keep


# ── engine object used by the GUI ──────────────────────────────────────────────

class CellposeDetector:
    """Caches inference (per channel key + model params) and the cleaned pool
    (+ outline params) so slider changes are fast; supports force-detect
    look-ups into tracks that the selection rules excluded."""

    def __init__(self):
        self._masks = {}     # ch_key -> (model_key, masks)
        self._pool  = {}     # ch_key -> (pool_key, tracks, pool)

    def clear(self):
        self._masks.clear(); self._pool.clear()

    def detect(self, vol, px_um, model_params, sel, roi_mask=None, ch_key=0, status=None):
        vol = np.asarray(vol)
        nz, ny, nx = vol.shape
        info = {'seconds_model': 0.0, 'cached': True}
        mkey = tuple(sorted((k, float(v)) for k, v in model_params.items()))
        if self._masks.get(ch_key, (None,))[0] != mkey:
            if status:
                status('Cellpose-SAM: segmenting {} planes…'.format(nz))
            t0 = time.time()
            masks = segment_planes(vol, model_params)
            info['seconds_model'] = time.time() - t0
            info['cached'] = False
            self._masks[ch_key] = (mkey, masks)
            self._pool.pop(ch_key, None)
        masks = self._masks[ch_key][1]

        pkey = (mkey, float(sel.get('open_um', 0)), float(sel.get('breadth_um', 0)),
                float(sel.get('min_area_um2', 0)), float(sel.get('max_area_um2', 0)),
                str(sel.get('rep_rule', 'brightest')))
        if self._pool.get(ch_key, (None,))[0] != pkey:
            if status:
                status('Cellpose-SAM: linking planes and cleaning outlines…')
            tracks = link_planes(masks)
            pool   = build_pool(masks, vol, px_um, tracks, sel)
            if status:
                status('Ranking cells…')
            score_pool(pool, vol, px_um, int(ch_key) if isinstance(ch_key, int) else 0)
            self._pool[ch_key] = (pkey, tracks, pool)
        pool = self._pool[ch_key][2]

        kept  = select_tracks(pool, sel, px_um, roi_mask, nz=nz)
        somas = [materialize(t, ny, nx) for t in kept]
        info.update({'n_tracks': len(pool), 'n_selected': len(somas)})
        return somas, info

    def lookup(self, ch_key, z, py, px, ny, nx, z_tol=2):
        """Force-detect: return the pooled track containing pixel (py, px) on
        a plane within ±z_tol of z, or None.  The soma is represented by the
        track's outline ON THE CLICKED PLANE (nearest available), not the
        brightest slice — the user chose that plane deliberately.  Used to
        re-add a cell the selection rules excluded, or to move a soma to the
        plane the user would have circled it on (reject, then force-detect)."""
        if ch_key not in self._pool:
            return None
        _, _, pool = self._pool[ch_key]
        best = None
        for t in pool:
            for s in t['slices']:
                if abs(s['z'] - z) > z_tol:
                    continue
                R0, C0, R1, C1 = s['bbox']
                if R0 <= py < R1 and C0 <= px < C1 and s['crop'][py - R0, px - C0]:
                    dz = abs(s['z'] - z)
                    if best is None or dz < best[0]:
                        best = (dz, t, s)
        return materialize(best[1], ny, nx, source='manual', rep=best[2]) if best else None
