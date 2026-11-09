"""
soma_area_and_intensity_gui_3.0.py
------------------------------------
v3.0 — Cellpose-SAM edition (the classical threshold/watershed engine lives
on in the v2.x series; it has been removed here).

  * Detection engine: cellpose_engine.py (Cellpose-SAM per plane -> z-linked
    tracks -> cleaned outlines -> per-CELL selection by contrast / z-extent /
    surface rule).  On the Brain-paper replication set (445 hand-circled somas,
    26 stacks) the defaults recover 94.6 % within 4 µm (classical: 63 %) with a
    median area ratio to the human outline of 1.00.
  * Control panel fits without scrolling: the selection sliders re-select from
    the cached model output instantly and carry curation over; rarely used
    settings (z-extent, area range, surface rule, representative plane, and the
    ⟳ model settings) live in "Advanced settings…".
  * Opening an image forces the File-info dialog (channel assignment by
    emission/excitation band, animal ID, mouse genotype), pre-filled from the
    filename and the previous image; OK starts detection.
  * Zoom is centred on the cursor, clamped to the image and resets to the full
    image when zoomed out; middle-drag pans.  Outline editing grabs the NEAREST
    handle in screen pixels, a click on the outline inserts a handle, right-
    click removes one.  Force-detect: toggle, or hold F while clicking; in
    Cellpose mode it returns the model's own outline on the plane you are on.

v2.1 additions (carried forward):

  1. Min/Max pixel statistics  — image window always shows % of pixels at
     zero (minimum) and at dtype-max (saturated) for the current view.

  2. Min/Max mask overlay  — toggle button overlays dark-blue on zero pixels
     and hot-pink on saturated pixels so clipping is immediately visible.

  3. Per-plane intensity scaling  — optional clip-and-rescale applied to
     every image plane and channel independently.  The clip percentile X is
     user-controlled (default 1 %).  Affects both the display and the scaled
     intensity export.

  4. Dual export  — "Export measurements" produces two Excel files:
       soma_measurements.xlsx               (raw pixel values, unchanged)
       intensity_scaled_soma_measurements.xlsx  (per-plane scaled values)
     (The previous histogram-matched export has been removed.)

  5. Click-to-jump IPC  — when launched alongside the Group Analysis GUI,
     a TCP listener on localhost:5099 accepts "open_cell" commands that
     auto-load the file, jump to the slice, and flash the matching ROI.

v2.0 improvements (carried forward):
  - Distance-transform peak fallback, per-slice adaptive sigma, evidence-
    based shape thresholds derived from 2055-row training dataset.
  - Force-detect mode, manual soma tracking, full diagnostic CSV export.
"""

import sys, os, warnings, threading, glob, re, json, time
import cell_jump_ipc
import cellpose_engine
import tempfile
import numpy as np
import czifile
import tifffile
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Polygon as MplPolygon
from skimage import measure, morphology
from skimage.draw import polygon as sk_polygon
# (skimage.exposure.match_histograms removed — histogram matching disabled)
from scipy import ndimage
import xml.etree.ElementTree as ET
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
import csv

warnings.filterwarnings('ignore')

# ── Constants ──────────────────────────────────────────────────────────────────
PIXEL_SIZE_UM_DEFAULT   = 0.4457
# Shape-filter defaults — calibrated against 7,949-soma training set
# (pS6-MLP n=5110 + pS6-Tariq n=2839, 72% accepted / 28% rejected).
# Conservative values: each threshold set no higher than the 5th percentile
# of accepted cells so that ≥95% of genuine somas (KO and WT alike) survive.
MIN_AREA_UM2            = 65     # 5th-pct of accepted: ~65 µm² (both channels)
MAX_AREA_UM2            = 400
MIN_SOLIDITY            = 0.80   # rejected median 0.83–0.84; accepted 25th-pct 0.87
MAX_ECCENTRICITY        = 0.93   # signal inconsistent between datasets — unchanged
MIN_COMPACTNESS         = 0.50   # rejected median 0.57–0.58; accepted 25th-pct 0.60
MIN_MINOR_AXIS_UM       = 7.5    # 5th-pct of accepted: ~7.4 µm; rounded to slider step
MIN_INSCRIBED_RADIUS_UM = 3.5    # modest discriminatory power — unchanged
BG_SIGMA_UM             = 6.68   # was BG_SIGMA_PX=15 @ 0.446 um/px
SMOOTH_SIGMA_UM         = 0.45   # was SMOOTH_SIGMA_PX=1
SUBPEAK_SMOOTH_UM       = 0.90   # was sigma=2.0 in peak-split smoothing
WS_MIN_DIST_UM          = 3.0    # was WS_MIN_DIST_PX=7
NMS_XY_RADIUS_UM        = 11.0   # was NMS_XY_RADIUS=25
NMS_Z_RADIUS            = 15     # in Z-slices (Z spacing already uniform)
CODET_Z_RADIUS          = 3   # tighter Z window for cross-channel co-detection

SENSITIVITY_MIN, SENSITIVITY_MAX = 0.5, 6.0
BREADTH_MIN,     BREADTH_MAX     = -2.0, 4.0

MANUAL_COLOUR = '#00E5FF'


# ── Detection helpers ──────────────────────────────────────────────────────────

def get_pixel_size(meta_xml):
    try:
        meta = ET.fromstring(meta_xml)
        for elem in meta.iter():
            if elem.get('Id') == 'X':
                val = float(elem.findtext('Value', '0'))
                if val > 0:
                    return val * 1e6
    except Exception:
        pass
    return PIXEL_SIZE_UM_DEFAULT


def _parse_detection_band(ch):
    """Parse a channel's emission detection band from its <Ranges> element
    (a nested descendant of the <Channel>, e.g. '461.9-543.2').
    Returns (lo_nm, hi_nm), or (0.0, 0.0) if unavailable."""
    el = ch.find('.//Ranges')
    txt = (el.text or '') if el is not None else ''
    m = re.search(r'(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)', txt)
    if not m:
        return 0.0, 0.0
    lo, hi = float(m.group(1)), float(m.group(2))
    return (lo, hi) if lo <= hi else (hi, lo)


def _channel_dict(ch):
    """Build a per-channel info dict from an acquired <Channel> element."""
    def _f(tag):
        s = (ch.findtext(tag, '') or '').strip()
        try:    return float(s)
        except ValueError: return 0.0
    def _f_deep(tag):
        el = ch.find('.//' + tag)
        s = (el.text or '').strip() if el is not None else ''
        try:    return float(s)
        except ValueError: return 0.0
    lo, hi = _parse_detection_band(ch)
    return {
        'name':     ch.get('Name', '') or '',
        'fluor':    ch.findtext('Fluor', '') or ch.findtext('DyeName', '') or '',
        'ex_nm':    _f('ExcitationWavelength'),
        'em_nm':    _f('EmissionWavelength'),
        'laser_nm': _f_deep('Wavelength'),   # actual laser line (nested descendant)
        'det_lo':   lo,
        'det_hi':   hi,
    }


def read_czi_channel_info(meta_xml):
    """Return per-channel info dicts for the ACQUIRED channels.

    Reads the Information/Image/Dimensions/Channels block, which lists exactly
    the acquired planes with real fluor names, excitation/emission peaks, the
    laser line, and the emission detection band (<Ranges>) — the physical
    hardware setting, independent of any (possibly wrong) channel label. Falls
    back to scanning every <Channel> element (deduped by name) for older files
    that lack that block.

    Each dict has: name, fluor, ex_nm, em_nm, laser_nm, det_lo, det_hi.
    """
    try:
        meta = ET.fromstring(meta_xml)
    except Exception:
        return []

    block = meta.find('.//Information/Image/Dimensions/Channels')
    if block is not None:
        chans = block.findall('Channel')
        if chans:
            return [_channel_dict(ch) for ch in chans]

    # Legacy fallback: every <Channel> in the document, deduped by name.
    channels, seen = [], set()
    for ch in meta.iter('Channel'):
        name = ch.get('Name', '')
        if name in seen:
            continue
        seen.add(name)
        channels.append(_channel_dict(ch))
    return channels


def format_channel_identifier(i, c):
    """Label-independent identifier for a physical channel, leading with the
    emission detection band (the true hardware setting), then the laser line."""
    lo, hi = c.get('det_lo', 0.0), c.get('det_hi', 0.0)
    ex = c.get('laser_nm', 0.0) or c.get('ex_nm', 0.0)
    if lo and hi:
        band = 'Em {:.0f}–{:.0f} nm'.format(lo, hi)
    elif lo:                       # long-pass filter (open upper bound)
        band = 'Em ≥{:.0f} nm'.format(lo)
    elif hi:                       # short-pass filter (open lower bound)
        band = 'Em ≤{:.0f} nm'.format(hi)
    elif c.get('em_nm', 0.0):
        band = 'Em ~{:.0f} nm'.format(c['em_nm'])
    else:
        band = c.get('name', '') or 'channel {}'.format(i)
    ex_txt = '  (Ex {:.0f})'.format(ex) if ex else ''
    return '{}: {}{}'.format(i, band, ex_txt)


def _parse_lsm_filter(name):
    """Parse an LSM emission filter name into (det_lo, det_hi, em_center) nm.

    Handles the common Zeiss forms: 'BP 505-550' / 'BP505-550' (band-pass),
    'BP 500/50' (centre/width), 'LP 650' (long-pass → (650, 0)), 'SP 550'
    (short-pass → (0, 550)). Returns (0.0, 0.0, 0.0) when nothing parses.
    """
    s = (name or '').upper().strip()
    if not s:
        return 0.0, 0.0, 0.0
    # Long-pass / short-pass first (so 'LP 650' isn't misread as a bare number).
    m = re.search(r'\bLP\s*(\d{3,4})', s)
    if m:
        lo = float(m.group(1))
        return lo, 0.0, lo
    m = re.search(r'\bSP\s*(\d{3,4})', s)
    if m:
        hi = float(m.group(1))
        return 0.0, hi, hi
    # Band-pass as an explicit range: 505-550.
    m = re.search(r'(\d{3,4})\s*[-–]\s*(\d{3,4})', s)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        if lo > hi:
            lo, hi = hi, lo
        return lo, hi, (lo + hi) / 2.0
    # Band-pass as centre/width: 500/50.
    m = re.search(r'(\d{3,4})\s*/\s*(\d{1,3})', s)
    if m:
        c, w = float(m.group(1)), float(m.group(2))
        return c - w / 2.0, c + w / 2.0, c
    return 0.0, 0.0, 0.0


# Floxed-allele genes recognized in filenames (shared vocabulary with the
# downstream group_analysis_gui genotype parser).
GENOTYPE_GENES = [('PTEN', 'Pten'), ('RAPTOR', 'Raptor'), ('RICTOR', 'Rictor'),
                  ('AKT1', 'Akt1'), ('AKT2', 'Akt2'), ('AKT3', 'Akt3')]


def genotype_to_richtext(runs):
    """Build an openpyxl CellRichText from (text, italic, superscript) runs,
    preserving italics and superscript for the Excel export. Returns None for
    empty input so the cell stays blank."""
    if not runs:
        return None
    blocks = []
    for text, italic, sup in runs:
        if not text:
            continue
        fnt = InlineFont(rFont='Arial', i=bool(italic))
        if sup:
            fnt.vertAlign = 'superscript'
        blocks.append(TextBlock(fnt, text))
    return CellRichText(blocks) if blocks else None


def guess_animal_id(filename):
    """Best-effort unique animal identifier from a filename.

    Recognizes the standard cage-card form (CC89603_111720#1) and a couple of
    common variants; falls back to a date-of-birth + sex + number token
    (e.g. 9-19-25-3F). Returns '' when nothing recognizable is found — the user
    then types the ID into the (always-editable) Animal ID field.
    """
    stem = os.path.splitext(os.path.basename(filename or ''))[0]
    # CC cage-card + date + animal#: CC89603_111720#1
    m = re.search(r'(CC\d+_\d+\s*#\s*\d+)', stem)
    if m:
        return m.group(1).replace(' ', '')
    m = re.search(r'(CC\d+)[ _-]+(\d{4,8})[ _-]*#?\s*(\d+)', stem)
    if m:
        return '{}_{}#{}'.format(m.group(1), m.group(2), m.group(3))
    m = re.search(r'(CC\d+)', stem)
    if m:
        return m.group(1)
    # Date-of-birth + number + sex, e.g. 9-19-25-3F
    m = re.search(r'(\d{1,2}-\d{1,2}-\d{2,4}-\d+[MF])', stem)
    if m:
        return m.group(1)
    return ''


def guess_genotype_runs(filename, allele='flx/flx'):
    """Guess a mouse genotype from a filename as formatted runs.

    Returns a list of (text, italic, superscript) tuples. For a filename whose
    genotype token contains 'PtenRictor':
        [('Pten', True, False), ('fl/fl', False, True), (';', False, False),
         ('Rictor', True, False), ('fl/fl', False, True)]
    Gene names are italicized; allele designations are superscript. Only the
    underscore-delimited token carrying a gene name is inspected, so reporter
    strings (e.g. 'pRubiGAKT1T2ACre') do not spawn spurious alleles.
    """
    stem = os.path.splitext(os.path.basename(filename or ''))[0]
    gt_tok = ''
    for tok in stem.split('_'):
        tu = tok.upper()
        if 'PTEN' in tu or 'RAPTOR' in tu or 'RICTOR' in tu:
            gt_tok = tu
            break
    if not gt_tok:
        if 'C57' in stem.upper():
            return [('C57BL/6', False, False)]
        return []
    genes = [disp for key, disp in GENOTYPE_GENES if key in gt_tok]
    if not genes:
        return []
    runs = []
    for i, g in enumerate(genes):
        if i:
            runs.append((';', False, False))
        runs.append((g, True, False))
        runs.append((allele, False, True))
    return runs


def auto_assign_channels(ch_info):
    """
    Map physical CZI channels to roles: 0=GFP-Cre, 1=mCherry, 2=pAntibody.
    Returns [gfp_idx, mcherry_idx, pantibody_idx] or None if assignment fails.
    """
    FLUOR_ROLES = {
        0: ['egfp', 'gfp', 'fitc', 'alexa488', 'af488'],
        1: ['dsred', 'mcherry', 'rfp', 'tritc', 'tdtomato', 'mruby',
            'alexa568', 'alexa555', 'rubi'],
        2: ['cy5', 'alexafluor647', 'af647', 'alexa647', 'a647'],
    }
    assigned = {}
    for pi, ch in enumerate(ch_info):
        fl = ch['fluor'].lower()
        ex = ch['ex_nm']
        role = None
        for r, kws in FLUOR_ROLES.items():
            if any(k in fl for k in kws):
                role = r
                break
        if role is None and ex > 0:
            if   470 <= ex <= 510: role = 0
            elif 540 <= ex <= 585: role = 1
            elif 620 <= ex <= 670: role = 2
        if role is not None and role not in assigned:
            assigned[role] = pi
    if set(assigned.keys()) >= {0, 1, 2}:
        return [assigned[0], assigned[1], assigned[2]]
    return None


def load_lsm(path):
    """
    Read a Zeiss LSM file and return (stack, px_um, ch_info).
    stack shape : (n_ch, nz, ny, nx)
    ch_info     : list of {name, fluor, ex_nm} — same format as read_czi_channel_info
    """
    with tifffile.TiffFile(path) as tif:
        data = tif.asarray()   # typically (Z, C, Y, X) for LSM

        px_um = PIXEL_SIZE_UM_DEFAULT
        lsm   = tif.lsm_metadata
        ch_info = []
        if lsm:
            vx = lsm.get('VoxelSizeX', 0)
            if vx and float(vx) > 0:
                px_um = float(vx) * 1e6

            # Build per-channel info from the acquisition tracks so the SAME
            # fields as read_czi_channel_info are available (name, fluor, ex_nm,
            # em_nm, laser_nm, det_lo, det_hi). LSM organises this per track:
            # excitation as illumination laser lines, emission as a filter name
            # (e.g. 'BP 505-550', 'LP 650'). ChannelNames is often absent, so we
            # must not key off it.
            cc = lsm.get('ChannelColors', {}) or {}
            color_names = cc.get('ColorNames', []) if isinstance(cc, dict) else []
            si_tracks = lsm.get('ScanInformation', {}).get('Tracks', []) or []
            try:
                for tr in si_tracks:
                    if not tr.get('Acquire', 1):
                        continue
                    # Primary excitation laser for the track (highest-power line).
                    laser, best_pow = 0.0, -1.0
                    for il in tr.get('IlluminationChannels', []):
                        wl = float(il.get('Wavelength', 0) or 0)
                        pw = float(il.get('Power', 0) or 0)
                        if wl > 0 and pw > best_pow:
                            best_pow, laser = pw, wl
                    if laser == 0.0:
                        for il in tr.get('IlluminationChannels', []):
                            wl = float(il.get('Wavelength', 0) or 0)
                            if wl > 0:
                                laser = wl
                                break
                    for dc in tr.get('DetectionChannels', []):
                        if not dc.get('Acquire', 1):
                            continue
                        lo, hi, em = _parse_lsm_filter(dc.get('FilterName', ''))
                        idx = len(ch_info)
                        nm = (color_names[idx] if idx < len(color_names)
                              else (dc.get('ChannelName') or 'Ch{}'.format(idx + 1)))
                        ch_info.append({
                            'name':     str(nm),
                            'fluor':    dc.get('DyeName', '') or '',
                            'ex_nm':    laser,
                            'em_nm':    em,
                            'laser_nm': laser,
                            'det_lo':   lo,
                            'det_hi':   hi,
                        })
            except Exception:
                ch_info = []

            # Fallback: no usable track info — at least emit channel names.
            if not ch_info:
                names = lsm.get('ChannelNames', []) or color_names
                for name in names:
                    ch_info.append({'name': str(name), 'fluor': str(name),
                                    'ex_nm': 0.0, 'em_nm': 0.0, 'laser_nm': 0.0,
                                    'det_lo': 0.0, 'det_hi': 0.0})

        axes = tif.series[0].axes.upper() if tif.series else ''

    # Remove time axis if present (take first time point)
    if 'T' in axes:
        t_idx = axes.index('T')
        data  = data.take(0, axis=t_idx)
        axes  = axes.replace('T', '')

    # Normalise to (C, Z, Y, X)
    if axes == 'CZYX':
        pass
    elif axes == 'ZCYX':
        data = data.transpose(1, 0, 2, 3)
    elif axes == 'ZYX':
        data = data[np.newaxis]           # single channel → (1, Z, Y, X)
    elif axes == 'CYX':
        data = data[:, np.newaxis]        # single slice   → (C, 1, Y, X)
    elif axes == 'YX':
        data = data[np.newaxis, np.newaxis]
    else:
        # Generic fallback: try to put C before Z based on index
        if data.ndim == 4:
            # Assume ZCYX (most common LSM layout)
            data = data.transpose(1, 0, 2, 3)
        elif data.ndim == 3:
            data = data[np.newaxis]

    # Pad ch_info to match actual channel count
    n_ch = data.shape[0]
    while len(ch_info) < n_ch:
        ch_info.append({'name': f'Ch{len(ch_info)+1}', 'fluor': '', 'ex_nm': 0.0,
                        'em_nm': 0.0, 'laser_nm': 0.0, 'det_lo': 0.0, 'det_hi': 0.0})
    ch_info = ch_info[:n_ch]

    return data, px_um, ch_info


# ── Diagnostics helpers ────────────────────────────────────────────────────────

def compute_soma_diagnostics(soma, stack, px_um, sigma_mult, ch_idx):
    mask2d = soma['mask2d']
    z      = soma['z']
    nz     = stack.shape[1]
    img    = stack[ch_idx, z].astype(np.float32)

    props = measure.regionprops(mask2d.astype(np.uint8), intensity_image=img)
    if not props:
        return {}
    rp  = props[0]
    per = float(rp.perimeter) if rp.perimeter > 0 else 1e-6
    comp = 4.0 * np.pi * rp.area / (per * per)
    dist = ndimage.distance_transform_edt(mask2d)

    smooth_px = max(0.5, SMOOTH_SIGMA_UM / px_um)
    bg_px     = max(2.0, BG_SIGMA_UM    / px_um)
    smoothed  = ndimage.gaussian_filter(img, sigma=smooth_px)
    bg_img    = ndimage.gaussian_filter(smoothed, sigma=bg_px)
    hp        = smoothed - bg_img
    hp      -= hp.min()
    hp_med   = float(np.median(hp))
    hp_mad   = float(np.median(np.abs(hp - hp_med)))
    if hp_mad < 1e-6:
        hp_mad = max(float(hp.std()) / 1.4826, 1e-6)
    thr = hp_med + sigma_mult * hp_mad * 1.4826

    cy_i = int(round(np.clip(soma['cy'], 0, img.shape[0] - 1)))
    cx_i = int(round(np.clip(soma['cx'], 0, img.shape[1] - 1)))
    hp_val = float(hp[cy_i, cx_i])

    bg_at_soma = float(bg_img[cy_i, cx_i])
    snr_raw    = round(float(img[mask2d].max()) / bg_at_soma, 3) \
                 if bg_at_soma > 0 else 0.0
    z_rel = round(z / max(nz - 1, 1), 3)

    return {
        'compactness':      round(comp, 4),
        'solidity':         round(float(rp.solidity), 4),
        'eccentricity':     round(float(rp.eccentricity), 4),
        'minor_axis_um':    round(float(rp.minor_axis_length) * px_um, 3),
        'major_axis_um':    round(float(rp.major_axis_length) * px_um, 3),
        'inscribed_r_um':   round(float(dist.max()) * px_um, 3),
        'peak_intensity':   round(float(img[mask2d].max()), 2),
        'mean_intensity':   round(float(rp.mean_intensity), 2),
        'bg_at_soma':       round(bg_at_soma, 2),
        'snr_raw':          snr_raw,
        'hp_mad':           round(hp_mad, 4),
        'hp_threshold':     round(thr, 4),
        'hp_val_at_center': round(hp_val, 4),
        'sigma_to_detect':  round(hp_val / (hp_mad * 1.4826) - hp_med / (hp_mad * 1.4826)
                                  if hp_mad > 0 else 0.0, 3),
        'z_rel':            z_rel,
    }


# ── Intensity scaling helper ───────────────────────────────────────────────────

def scale_plane_intensity(img_plane, pct, dtype_max=65535.0):
    """
    Clip the bottom and top `pct` percent of pixel values in a single 2-D
    plane, then rescale linearly so the surviving range maps to [0, dtype_max].
    Each plane is scaled independently, so a dim mid-stack plane gets the same
    output range as a bright surface plane.
    Returns float32.  Input may be any numeric dtype.
    """
    img    = img_plane.astype(np.float32)
    p_low  = float(np.percentile(img, pct))
    p_high = float(np.percentile(img, 100.0 - pct))
    if p_high <= p_low:
        return img
    img = np.clip(img, p_low, p_high)
    return ((img - p_low) / (p_high - p_low) * dtype_max).astype(np.float32)


# ── Standard Excel export ──────────────────────────────────────────────────────
# (Histogram matching removed — raw and intensity-scaled exports only.)

def batch_header_channel_names(pairs, default_names):
    """Column-header labels for a batch of (image, npz) pairs: per channel, the
    names found in the NPZ metadata across all files, sorted and joined with
    ' / ' when they differ (a file without metadata contributes default_names)."""
    names_by_ch = [[] for _ in range(len(default_names))]
    for _img, npz_path in pairs:
        names = None
        try:
            d = np.load(npz_path, allow_pickle=False)
            if 'meta_json' in d:
                names = json.loads(str(d['meta_json'])).get('ch_names')
        except Exception:
            names = None
        if not names or len(names) < len(default_names):
            names = default_names
        for i, n in enumerate(names[:len(default_names)]):
            n = str(n).strip()
            if n and n not in names_by_ch[i]:
                names_by_ch[i].append(n)
    return [' / '.join(sorted(ns)) if ns else default_names[i]
            for i, ns in enumerate(names_by_ch)]


def measure_and_export(stack, px_um, czi_fname, soma_by_ch, rejected_ids,
                       ch_names, out_path, scale_pct=None, genotype_runs=None,
                       animal_id=None, header_ch_names=None):
    """ch_names are THIS file's channel names (they go into the per-cell
    'Channel detected' / 'Source channel' values).  header_ch_names, when
    given, are the labels used for the workbook's channel columns instead:
    a batch directory can mix files whose channel names differ (GFP-AKT1-Cre
    vs GFP-AKT3-Cre), and the shared column headers must not pretend that
    one file's names apply to all of them."""
    n_ch, nz, ny, nx = stack.shape
    hdr_names = list(header_ch_names) if header_ch_names else list(ch_names)
    # Native dtype max — 255 for uint8, 65535 for uint16, etc.
    dtype_max = float(np.iinfo(stack.dtype).max) \
                if np.issubdtype(stack.dtype, np.integer) else 65535.0

    combined_mask = np.zeros((nz, ny, nx), bool)
    all_cells     = []

    # Cross-channel co-detection: (ch_idx, si) -> set of partner ch names
    # and (ch_idx, si) -> (partner_ch_idx, partner_si) for soma-number sharing.
    # Uses CODET_Z_RADIUS (3 slices) — tighter than the within-channel NMS_Z_RADIUS.
    _codet         = {}
    _codet_partner = {}
    _ch_keys = list(soma_by_ch.keys())
    for _i, _ch_a in enumerate(_ch_keys):
        for _ch_b in _ch_keys[_i + 1:]:
            for _sa_i, _sa in enumerate(soma_by_ch[_ch_a]):
                if (_ch_a, _sa_i) in rejected_ids:
                    continue
                for _sb_i, _sb in enumerate(soma_by_ch[_ch_b]):
                    if (_ch_b, _sb_i) in rejected_ids:
                        continue
                    if (np.hypot(_sa['cx'] - _sb['cx'], _sa['cy'] - _sb['cy']) < NMS_XY_RADIUS_UM / px_um
                            and abs(_sa['z'] - _sb['z']) <= CODET_Z_RADIUS):
                        _codet.setdefault((_ch_a, _sa_i), set()).add(ch_names[_ch_b])
                        _codet.setdefault((_ch_b, _sb_i), set()).add(ch_names[_ch_a])
                        _codet_partner[(_ch_a, _sa_i)] = (_ch_b, _sb_i)
                        _codet_partner[(_ch_b, _sb_i)] = (_ch_a, _sa_i)

    # Determine starting soma number — continues from max in existing workbook
    # so every soma across all files in the workbook has a unique number.
    soma_n = 1
    if os.path.isfile(out_path):
        try:
            _wb_peek = openpyxl.load_workbook(out_path, read_only=True, data_only=True)
            if 'Cells' in _wb_peek.sheetnames:
                _ws_peek = _wb_peek['Cells']
                # Locate the 'Soma #' column from the header row (row 2) rather
                # than assuming a fixed position — the column layout can change.
                _soma_ci = None
                for _hrow in _ws_peek.iter_rows(min_row=2, max_row=2,
                                                values_only=True):
                    for _ci, _hv in enumerate(_hrow):
                        if str(_hv).strip() == 'Soma #':
                            _soma_ci = _ci
                            break
                    break
                if _soma_ci is not None:
                    for _r in _ws_peek.iter_rows(min_row=3, values_only=True):
                        if _soma_ci < len(_r):
                            _v = _r[_soma_ci]
                            if isinstance(_v, (int, float)) and int(_v) >= soma_n:
                                soma_n = int(_v) + 1
            _wb_peek.close()
        except Exception:
            soma_n = 1

    # Pre-assign soma numbers: co-detected pairs share one number,
    # each solo soma gets its own.
    _soma_n_map = {}
    for ch_idx, soma_list in soma_by_ch.items():
        for si in range(len(soma_list)):
            if (ch_idx, si) in rejected_ids:
                continue
            if (ch_idx, si) in _soma_n_map:
                continue  # already assigned as the partner of an earlier soma
            partner = _codet_partner.get((ch_idx, si))
            if partner and partner in _soma_n_map:
                _soma_n_map[(ch_idx, si)] = _soma_n_map[partner]
            else:
                _soma_n_map[(ch_idx, si)] = soma_n
                soma_n += 1

    for ch_idx, soma_list in soma_by_ch.items():
        for si, s in enumerate(soma_list):
            if (ch_idx, si) in rejected_ids:
                continue
            z, mask = s['z'], s['mask2d']
            combined_mask[z] |= mask
            _partners = _codet.get((ch_idx, si), set())
            ch_det = ' + '.join(sorted([ch_names[ch_idx]] + list(_partners))) \
                     if _partners else ch_names[ch_idx]
            row = {
                'file':        czi_fname,
                'ch_detected': ch_det,
                'source_ch':   ch_names[ch_idx],
                'soma_id':     _soma_n_map[(ch_idx, si)],
                'cx_um':       round(s['cx'] * px_um, 2),
                'cy_um':       round(s['cy'] * px_um, 2),
                'area_um2':    round(s['area_um2'], 3),
                'slice':       z + 1,
            }
            for ch in range(n_ch):
                if scale_pct is not None:
                    plane = scale_plane_intensity(stack[ch, z], scale_pct, dtype_max)
                else:
                    plane = stack[ch, z].astype(np.float32)
                row['mean_ch{}'.format(ch + 1)] = round(float(plane[mask].mean()), 3)
                # Slice dynamic range always from raw data (p99 - p1), regardless of scaling
                raw_sl = stack[ch, z].astype(np.float32)
                row['range_ch{}'.format(ch + 1)] = round(
                    float(np.percentile(raw_sl, 99.0) - np.percentile(raw_sl, 1.0)), 1)
            all_cells.append(row)

    # Background noise floor: pixels at or below 1 % of the dtype maximum are
    # considered true-dark (empty glass, zero-offset regions) and excluded from
    # the tissue background estimate.  Using the median of the remaining pixels
    # rather than the mean makes the estimate robust across different microscope
    # gain/offset settings (e.g. new CZI files vs. old LSM files).
    bg_noise_floor = dtype_max * 0.01

    all_bgs = []
    for z in range(nz):
        not_soma = ~combined_mask[z]
        if not_soma.sum() < 500:
            continue
        row = {'file': czi_fname, 'slice': z + 1}
        for ch in range(n_ch):
            if scale_pct is not None:
                plane = scale_plane_intensity(stack[ch, z], scale_pct, dtype_max)
            else:
                plane = stack[ch, z].astype(np.float32)
            bg_pix = plane[not_soma]
            tissue_pix = bg_pix[bg_pix > bg_noise_floor]
            if tissue_pix.size >= 100:
                bg_val = float(np.median(tissue_pix))
            else:
                bg_val = float(np.median(bg_pix))
            row['bg_ch{}'.format(ch + 1)] = round(bg_val, 3)
        all_bgs.append(row)

    hdr   = Font(name='Arial', bold=True, color='FFFFFF', size=11)
    sub   = Font(name='Arial', bold=True, size=10)
    dat   = Font(name='Arial', size=10)
    fills = {
        'meta': PatternFill('solid', start_color='4A4A4A'),
        'ch1':  PatternFill('solid', start_color='1F7D1F'),
        'ch2':  PatternFill('solid', start_color='B22222'),
        'ch3':  PatternFill('solid', start_color='1A5276'),
        'sub':  PatternFill('solid', start_color='D9D9D9'),
    }
    ctr  = Alignment(horizontal='center', vertical='center', wrap_text=True)
    thin = Side(style='thin',   color='AAAAAA')
    med  = Side(style='medium', color='555555')
    brd  = Border(left=thin, right=thin, top=thin, bottom=thin)
    brdR = Border(left=thin, right=med,  top=thin, bottom=thin)

    # ── Column layout (name-driven so columns can be added without renumbering) ─
    meta_hdrs = ['File', 'AnimalID', 'Genotype', 'Channel detected',
                 'Source channel', 'Soma #', 'Area (um2)', 'Slice',
                 'Centroid X (µm)', 'Centroid Y (µm)']
    n_meta   = len(meta_hdrs)
    col_hdrs = list(meta_hdrs)
    for n in hdr_names:
        col_hdrs += ['Mean {}'.format(n), 'BG {}'.format(n),
                     'F/F0 {}'.format(n), 'Range {}'.format(n)]
    name_to_col    = {h: i for i, h in enumerate(col_hdrs, 1)}
    _GENO_COL      = name_to_col['Genotype']
    _ANIMAL_COL    = name_to_col['AnimalID']
    _centroid_cols = {name_to_col['Centroid X (µm)'],
                      name_to_col['Centroid Y (µm)']}
    _ff0_cols   = {name_to_col['F/F0 {}'.format(n)] for n in hdr_names}
    _range_cols = {name_to_col['Range {}'.format(n)] for n in hdr_names}
    # Medium-right border after each logical group boundary.
    _brdR_cols  = ({_GENO_COL, name_to_col['Slice'], max(_centroid_cols)}
                   | _range_cols)

    appending = os.path.isfile(out_path)
    if appending:
        # rich_text=True so previously-written genotype cells keep their
        # italic/superscript formatting across the load/save cycle.
        wb  = openpyxl.load_workbook(out_path, rich_text=True)
        ws  = wb['Cells']
        ws2 = wb['Backgrounds']
        # Remove any existing rows for this file so re-export replaces rather
        # than duplicates them.  Iterate in reverse to keep row indices stable.
        for sheet, hdr_rows in [(ws, 2), (ws2, 1)]:
            to_delete = [r for r in sheet.iter_rows(min_row=hdr_rows + 1)
                         if r[0].value == czi_fname]
            for row in reversed(to_delete):
                sheet.delete_rows(row[0].row)
        data_start = ws.max_row + 1
        bg_start   = ws2.max_row + 1
    else:
        wb  = openpyxl.Workbook()
        ws  = wb.active;  ws.title = 'Cells'
        ws2 = wb.create_sheet('Backgrounds')
        data_start = 3
        bg_start   = 2

        # Cell Info band over the meta columns (File \u2026 Centroid Y).
        ws.merge_cells('A1:{}1'.format(get_column_letter(n_meta)))
        ws['A1'].value = 'Cell Info'
        ws['A1'].font = hdr; ws['A1'].fill = fills['meta']
        ws['A1'].alignment = ctr
        # 4 columns per channel: Mean, BG, F/F0, Range (p99-p1 of raw slice)
        ch_fkeys = ['ch1', 'ch2', 'ch3']
        for i, n in enumerate(hdr_names):
            s = n_meta + 1 + i * 4
            span = '{}1:{}1'.format(get_column_letter(s), get_column_letter(s + 3))
            ws.merge_cells(span)
            c = ws[span.split(':')[0]]
            c.value = 'Ch{} - {}'.format(i + 1, n)
            if i < len(ch_names) and n != ch_names[i]:
                # Name varies between the files of this batch.
                c.value += '  (per file: see Channel detected)'
            c.font = hdr; c.fill = fills[ch_fkeys[i] if i < 3 else 'ch3']
            c.alignment = ctr

        for ci, lbl in enumerate(col_hdrs, 1):
            c = ws.cell(row=2, column=ci, value=lbl)
            c.font = sub; c.fill = fills['sub']
            c.alignment = ctr; c.border = brd

        meta_w = {'File': 34, 'AnimalID': 20, 'Genotype': 22,
                  'Channel detected': 16, 'Source channel': 14, 'Soma #': 8,
                  'Area (um2)': 11, 'Slice': 7,
                  'Centroid X (\u00b5m)': 9, 'Centroid Y (\u00b5m)': 9}
        widths = ([meta_w.get(h, 12) for h in meta_hdrs]
                  + [11, 11, 9, 10] * len(ch_names))
        for ci, w in enumerate(widths, 1):
            ws.column_dimensions[get_column_letter(ci)].width = w
        ws.row_dimensions[1].height = 22
        ws.row_dimensions[2].height = 34
        ws.freeze_panes = 'A3'

        bg_hdrs = ['File', 'Slice'] + ['BG {}'.format(n) for n in hdr_names]
        for ci, lbl in enumerate(bg_hdrs, 1):
            c = ws2.cell(row=1, column=ci, value=lbl)
            c.font = sub; c.fill = fills['sub']; c.alignment = ctr
        for ci, w in enumerate([34, 7, 18, 18, 18], 1):
            ws2.column_dimensions[get_column_letter(ci)].width = w
        ws2.freeze_panes = 'A2'

    _animal = '' if animal_id is None else str(animal_id)
    bg_lut = {(r['file'], r['slice']): r for r in all_bgs}
    for ri, row in enumerate(all_cells, start=data_start):
        bg   = bg_lut.get((row['file'], row['slice']), {})
        # Values keyed by header name; Genotype is written separately (rich text).
        vals = {'File':             row['file'],
                'AnimalID':         _animal,
                'Channel detected': row['ch_detected'],
                'Source channel':   row.get('source_ch', ''),
                'Soma #':           row['soma_id'],
                'Area (um2)':       row['area_um2'],
                'Slice':            row['slice'],
                'Centroid X (µm)':  row.get('cx_um', ''),
                'Centroid Y (µm)':  row.get('cy_um', '')}
        for ch in range(n_ch):
            n = hdr_names[ch]
            m   = row.get('mean_ch{}'.format(ch + 1),  '')
            b   = bg.get('bg_ch{}'.format(ch + 1),     '')
            ff  = m / b if (b and b != 0) else ''
            rng = row.get('range_ch{}'.format(ch + 1), '')
            vals['Mean {}'.format(n)] = m
            vals['BG {}'.format(n)]   = b
            vals['F/F0 {}'.format(n)] = ff
            vals['Range {}'.format(n)] = rng
        for ci, hname in enumerate(col_hdrs, 1):
            if ci == _GENO_COL:
                c = ws.cell(row=ri, column=ci,
                            value=genotype_to_richtext(genotype_runs))
            else:
                val = vals.get(hname, '')
                c = ws.cell(row=ri, column=ci, value=val if val != '' else None)
                if isinstance(val, float):
                    if   ci in _ff0_cols:       c.number_format = '0.00'
                    elif ci in _range_cols:     c.number_format = '0.0'
                    elif ci in _centroid_cols:  c.number_format = '0.00'
                    else:                       c.number_format = '0.000'
            c.font = dat
            c.alignment = Alignment(horizontal='center', vertical='center')
            c.border = brdR if ci in _brdR_cols else brd

    for ri, row in enumerate(all_bgs, start=bg_start):
        ws2.append([row['file'], row['slice']] +
                   [row.get('bg_ch{}'.format(ch + 1), '') for ch in range(n_ch)])

    wb.save(out_path)
    return len(all_cells), len(all_bgs)


# ── Refiner CSV export ─────────────────────────────────────────────────────────

def export_refiner_data(stack, px_um, czi_fname, soma_by_ch, rejected_ids,
                        sens_by_ch, ch_names, out_path):
    n_ch, nz = stack.shape[:2]

    bg_by_z = {}
    for z in range(nz):
        row = {}
        for ch in range(n_ch):
            row[ch] = float(stack[ch, z].mean())
        bg_by_z[z] = row

    fieldnames = [
        'file', 'ch_detected', 'soma_id', 'source', 'shape_edited', 'status', 'slice',
        'area_um2',
        'area_pre_edit_um2', 'solidity_pre_edit',
        'compactness_pre_edit', 'minor_axis_pre_edit_um',
        'mean_ch1', 'mean_ch2', 'mean_ch3',
        'bg_ch1',   'bg_ch2',   'bg_ch3',
        'ff0_ch1',  'ff0_ch2',  'ff0_ch3',
        'compactness', 'solidity', 'eccentricity',
        'minor_axis_um', 'major_axis_um', 'inscribed_r_um',
        'peak_intensity', 'mean_intensity', 'bg_at_soma', 'snr_raw',
        'hp_mad', 'hp_threshold', 'hp_val_at_center', 'sigma_to_detect',
        'z_rel',
    ]

    rows   = []
    soma_n = 1
    for ch_idx, soma_list in soma_by_ch.items():
        sigma_mult = sens_by_ch.get(ch_idx, 2.5)
        for si, s in enumerate(soma_list):
            z    = s['z']
            mask = s['mask2d']
            src  = s.get('source', 'auto')
            stat = 'rejected' if (ch_idx, si) in rejected_ids else 'accepted'

            means = [round(float(stack[ch, z][mask].mean()), 3) for ch in range(n_ch)]
            bgs   = [bg_by_z[z][ch] for ch in range(n_ch)]
            ff0s  = [m / b if b > 0 else '' for m, b in zip(means, bgs)]

            diag  = compute_soma_diagnostics(s, stack, px_um, sigma_mult, ch_idx)

            row = {
                'file':                    czi_fname,
                'ch_detected':             ch_names[ch_idx],
                'soma_id':                 soma_n,
                'source':                  src,
                'shape_edited':            s.get('shape_edited', False),
                'status':                  stat,
                'slice':                   z + 1,
                'area_um2':                round(s['area_um2'], 3),
                'area_pre_edit_um2':       s.get('area_pre_edit_um2', ''),
                'solidity_pre_edit':       s.get('solidity_pre_edit', ''),
                'compactness_pre_edit':    s.get('compactness_pre_edit', ''),
                'minor_axis_pre_edit_um':  s.get('minor_axis_pre_edit_um', ''),
            }
            for i in range(n_ch):
                row['mean_ch{}'.format(i + 1)] = means[i]
                row['bg_ch{}'.format(i + 1)]   = round(bgs[i], 3)
                row['ff0_ch{}'.format(i + 1)]  = round(ff0s[i], 4) if ff0s[i] != '' else ''
            row.update(diag)
            rows.append(row)
            soma_n += 1

    append = False
    if os.path.isfile(out_path):
        try:
            with open(out_path, 'r', newline='', encoding='utf-8') as fcheck:
                existing_header = fcheck.readline().strip().split(',')
            append = (existing_header == fieldnames)
        except Exception:
            append = False

    mode = 'a' if append else 'w'
    with open(out_path, mode, newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not append:
            writer.writeheader()
        writer.writerows(rows)

    return len(rows)


# ── GUI ────────────────────────────────────────────────────────────────────────

class SomaDetectorGUI:
    def __init__(self, root, czi_path=None):
        self.root = root
        self.root.title('Soma Detector v3.0 — model: ' + cellpose_engine.model_name())
        self.root.protocol('WM_DELETE_WINDOW', self._on_close)

        self.czi_path     = None
        self.stack        = None
        self._stack_orig  = None   # raw array before channel reorder
        self._czi_ch_info = []     # [{name, fluor, ex_nm}, ...]
        self.px_um        = PIXEL_SIZE_UM_DEFAULT
        self.soma_by_ch   = {0: [], 1: []}
        self.rejected     = set()
        self.current_z  = -1
        self.current_ch = 0
        self._detecting = False
        self._patches   = []
        self._cp_detector = cellpose_engine.CellposeDetector()   # cached Cellpose inference / track pool
        self._cp_note     = ''

        self._drawing_roi = False
        self._roi_verts   = []
        self._roi_polygon = None
        self._roi_mask    = None

        self._edit_mode  = False
        self._edit_key   = None
        self._edit_verts = []
        self._drag_idx   = -1

        self._zoom_xlim = None
        self._zoom_ylim = None
        self._pan_last  = None   # (display_x, display_y) while middle-button held

        self._force_detect_mode  = False
        self._show_minmax_mask   = False
        self._intensity_scaling  = False
        self._scale_pct_var      = None  # initialised in _build_ui

        self._f_held        = False    # F key held → click force-detects
        self._f_clear_job   = None
        self._rerun_pending = False
        self._applied_order = None
        self._jump_target   = None     # persistent highlight set by click-to-jump
        self._build_ui()
        root.bind_all('<Escape>', lambda e: self._clear_jump_highlight())
        for seq in ('<KeyPress-f>', '<KeyPress-F>'):
            root.bind_all(seq, self._on_f_press)
        for seq in ('<KeyRelease-f>', '<KeyRelease-F>'):
            root.bind_all(seq, self._on_f_release)
        # If the window loses focus while F is held, its release never arrives
        # and every later click would force-detect: drop the hold on focus loss.
        root.bind('<FocusOut>', lambda e: self._clear_f() if self._f_held else None, add='+')

        if czi_path and os.path.isfile(czi_path):
            root.after(200, lambda: self._load_file(czi_path))   # window first, then the File-info dialog

        # ── Click-to-jump IPC server ──
        # Listens for jump commands from the Group Analysis GUI. Cached
        # search directories accumulate as files are resolved so subsequent
        # jumps within the same session are instant.
        self._jump_search_cache = []
        self._jump_server = cell_jump_ipc.JumpServer(
            self.root, on_message=self._handle_jump_command)
        ok, err = self._jump_server.start()
        if not ok:
            # Most common cause: port already bound by another soma GUI.
            # Print to stderr but keep running; the GUI is still usable.
            sys.stderr.write(
                'cell-jump IPC unavailable on port {}: {}\n'
                'The Group Analysis GUI will fall back to clipboard copy.\n'
                .format(cell_jump_ipc.DEFAULT_PORT, err))

    # ── UI ────────────────────────────────────────────────────────────────────

    HANDLE_GRAB_PX    = 10   # outline-edit: grab radius for a handle (screen px)
    EDGE_GRAB_PX      = 7    # outline-edit: click this close to the outline inserts a handle
    HANDLE_SPACING_UM = 3.0  # outline-edit: one handle per this much outline (~12 on a 12-µm soma)
    HANDLE_MIN, HANDLE_MAX = 8, 24

    def _build_ui(self):
        root = self.root

        # ── top bar ───────────────────────────────────────────────────────────
        top = tk.Frame(root, pady=4)
        top.pack(side=tk.TOP, fill=tk.X, padx=8)
        tk.Button(top, text='Open image', command=self._open_file,
                  width=11).pack(side=tk.LEFT)
        tk.Button(top, text='Load ROI Set', command=self._load_rois,
                  bg='#455A64', fg='white', font=('Arial', 9, 'bold'),
                  width=12).pack(side=tk.LEFT, padx=(6, 0))
        self.lbl_file = tk.Label(top, text='No file loaded', anchor='w',
                                 fg='#555', font=('Arial', 9))
        self.lbl_file.pack(side=tk.LEFT, padx=8)
        self.lbl_pxum = tk.Label(top, text='', anchor='w',
                                 fg='#888', font=('Arial', 8))
        self.lbl_pxum.pack(side=tk.LEFT, padx=4)

        # ── per-file metadata lives in the File-info dialog ───────────────────
        self.ch_name_var   = [tk.StringVar(value='GFP-Cre'), tk.StringVar(value='mCherry')]
        self.ch3_name_var  = tk.StringVar(value='pAntibody')
        self.animal_id_var = tk.StringVar(value='')
        self._meta_result  = tk.IntVar(value=0)
        self._build_meta_dialog()

        # ── file-info strip: what the File-info dialog confirmed ──────────────
        strip = tk.Frame(root, padx=8)
        strip.pack(side=tk.TOP, fill=tk.X, pady=(0, 3))
        tk.Button(strip, text='Edit file info…', font=('Arial', 8),
                  command=lambda: self._show_meta_dialog(forced=False)).pack(side=tk.LEFT)
        self.lbl_info_ch = tk.Label(strip, text='', font=('Arial', 8), fg='#333',
                                    anchor='w', justify='left', wraplength=720)
        self.lbl_info_ch.pack(side=tk.LEFT, padx=(8, 0))
        tk.Label(strip, text='Animal:', font=('Arial', 8, 'bold'),
                 fg='#333').pack(side=tk.LEFT, padx=(12, 2))
        self.lbl_info_animal = tk.Label(strip, text='—', font=('Arial', 8), fg='#333')
        self.lbl_info_animal.pack(side=tk.LEFT)
        tk.Label(strip, text='Genotype:', font=('Arial', 8, 'bold'),
                 fg='#333').pack(side=tk.LEFT, padx=(12, 2))
        self.txt_info_gt = tk.Text(strip, height=1, width=30, font=('Arial', 9),
                                   relief='flat', bg=root.cget('bg'),
                                   state='disabled', cursor='arrow')
        self.txt_info_gt.pack(side=tk.LEFT)
        self.txt_info_gt.tag_configure('gt_ital', font=('Arial', 9, 'italic'))
        self.txt_info_gt.tag_configure('gt_sup', font=('Arial', 7), offset=3)
        self.txt_info_gt.tag_configure('gt_italsup', font=('Arial', 7, 'italic'), offset=3)
        for v in (self.ch_name_var[0], self.ch_name_var[1], self.ch3_name_var, self.animal_id_var):
            v.trace_add('write', lambda *a: self._update_info_strip())

        # ── left control panel: fixed height, nothing scrolls ─────────────────
        left = tk.Frame(root, width=300, padx=8, pady=4)
        left.pack(side=tk.LEFT, fill=tk.Y)
        left.pack_propagate(False)

        def sep():
            ttk.Separator(left, orient='horizontal').pack(fill=tk.X, pady=5)

        def scale(parent, label, var, lo, hi, res):
            tk.Label(parent, text=label, font=('Arial', 8), fg='#444').pack(anchor='w')
            s = tk.Scale(parent, from_=lo, to=hi, resolution=res, orient=tk.HORIZONTAL,
                         variable=var, length=270, sliderlength=14, font=('Arial', 7))
            s.pack(anchor='w')
            s.bind('<ButtonRelease-1>', self._on_sel_release)
            return s
        self._make_scale = scale

        # View
        vf = tk.Frame(left); vf.pack(fill=tk.X)
        tk.Label(vf, text='View', font=('Arial', 9, 'bold')).pack(side=tk.LEFT)
        self.ch_var = tk.IntVar(value=0)
        self._ch_radios = []
        for val, txt in [(0, 'Ch1'), (1, 'Ch2'), (2, 'Both')]:
            rb = tk.Radiobutton(vf, text=txt, variable=self.ch_var, value=val,
                                command=self._redraw, font=('Arial', 8))
            rb.pack(side=tk.LEFT, padx=(6, 0))
            self._ch_radios.append(rb)

        def _update_radio(*_, i=0):
            self._ch_radios[i].config(
                text='Ch{} {}'.format(i + 1, self.ch_name_var[i].get())[:16])
        for i in (0, 1):
            self.ch_name_var[i].trace_add('write', lambda *a, i=i: _update_radio(i=i))
            _update_radio(i=i)
        sep()

        # Detect
        cp_ok = cellpose_engine.available()
        self.btn_detect = tk.Button(left, text='Detect somas   (Cellpose-SAM)',
                                    bg='#1976D2', fg='white', font=('Arial', 10, 'bold'),
                                    command=self._redetect,
                                    state=tk.NORMAL if cp_ok else tk.DISABLED)
        self.btn_detect.pack(fill=tk.X, pady=(0, 2))
        eng_txt = ('GPU: ' + cellpose_engine.device_name()) if cp_ok else \
                  'Cellpose not available in this Python — launch with Soma_Detector_3.0.bat'
        tk.Label(left, text=eng_txt, font=('Arial', 7), fg='#777',
                 wraplength=280, justify='left').pack(anchor='w')
        self.lbl_status = tk.Label(left, text='', font=('Arial', 8), fg='#333',
                                   wraplength=280, justify='left')
        self.lbl_status.pack(anchor='w')
        sep()

        # Selection — instant re-selection from the cached model output
        tk.Label(left, text='Which cells are shown   (re-selects instantly)',
                 font=('Arial', 9, 'bold')).pack(anchor='w')
        SD, MD = cellpose_engine.SELECT_DEFAULTS, cellpose_engine.MODEL_DEFAULTS
        self._rep_rule_labels = dict(cellpose_engine.REP_RULES)
        self.cp_vars = {
            'min_slices':      tk.IntVar(value=SD['min_slices']),
            'open_um':         tk.DoubleVar(value=SD['open_um']),
            'min_area':        tk.DoubleVar(value=SD['min_area_um2']),
            'max_area':        tk.DoubleVar(value=SD['max_area_um2']),
            'exclude_z_edges': tk.BooleanVar(value=SD['exclude_z_edges']),
            'rep_rule':        tk.StringVar(value=self._rep_rule_labels[SD['rep_rule']]),
            'cellprob':        tk.DoubleVar(value=MD['cellprob']),
            'diam':            tk.DoubleVar(value=MD['diam']),
        }
        self.cp_contrast_var = [tk.DoubleVar(value=SD['min_contrast']),
                                tk.DoubleVar(value=SD['min_contrast'])]
        _breadth = tk.DoubleVar(value=SD['breadth_um'])
        self.breadth_var = [_breadth, _breadth]          # one breadth for both channels
        pick_ok = cellpose_engine.pick_model() is not None
        self.cp_vars['min_pick'] = tk.DoubleVar(value=SD['min_pick'] if pick_ok else 0.0)
        scale(left, 'Selectivity   (ranking model: cells a person would circle)' if pick_ok
              else 'Selectivity   (no ranking model installed — inactive)',
              self.cp_vars['min_pick'], 0.0, 0.95, 0.05)
        scale(left, 'Min contrast Ch1   (more ←→ fewer cells)', self.cp_contrast_var[0], 1.0, 30.0, 0.5)
        scale(left, 'Min contrast Ch2', self.cp_contrast_var[1], 1.0, 30.0, 0.5)
        scale(left, 'Outline breadth µm   (tighter ←→ broader)', _breadth, -2.0, 2.0, 0.25)
        scale(left, 'Outline smoothing µm   (removes dendrite stubs)', self.cp_vars['open_um'], 0.0, 3.0, 0.3)
        tk.Button(left, text='Advanced settings…', font=('Arial', 8),
                  command=self._show_advanced).pack(anchor='w', pady=(2, 0))
        sep()

        # Curate
        tk.Label(left, text='Curate   (click an outline to reject / restore)',
                 font=('Arial', 9, 'bold')).pack(anchor='w')
        cf = tk.Frame(left); cf.pack(fill=tk.X)
        self.btn_force = tk.Button(cf, text='Force-detect: OFF', width=17,
                                   bg='#00838F', fg='white', font=('Arial', 8, 'bold'),
                                   command=self._toggle_force_detect)
        self.btn_force.pack(side=tk.LEFT)
        tk.Label(cf, text='or hold  F  and click', font=('Arial', 8),
                 fg='#00838F').pack(side=tk.LEFT, padx=6)
        self.lbl_force = tk.Label(left, text='', font=('Arial', 7), fg='#00838F',
                                  wraplength=280, justify='left')
        self.lbl_force.pack(anchor='w')
        ef = tk.Frame(left); ef.pack(fill=tk.X, pady=(2, 0))
        self.btn_edit_mode = tk.Button(ef, text='Edit outline: OFF', width=17,
                                       bg='#4A148C', fg='white', font=('Arial', 8, 'bold'),
                                       command=self._toggle_edit_mode)
        self.btn_edit_mode.pack(side=tk.LEFT)
        tk.Label(ef, text='drag handles · click edge adds · right-click removes',
                 font=('Arial', 7), fg='#7B1FA2', wraplength=150, justify='left').pack(side=tk.LEFT, padx=6)
        self.lbl_edit = tk.Label(left, text='', font=('Arial', 7), fg='#9C27B0',
                                 wraplength=280, justify='left')
        self.lbl_edit.pack(anchor='w')
        nf = tk.Frame(left); nf.pack(fill=tk.X)
        self.lbl_rejected = tk.Label(nf, text='Rejected: 0', font=('Arial', 8), fg='#c00')
        self.lbl_rejected.pack(side=tk.LEFT)
        self.lbl_manual = tk.Label(nf, text='Manual: 0', font=('Arial', 8), fg='#00838F')
        self.lbl_manual.pack(side=tk.LEFT, padx=12)
        sep()

        # Region
        rf = tk.Frame(left); rf.pack(fill=tk.X)
        tk.Label(rf, text='Region', font=('Arial', 9, 'bold')).pack(side=tk.LEFT)
        self.btn_draw_roi = tk.Button(rf, text='Draw ROI', width=9, bg='#E65100',
                                      fg='white', font=('Arial', 8, 'bold'),
                                      command=self._toggle_draw_roi)
        self.btn_draw_roi.pack(side=tk.LEFT, padx=(8, 4))
        tk.Button(rf, text='Clear', width=6, font=('Arial', 8),
                  command=self._clear_roi).pack(side=tk.LEFT)
        self.lbl_roi = tk.Label(left, text='No ROI set (whole image) — click vertices, double-click to close',
                                font=('Arial', 7), fg='#555', wraplength=280, justify='left')
        self.lbl_roi.pack(anchor='w')
        sep()

        # Display
        df = tk.Frame(left); df.pack(fill=tk.X)
        tk.Label(df, text='Display', font=('Arial', 9, 'bold')).pack(side=tk.LEFT)
        self._scale_pct_var = tk.DoubleVar(value=0.05)
        self.btn_scaling = tk.Button(df, text='Intensity Scaling: OFF', width=19,
                                     bg='#4E342E', fg='white', font=('Arial', 8, 'bold'),
                                     command=self._toggle_intensity_scaling)
        self.btn_scaling.pack(side=tk.LEFT, padx=(8, 4))
        tk.Label(df, text='clip %', font=('Arial', 7)).pack(side=tk.LEFT)
        pct_spin = tk.Spinbox(df, from_=0.01, to=25.0, increment=0.01,
                              textvariable=self._scale_pct_var, width=5,
                              font=('Arial', 7), format='%.2f',
                              command=self._on_scale_pct_change)
        pct_spin.pack(side=tk.LEFT, padx=2)
        pct_spin.bind('<Return>', lambda e: self._on_scale_pct_change())
        df2 = tk.Frame(left); df2.pack(fill=tk.X, pady=(2, 0))
        self.btn_minmax_mask = tk.Button(df2, text='Min/Max Mask: OFF', width=19,
                                         bg='#37474F', fg='white', font=('Arial', 8, 'bold'),
                                         command=self._toggle_minmax_mask)
        self.btn_minmax_mask.pack(side=tk.LEFT)
        self.lbl_pix_stats = tk.Label(df2, text='Zero: --   Sat: --',
                                      font=('Arial', 7), fg='#444')
        self.lbl_pix_stats.pack(side=tk.LEFT, padx=6)
        self.lbl_scaling = tk.Label(left, text='', font=('Arial', 7), fg='#BF360C',
                                    wraplength=280, justify='left')
        self.lbl_scaling.pack(anchor='w')
        sep()

        # Export
        self.btn_export = tk.Button(left, text='Export measurements',
                                    bg='#2E7D32', fg='white', font=('Arial', 10, 'bold'),
                                    command=self._export)
        self.btn_export.pack(fill=tk.X, pady=(0, 2))
        tk.Label(left, text='raw + scaled xlsx, refiner CSV and the *_ROIs.npz set',
                 font=('Arial', 7), fg='#777').pack(anchor='w')
        self.btn_batch = tk.Button(left, text='Batch process directory  (saved ROI sets)',
                                   bg='#4A148C', fg='white', font=('Arial', 9, 'bold'),
                                   command=self._batch_process)
        self.btn_batch.pack(fill=tk.X, pady=(6, 0))
        self.lbl_batch = tk.Label(left, text='', font=('Arial', 7), fg='#7B1FA2',
                                  wraplength=280, justify='left')
        self.lbl_batch.pack(anchor='w')

        # ── right panel ───────────────────────────────────────────────────────
        right = tk.Frame(root)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.fig, self.ax = plt.subplots(figsize=(7, 7))
        self.fig.patch.set_facecolor('#1a1a1a')
        self.ax.set_facecolor('#1a1a1a')
        self.ax.axis('off')
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect('button_press_event',   self._on_click)
        self.canvas.mpl_connect('motion_notify_event',  self._on_mouse_move)
        self.canvas.mpl_connect('button_release_event', self._on_release)
        self.canvas.mpl_connect('scroll_event',         self._on_scroll)

        bot = tk.Frame(right); bot.pack(fill=tk.X)
        tk.Label(bot, text='Slice:', font=('Arial', 10)).pack(side=tk.LEFT, padx=4)
        self.slice_var = tk.IntVar(value=0)
        self.sl_slice  = tk.Scale(bot, from_=0, to=1, orient=tk.HORIZONTAL,
                                  variable=self.slice_var, length=500,
                                  command=self._on_slice)
        self.sl_slice.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.lbl_slice = tk.Label(bot, text='Z-proj', font=('Arial', 10), width=6)
        self.lbl_slice.pack(side=tk.LEFT, padx=4)
        self._update_info_strip()

    # ── File-info dialog (channel assignment, animal ID, genotype) ────────────

    def _build_meta_dialog(self):
        """A persistent, hidden Toplevel: shown (modally) after every image
        load and from the File info… button.  Its widgets are the same objects
        the rest of the code reads (combos, animal_id_var, genotype_text)."""
        w = tk.Toplevel(self.root)
        w.withdraw()
        w.title('File info — channels, animal, genotype')
        w.resizable(False, False)
        w.protocol('WM_DELETE_WINDOW', self._meta_cancel)
        self._meta_win = w
        body = tk.Frame(w, padx=10, pady=8); body.pack(fill=tk.BOTH, expand=True)

        ch_frm = tk.LabelFrame(body, text='Channel Assignment',
                               font=('Arial', 9, 'bold'), padx=6, pady=4)
        ch_frm.pack(fill=tk.X)
        tk.Label(ch_frm,
                 text='Channel: emission and excitation λ  —  set the name to define virus and antibody',
                 font=('Arial', 8), fg='#555', justify='left').pack(anchor='w', pady=(0, 3))
        self._ch_assign_vars = []
        self._ch_assign_cbs  = []
        for lbl, name_var in [('Ch1 – detect (green):', self.ch_name_var[0]),
                              ('Ch2 – detect (red):',   self.ch_name_var[1]),
                              ('Ch3 – measure only:',   self.ch3_name_var)]:
            row = tk.Frame(ch_frm); row.pack(fill=tk.X, pady=1)
            tk.Label(row, text=lbl, font=('Arial', 8), width=19, anchor='w').pack(side=tk.LEFT)
            var = tk.StringVar(value='—')
            cb  = ttk.Combobox(row, textvariable=var, state='readonly', width=30, font=('Arial', 8))
            cb.pack(side=tk.LEFT)
            tk.Label(row, text='  name:', font=('Arial', 8), fg='#555').pack(side=tk.LEFT)
            tk.Entry(row, textvariable=name_var, font=('Arial', 8), width=14).pack(side=tk.LEFT)
            self._ch_assign_vars.append(var)
            self._ch_assign_cbs.append(cb)

        aid_frm = tk.LabelFrame(body, text='Animal ID   (unique per mouse, e.g. CC89603_111720#1 — guessed from the filename)',
                                font=('Arial', 9, 'bold'), padx=6, pady=4)
        aid_frm.pack(fill=tk.X, pady=(6, 0))
        aid_row = tk.Frame(aid_frm); aid_row.pack(fill=tk.X)
        tk.Entry(aid_row, textvariable=self.animal_id_var,
                 font=('Arial', 10), width=30).pack(side=tk.LEFT)
        tk.Button(aid_row, text='Re-guess', font=('Arial', 8),
                  command=self._fill_animal_guess).pack(side=tk.LEFT, padx=(6, 0))

        gt_frm = tk.LabelFrame(body, text='Mouse Genotype   (guessed from the filename — select text, then Italic / x²)',
                               font=('Arial', 9, 'bold'), padx=6, pady=4)
        gt_frm.pack(fill=tk.X, pady=(6, 0))
        gt_btns = tk.Frame(gt_frm); gt_btns.pack(anchor='w', pady=(0, 2))
        tk.Button(gt_btns, text='Italic', font=('Arial', 8, 'italic'), width=6,
                  command=lambda: self._toggle_gt_style('ital')).pack(side=tk.LEFT)
        tk.Button(gt_btns, text='x²', font=('Arial', 8, 'bold'), width=4,
                  command=lambda: self._toggle_gt_style('sup')).pack(side=tk.LEFT, padx=3)
        tk.Button(gt_btns, text='Re-guess', font=('Arial', 8), width=8,
                  command=self._fill_genotype_guess).pack(side=tk.LEFT, padx=(3, 0))
        self.genotype_text = tk.Text(gt_frm, height=2, width=60, font=('Arial', 10),
                                     wrap='word', relief='solid', borderwidth=1)
        self.genotype_text.pack(fill=tk.X)
        self.genotype_text.tag_configure('gt_ital', font=('Arial', 10, 'italic'))
        self.genotype_text.tag_configure('gt_sup', font=('Arial', 8), offset=4)
        self.genotype_text.tag_configure('gt_italsup', font=('Arial', 8, 'italic'), offset=4)

        bf = tk.Frame(body); bf.pack(fill=tk.X, pady=(8, 0))
        tk.Label(bf, text='Values carry over to the next image; the filename guess fills in where it can.',
                 font=('Arial', 7), fg='#777').pack(side=tk.LEFT)
        tk.Button(bf, text='Cancel', width=8, command=self._meta_cancel).pack(side=tk.RIGHT)
        self._meta_ok_btn = tk.Button(bf, text='OK — detect', width=12, bg='#1565C0', fg='white',
                                      font=('Arial', 9, 'bold'), command=self._meta_ok)
        self._meta_ok_btn.pack(side=tk.RIGHT, padx=6)

    def _show_meta_dialog(self, forced=True):
        """Modal File-info dialog.  Returns True on OK.  From the button (not
        forced) OK re-applies the channel order only if it changed."""
        if self.stack is None:
            if not forced:
                messagebox.showinfo('No file', 'Open an image first.')
            return False
        w = self._meta_win
        self._meta_ok_btn.config(text='OK — detect' if forced else 'OK')
        w.transient(self.root)
        w.deiconify(); w.lift()
        w.update_idletasks()
        w.geometry('+{}+{}'.format(self.root.winfo_rootx() + 80, self.root.winfo_rooty() + 80))
        try:
            w.grab_set()
        except tk.TclError:
            pass
        self._meta_result.set(-1)
        self.root.wait_variable(self._meta_result)
        try:
            w.grab_release()
        except tk.TclError:
            pass
        w.withdraw()
        self._update_info_strip()
        ok = self._meta_result.get() == 1
        if ok and not forced:
            if self._current_order() != self._applied_order:
                self._apply_ch_assignment()
                self._redetect()
            else:
                self._redraw()
        return ok

    def _meta_ok(self):
        self._meta_result.set(1)

    def _meta_cancel(self):
        self._meta_result.set(0)

    def _current_order(self):
        try:
            return [int(v.get().split(':')[0]) for v in self._ch_assign_vars]
        except (ValueError, IndexError):
            return None

    def _update_info_strip(self):
        """Show the confirmed channel assignment/names, animal ID and the
        formatted genotype in the strip under the top bar."""
        if not hasattr(self, 'txt_info_gt'):
            return
        if self.stack is None:
            self.lbl_info_ch.config(text='No file loaded')
        else:
            names = [self.ch_name_var[0].get(), self.ch_name_var[1].get(), self.ch3_name_var.get()]
            parts = []
            for i, nm in enumerate(names):
                phys = self._ch_assign_vars[i].get() if i < len(self._ch_assign_vars) else '—'
                phys = phys.split(': ', 1)[1] if ': ' in phys else phys      # drop the index
                parts.append('Ch{} {} = {}'.format(i + 1, nm.strip() or '?', ' '.join(phys.split())))
            self.lbl_info_ch.config(text='   ·   '.join(parts))
        self.lbl_info_animal.config(text=self.animal_id_var.get().strip() or '—')
        t = self.txt_info_gt
        t.config(state='normal')
        t.delete('1.0', 'end')
        runs = self._read_genotype_runs() if self.stack is not None else []
        if runs:
            for text, ital, sup in runs:
                tag = 'gt_italsup' if (ital and sup) else ('gt_ital' if ital else ('gt_sup' if sup else None))
                if tag:
                    t.insert('end', text, tag)
                else:
                    t.insert('end', text)
        else:
            t.insert('end', '—')
        t.config(state='disabled')

    # ── Advanced settings window ───────────────────────────────────────────────

    def _show_advanced(self):
        w = getattr(self, '_adv_win', None)
        if w is not None and w.winfo_exists():
            w.deiconify(); w.lift()
            return
        w = tk.Toplevel(self.root)
        w.title('Advanced settings')
        w.resizable(False, False)
        self._adv_win = w
        body = tk.Frame(w, padx=10, pady=8); body.pack()
        tk.Label(body, text='Selection rules   (re-select instantly, curation kept)',
                 font=('Arial', 9, 'bold')).pack(anchor='w')
        self._make_scale(body, 'Min z-slices per soma', self.cp_vars['min_slices'], 1, 8, 1)
        self._make_scale(body, 'Min area µm²', self.cp_vars['min_area'], 10, 100, 5)
        self._make_scale(body, 'Max area µm²', self.cp_vars['max_area'], 200, 1000, 25)
        tk.Checkbutton(body, text='Exclude surface-cut somas (largest cross-section on the top/bottom plane)',
                       variable=self.cp_vars['exclude_z_edges'], font=('Arial', 8),
                       command=self._on_sel_release).pack(anchor='w', pady=(4, 0))
        tk.Label(body, text='Representative plane of each soma', font=('Arial', 8), fg='#444').pack(anchor='w', pady=(4, 0))
        cb = ttk.Combobox(body, textvariable=self.cp_vars['rep_rule'],
                          values=list(self._rep_rule_labels.values()),
                          state='readonly', width=40, font=('Arial', 8))
        cb.pack(anchor='w')
        cb.bind('<<ComboboxSelected>>', self._on_sel_release)
        ttk.Separator(body, orient='horizontal').pack(fill=tk.X, pady=8)
        tk.Label(body, text='⟳ Model settings   (re-run the model, ~30 s per channel; curation is reset)',
                 font=('Arial', 9, 'bold')).pack(anchor='w')
        self._make_scale(body, '⟳ Cell probability threshold   (lower = more cells)', self.cp_vars['cellprob'], -4.0, 4.0, 0.5)
        self._make_scale(body, '⟳ Cell diameter px   (0 = automatic)', self.cp_vars['diam'], 0, 60, 5)
        tk.Button(body, text='Close', width=8, command=w.withdraw).pack(anchor='e', pady=(8, 0))

    # ── F key: hold to force-detect ────────────────────────────────────────────

    def _on_f_press(self, event):
        if event.widget.winfo_class() in ('Entry', 'Text', 'TEntry', 'TCombobox', 'Spinbox'):
            return
        if self._f_clear_job is not None:
            self.root.after_cancel(self._f_clear_job)
            self._f_clear_job = None
        if not self._f_held:
            self._f_held = True
            self.canvas.get_tk_widget().config(cursor='crosshair')

    def _on_f_release(self, event):
        # Windows auto-repeat fires release/press pairs while the key is held;
        # only a release NOT followed by a press within 120 ms ends the hold.
        if self._f_clear_job is not None:
            self.root.after_cancel(self._f_clear_job)
        self._f_clear_job = self.root.after(120, self._clear_f)

    def _clear_f(self):
        self._f_held = False
        self._f_clear_job = None
        self.canvas.get_tk_widget().config(cursor='')

    # ── Channel assignment ────────────────────────────────────────────────────

    def _apply_czi_ch_assignment_on_load(self):
        """Called after loading a CZI: populate dropdowns, auto-assign if possible."""
        n_ch = self._stack_orig.shape[0]
        info = self._czi_ch_info[:n_ch]

        opts = [format_channel_identifier(i, c) for i, c in enumerate(info)]
        # pad to 3 if fewer channels
        while len(opts) < 3:
            opts.append('{}: (n/a)'.format(len(opts)))

        for cb in self._ch_assign_cbs:
            cb['values'] = opts

        order = auto_assign_channels(info) if info else None

        if order and max(order) < n_ch:
            self.stack = self._stack_orig[order]
            for var, idx in zip(self._ch_assign_vars, order):
                var.set(opts[idx] if idx < len(opts) else opts[0])
            return True   # auto-assigned
        else:
            # Default: sequential
            for i, var in enumerate(self._ch_assign_vars):
                var.set(opts[i] if i < len(opts) else opts[0])
            return False  # not auto-assigned

    def _apply_ch_assignment(self):
        """Apply the user-selected channel order and clear detections."""
        if self._stack_orig is None:
            return
        try:
            order = []
            for var in self._ch_assign_vars:
                idx = int(var.get().split(':')[0])
                order.append(idx)
        except (ValueError, IndexError):
            messagebox.showwarning('Invalid',
                'Could not read channel selection. Load a file first.')
            return
        if len(set(order)) != 3:
            messagebox.showwarning('Invalid',
                'Each role must be assigned a different physical channel.')
            return
        if max(order) >= self._stack_orig.shape[0]:
            messagebox.showwarning('Invalid',
                'Selected channel index out of range.')
            return

        self.stack      = self._stack_orig[order]
        self.soma_by_ch = {0: [], 1: []}
        self.rejected.clear()
        self._cp_detector.clear()
        self._cp_note = ''
        self.lbl_status.config(text='Channel assignment applied.')
        self._update_info_strip()
        self._redraw()

    # ── Mouse genotype (rich-text: italic + superscript) ───────────────────────
    _GT_STYLE_TAGS = ('gt_ital', 'gt_sup', 'gt_italsup')

    def _apply_gt_style(self, start, end, italic, sup):
        """Set exactly one style tag (or none) on a text range."""
        t = self.genotype_text
        for tag in self._GT_STYLE_TAGS:
            t.tag_remove(tag, start, end)
        if italic and sup:
            t.tag_add('gt_italsup', start, end)
        elif italic:
            t.tag_add('gt_ital', start, end)
        elif sup:
            t.tag_add('gt_sup', start, end)

    def _char_gt_style(self, idx):
        """Return (italic, superscript) booleans for the character at idx."""
        names = set(self.genotype_text.tag_names(idx))
        italic = 'gt_ital' in names or 'gt_italsup' in names
        sup    = 'gt_sup'  in names or 'gt_italsup' in names
        return italic, sup

    def _toggle_gt_style(self, which):
        """Toggle italic or superscript on the current text selection."""
        t = self.genotype_text
        try:
            start = t.index('sel.first')
            end   = t.index('sel.last')
        except tk.TclError:
            return   # nothing selected
        idx = start
        while t.compare(idx, '<', end):
            nxt = t.index('{}+1c'.format(idx))
            italic, sup = self._char_gt_style(idx)
            if which == 'ital':
                italic = not italic
            else:
                sup = not sup
            self._apply_gt_style(idx, nxt, italic, sup)
            idx = nxt

    def _set_genotype_runs(self, runs):
        """Replace the genotype editor contents with the given runs."""
        if not hasattr(self, 'genotype_text'):
            return
        t = self.genotype_text
        t.delete('1.0', 'end')
        for text, italic, sup in (runs or []):
            start = t.index('end-1c')
            t.insert('end', text)
            self._apply_gt_style(start, t.index('end-1c'), italic, sup)
        self._update_info_strip()

    def _fill_animal_guess(self):
        """Populate the Animal ID field from the current filename's guess."""
        if hasattr(self, 'animal_id_var'):
            self.animal_id_var.set(
                guess_animal_id(getattr(self, 'czi_path', '') or ''))

    def _fill_genotype_guess(self):
        """Populate the genotype editor from the current filename's guess."""
        self._set_genotype_runs(
            guess_genotype_runs(getattr(self, 'czi_path', '') or ''))

    def _read_genotype_runs(self):
        """Serialize the editor to a list of (text, italic, superscript) runs."""
        if not hasattr(self, 'genotype_text'):
            return []
        content = self.genotype_text.get('1.0', 'end-1c')
        runs, cur, cur_style = [], '', None
        for i, ch in enumerate(content):
            style = self._char_gt_style('1.0+{}c'.format(i))
            if cur and style != cur_style:
                runs.append((cur, cur_style[0], cur_style[1]))
                cur = ''
            cur_style = style
            cur += ch
        if cur:
            runs.append((cur, cur_style[0], cur_style[1]))
        return runs

    def genotype_plain(self):
        """Plain-text genotype (formatting stripped)."""
        return ''.join(txt for txt, _i, _s in self._read_genotype_runs())

    # ── File loading ──────────────────────────────────────────────────────────

    def _open_file(self):
        path = filedialog.askopenfilename(
            title='Select .czi or .lsm file',
            filetypes=[('All files', '*.*'),
                       ('Zeiss images', '*.czi *.lsm'),
                       ('Zeiss CZI', '*.czi'),
                       ('Zeiss LSM', '*.lsm')])
        if path:
            self._load_file(path)

    def _load_file(self, path, prompt=True):
        """Load an image.  prompt=True (interactive) shows the File-info dialog
        and starts detection on OK; the click-to-jump path passes prompt=False
        because its metadata comes from the saved ROI set."""
        self.lbl_status.config(text='Loading...')
        self.root.update()
        try:
            ext = os.path.splitext(path)[1].lower()
            if ext == '.lsm':
                self._stack_orig, self.px_um, self._czi_ch_info = load_lsm(path)
                self.stack = self._stack_orig.copy()
            else:
                with czifile.CziFile(path) as czi:
                    self._stack_orig  = czi.asarray()
                    self.px_um        = get_pixel_size(czi.metadata())
                    self._czi_ch_info = read_czi_channel_info(czi.metadata())
                self.stack = self._stack_orig.copy()
            self.czi_path = path
            # Filename guesses; the previous image's values stay when a guess is empty.
            prev_animal = self.animal_id_var.get()
            prev_runs   = self._read_genotype_runs()
            self._fill_animal_guess()
            if not self.animal_id_var.get().strip():
                self.animal_id_var.set(prev_animal)
            self._fill_genotype_guess()
            if not self.genotype_plain().strip() and prev_runs:
                self._set_genotype_runs(prev_runs)
            nz = self.stack.shape[1]
            self.sl_slice.config(from_=0, to=nz)
            self.slice_var.set(0)
            self.current_z  = -1
            self.soma_by_ch = {0: [], 1: []}
            self.rejected.clear()
            self._cp_detector.clear()
            self._cp_note = ''
            self._zoom_xlim = self._zoom_ylim = None
            self._jump_target = None
            self.lbl_file.config(text=os.path.basename(path))
            dtype = self.stack.dtype
            if np.issubdtype(dtype, np.integer):
                bits   = np.iinfo(dtype).bits
                d_max  = np.iinfo(dtype).max
                depth_str = '  |  {}-bit (max {})'.format(bits, d_max)
            else:
                depth_str = '  |  float'
            self.lbl_pxum.config(
                text='  |  {:.4f} µm/px{}'.format(self.px_um, depth_str))
            if not (0.1 <= self.px_um <= 1.2):
                messagebox.showwarning(
                    'Check pixel size',
                    'This file reports {:.3f} µm/px — unusual for 40× soma imaging '
                    '(typically 0.3–0.7 µm/px).\n\nThe objective recorded in the '
                    'acquisition metadata may be wrong (e.g. a 10× overview saved '
                    'under a 40× name). Areas and detection would be off by the '
                    'square of the error.'.format(self.px_um))
            self._apply_czi_ch_assignment_on_load()      # fills the dialog's channel combos
            self._applied_order = self._current_order()
            self._redraw()
            if prompt:
                if self._show_meta_dialog(forced=True):
                    self._apply_ch_assignment()
                    self._redetect()
                else:
                    self.lbl_status.config(
                        text='File info not confirmed — use File info…, then Detect.')
            else:
                self.lbl_status.config(text='Loaded.')
            self._update_info_strip()
        except Exception as e:
            messagebox.showerror('Load error', str(e))
            self.lbl_status.config(text='Error loading file.')

    # ── Detection (Cellpose-SAM) ──────────────────────────────────────────────

    def _cp_params(self, ch):
        """(model params, selection params) for one channel from the controls."""
        model = {'cellprob': self.cp_vars['cellprob'].get(),
                 'flow':     cellpose_engine.MODEL_DEFAULTS['flow'],
                 'diam':     self.cp_vars['diam'].get()}
        sel = dict(cellpose_engine.SELECT_DEFAULTS)
        sel.update({'min_contrast':    self.cp_contrast_var[ch].get(),
                    'min_slices':      self.cp_vars['min_slices'].get(),
                    'open_um':         self.cp_vars['open_um'].get(),
                    'breadth_um':      self.breadth_var[ch].get(),
                    'min_area_um2':    self.cp_vars['min_area'].get(),
                    'max_area_um2':    self.cp_vars['max_area'].get(),
                    'exclude_z_edges': bool(self.cp_vars['exclude_z_edges'].get()),
                    'rep_rule':        self._rep_rule_key(),
                    'min_pick':        self.cp_vars['min_pick'].get()})
        return model, sel

    def _rep_rule_key(self):
        label = self.cp_vars['rep_rule'].get()
        for key, lbl in self._rep_rule_labels.items():
            if lbl == label:
                return key
        return cellpose_engine.SELECT_DEFAULTS['rep_rule']

    def _on_sel_release(self, _event=None):
        """A selection control changed: re-select from the cached model output
        (instant), or — for the ⟳ model settings — re-run the model.  Runs
        only once something has been detected."""
        if self.stack is None or not self._cp_detector._masks:
            return
        self._redetect()

    def _redetect(self, ch_idx=None):
        if self.stack is None:
            messagebox.showinfo('No file', 'Open an image first.')
            return
        if not cellpose_engine.available():
            messagebox.showinfo(
                'Cellpose unavailable',
                'Cellpose / PyTorch are not importable in this Python.\n'
                'Launch the detector with Soma_Detector_3.0.bat (python312 or .venv-cellpose).\n\n'
                + cellpose_engine.import_error())
            return
        if self._detecting:
            self._rerun_pending = True       # settings changed mid-run: run again when done
            return
        self._detecting = True
        self._rerun_pending = False
        if self._force_detect_mode:
            self._toggle_force_detect()
        self.btn_detect.config(state=tk.DISABLED)
        self.lbl_status.config(text='Detecting…')
        self.root.update()

        stack, px_um, roi_mask = self.stack, self.px_um, self._roi_mask
        cp_params = {ch: self._cp_params(ch) for ch in (0, 1) if ch < stack.shape[0]}
        det = self._cp_detector

        def _status(msg):
            self.root.after(0, lambda m=msg: self.lbl_status.config(text=m))

        def worker():
            try:
                result, notes = {}, []
                for ch, (model_p, sel_p) in cp_params.items():
                    somas, info = det.detect(stack[ch], px_um, model_p, sel_p,
                                             roi_mask=roi_mask, ch_key=ch, status=_status)
                    result[ch] = somas
                    notes.append('Ch{}: {} of {} cells pass{}'.format(
                        ch + 1, info['n_selected'], info['n_tracks'],
                        '' if info['cached'] else '  ({:.0f} s model)'.format(info['seconds_model'])))
                self._cp_note = '   '.join(notes)
                self.root.after(0, lambda: self._detection_done(result))
            except Exception as e:
                self.root.after(0, lambda: self._detection_error(str(e)))

        threading.Thread(target=worker, daemon=True).start()

    @staticmethod
    def _soma_key(s):
        """Identity of a soma across re-selections of the same model output."""
        return (int(round(s['cx'])), int(round(s['cy'])),
                int(s.get('z_min', s['z'])), int(s.get('z_max', s['z'])))

    def _detection_done(self, new_soma, ch_idx=None):
        """Install a (re-)selection, carrying curation over by soma identity:
        rejections, manual somas and hand-edited outlines survive slider
        changes (a ⟳ model re-run changes identities and starts clean)."""
        for ch, new_list in new_soma.items():
            old = self.soma_by_ch.get(ch, [])
            old_by_key = {self._soma_key(s): (s, (ch, i) in self.rejected)
                          for i, s in enumerate(old)}
            merged, rej_keys = [], set()
            for s in new_list:
                s['source'] = 'auto'
                k = self._soma_key(s)
                if k in old_by_key:
                    prev, was_rej = old_by_key[k]
                    if prev.get('shape_edited'):
                        s = prev
                    if was_rej:
                        rej_keys.add(k)
                merged.append(s)
            for s in old:
                if s.get('source') == 'manual':
                    merged.append(s)
                    if old_by_key[self._soma_key(s)][1]:
                        rej_keys.add(self._soma_key(s))
            self.soma_by_ch[ch] = merged
            self.rejected = {k for k in self.rejected if k[0] != ch} | \
                            {(ch, i) for i, s in enumerate(merged) if self._soma_key(s) in rej_keys}

        self._update_counts()
        self._detecting = False
        self.btn_detect.config(state=tk.NORMAL)
        self._redraw()
        if self._rerun_pending:
            self._rerun_pending = False
            self.root.after(50, self._redetect)

    def _detection_error(self, msg, ch_idx=None):
        self._detecting = False
        self._rerun_pending = False
        self.btn_detect.config(state=tk.NORMAL)
        messagebox.showerror('Detection error', msg)

    def _update_counts(self):
        n0  = len(self.soma_by_ch.get(0, []))
        n1  = len(self.soma_by_ch.get(1, []))
        nm  = sum(1 for sl in self.soma_by_ch.values()
                  for s in sl if s.get('source') == 'manual')
        rej = len(self.rejected)
        self.lbl_status.config(
            text='Ch1: {}  Ch2: {}  somas{}'.format(
                n0, n1, ('\n' + self._cp_note) if self._cp_note else ''))
        self.lbl_rejected.config(text='Rejected: {}'.format(rej))
        self.lbl_manual.config(text='Manual: {}'.format(nm))

    # ── Force-detect ──────────────────────────────────────────────────────────

    def _toggle_force_detect(self):
        self._force_detect_mode = not self._force_detect_mode
        if self._force_detect_mode:
            if self._edit_mode:
                self._toggle_edit_mode()
            if self._drawing_roi:
                self._toggle_draw_roi()
            self.btn_force.config(text='Force-detect: ON', bg='#006064')
            self.lbl_force.config(
                text='Click on a missed soma to force-segment it.\n'
                     'Go to the slice where it is most visible first.')
        else:
            self.btn_force.config(text='Force-detect: OFF', bg='#00838F')
            self.lbl_force.config(text='')

    def _force_detect_at(self, ex, ey):
        """Manual soma under the click: from the Cellpose pool on the CURRENT
        plane (a cell the selection rules excluded, or a cell wanted on this
        plane rather than the model's), otherwise by region growing."""
        ch_idx = min(self.current_ch, 1)
        nz     = self.stack.shape[1]
        if self.current_z < 0:
            col = int(round(np.clip(ex, 0, self.stack.shape[3] - 1)))
            row = int(round(np.clip(ey, 0, self.stack.shape[2] - 1)))
            z   = int(np.argmax(self.stack[ch_idx, :, row, col]))
        else:
            z = self.current_z
        z   = int(np.clip(z, 0, nz - 1))
        img = self.stack[ch_idx, z].astype(np.float32)
        ny, nx = img.shape
        px = int(round(np.clip(ex, 0, nx - 1)))
        py = int(round(np.clip(ey, 0, ny - 1)))
        cand = self._cp_detector.lookup(ch_idx, z, py, px, ny, nx)
        if cand is not None:
            return ch_idx, cand
        return ch_idx, self._grow_from_point(img, px, py, z)

    def _grow_from_point(self, img, px, py, z):
        smooth_px = max(0.5, SMOOTH_SIGMA_UM / self.px_um)
        bg_px     = max(2.0, BG_SIGMA_UM    / self.px_um)
        smoothed  = ndimage.gaussian_filter(img, sigma=smooth_px)
        bg        = ndimage.gaussian_filter(smoothed, sigma=bg_px)
        hp        = smoothed - bg
        hp      -= hp.min()

        click_val = float(hp[py, px])
        if click_val <= 0:
            return None

        bw      = hp > (click_val * 0.20)
        labeled, _ = ndimage.label(bw)
        lbl     = int(labeled[py, px])
        if lbl == 0:
            return None

        mask     = labeled == lbl
        expand_px = max(1, int(round(0.5 / self.px_um)))
        mask      = morphology.binary_dilation(mask, morphology.disk(expand_px))
        mask      = ndimage.binary_fill_holes(mask)

        rp_list  = measure.regionprops(mask.astype(np.uint8), intensity_image=img)
        if not rp_list:
            return None
        rp = rp_list[0]

        area_um2 = float(rp.area) * self.px_um ** 2
        min_area = max(15.0, self.cp_vars['min_area'].get() * 0.5)
        max_area = self.cp_vars['max_area'].get() * 2.0
        if area_um2 < min_area or area_um2 > max_area:
            self.lbl_force.config(
                text='Region found but area {:.0f} µm² out of range.'.format(area_um2))
            return None

        return {
            'z': z, 'z_min': z, 'z_max': z,
            'cy': float(rp.centroid[0]), 'cx': float(rp.centroid[1]),
            'area_um2':       area_um2,
            'cell_mean':      float(rp.mean_intensity),
            'peak_intensity': float(img[mask].max()),
            'mask2d':         mask,
            'source':         'manual',
        }

    def _n_manual(self):
        return sum(1 for sl in self.soma_by_ch.values()
                   for s in sl if s.get('source') == 'manual')

    # ── Drawing ───────────────────────────────────────────────────────────────

    def _on_slice(self, _=None):
        val = self.slice_var.get()
        if val == 0:
            self.current_z = -1
            self.lbl_slice.config(text='Z-proj')
        else:
            self.current_z = val - 1
            self.lbl_slice.config(text='z={}'.format(val))
        self._redraw()

    def _redraw(self, preview_pt=None):
        if self.stack is None:
            return
        self.current_ch = self.ch_var.get()
        z   = self.current_z
        sel = self.current_ch

        self.ax.cla()
        self.ax.axis('off')
        self._patches = []

        display_ch = 2 if sel == 2 else sel
        if z == -1:
            raw = self.stack[display_ch].max(axis=0).astype(np.float32)
        else:
            raw = self.stack[display_ch, z].astype(np.float32)

        # Apply per-plane display scaling if active
        disp = self._apply_display_scaling(raw)
        dmax = self._dtype_max()

        # When scaling is ON: show stats/mask on the scaled image so the user
        # can confirm each plane has been independently normalised (top X% will
        # appear as hot pink, bottom X% as dark blue).
        # When scaling is OFF: operate on raw data to show true saturation.
        stat_img = disp if self._intensity_scaling else raw
        pct_min, pct_max = self._compute_minmax_stats(stat_img)
        self.lbl_pix_stats.config(
            text='Zero px: {:.2f}%   Saturated: {:.2f}%'.format(pct_min, pct_max))

        # When scaling is ON use the full native range so the rescaling is
        # shown as-is; otherwise use the 99.5th-percentile autoscale.
        if self._intensity_scaling:
            vmax = dmax
        else:
            vmax = float(np.percentile(disp, 99.5))
            if vmax < 1:
                vmax = 1.0
        self.ax.imshow(disp, cmap='gray', vmin=0, vmax=vmax,
                       interpolation='nearest', aspect='equal')

        # Min/max mask overlay
        if self._show_minmax_mask:
            overlay = np.zeros((raw.shape[0], raw.shape[1], 4), dtype=np.float32)
            overlay[stat_img == 0]    = [0.00, 0.00, 0.65, 0.85]  # dark blue
            overlay[stat_img >= dmax] = [1.00, 0.07, 0.57, 0.85]  # hot pink
            self.ax.imshow(overlay, interpolation='nearest', aspect='equal')

        # Stats text pinned to bottom-left of the image
        scale_tag = '  [scaled {:.2f}%]'.format(
            self._scale_pct_var.get()) if self._intensity_scaling else ''
        self.ax.text(
            0.01, 0.01,
            'Zero: {:.2f}%   Sat: {:.2f}%{}'.format(pct_min, pct_max, scale_tag),
            transform=self.ax.transAxes, fontsize=8, color='white',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.55),
            verticalalignment='bottom')
        if self._zoom_xlim is not None:
            self.ax.set_xlim(self._zoom_xlim)
            self.ax.set_ylim(self._zoom_ylim)

        if self._roi_polygon is not None:
            self.ax.add_patch(plt.Polygon(
                self._roi_polygon, closed=True,
                facecolor='#FF6600', alpha=0.12,
                edgecolor='#FF6600', linewidth=2.0))
        elif self._drawing_roi and len(self._roi_verts) > 0:
            pts = list(self._roi_verts)
            if preview_pt and preview_pt[0] is not None:
                pts_draw = pts + [preview_pt]
            else:
                pts_draw = pts
            self.ax.plot([p[0] for p in pts_draw], [p[1] for p in pts_draw],
                         color='#FF6600', lw=1.5, alpha=0.9)
            self.ax.scatter([p[0] for p in pts], [p[1] for p in pts],
                            color='#FF6600', s=30, zorder=5)

        draw_channels = [(0, '#00FF00'), (1, '#FF4444')] if sel == 2 \
                        else [(sel, '#00FF00' if sel == 0 else '#FF4444')]

        n_accepted = 0
        for ch, auto_col in draw_channels:
            for si, s in enumerate(self.soma_by_ch.get(ch, [])):
                if z == -1:
                    show    = True
                    is_best = True
                else:
                    z_min = s.get('z_min', s['z'])
                    z_max = s.get('z_max', s['z'])
                    show    = (z_min <= z <= z_max)
                    is_best = (z == s['z'])
                if not show:
                    continue

                rejected = (ch, si) in self.rejected
                src      = s.get('source', 'auto')

                if rejected:
                    col = '#888888'
                elif src == 'manual':
                    col = MANUAL_COLOUR
                else:
                    col = auto_col

                alpha = (0.45 if rejected else 0.85) if is_best else \
                        (0.25 if rejected else 0.45)
                lw    = 1.8 if is_best else 1.0

                cnts = self._soma_contours(s)
                for cnt in cnts:
                    patch = MplPolygon(cnt[:, ::-1], closed=True,
                                       fill=False, edgecolor=col,
                                       linewidth=lw, alpha=alpha, picker=5)
                    patch._soma_key = (ch, si)
                    self.ax.add_patch(patch)
                    self._patches.append(patch)
                if not rejected:
                    n_accepted += 1

        if self._jump_target is not None:
            self._draw_jump_highlight(z)

        if self._edit_mode and self._edit_key is not None and self._edit_verts:
            xs = [v[0] for v in self._edit_verts]
            ys = [v[1] for v in self._edit_verts]
            self.ax.plot(xs + [xs[0]], ys + [ys[0]],
                         color='#CE93D8', lw=1.5, alpha=0.9, zorder=6)
            self.ax.scatter(xs, ys, color='#CE93D8', s=42,
                            zorder=7, edgecolors='white', linewidths=0.8)

        ch_labels = {0: 'Ch1 {}'.format(self.ch_name_var[0].get()),
                     1: 'Ch2 {}'.format(self.ch_name_var[1].get()),
                     2: 'Ch3 {}'.format(self.ch3_name_var.get())}

        mode_tag = ''
        if self._edit_mode:
            mode_tag = '  [EDIT]'
        elif self._force_detect_mode:
            mode_tag = '  [FORCE-DETECT — click missed soma]'
        if self._jump_target is not None:
            mode_tag += '  [JUMP → {}  (Esc clears)]'.format(self._jump_target['label'])

        title_col = '#CE93D8' if self._edit_mode else \
                    ('#00E5FF' if self._force_detect_mode else 'white')
        self.ax.set_title(
            '{}  |  {}  |  {} accepted{}'.format(
                'Z-proj' if z == -1 else 'z={}'.format(z + 1),
                ch_labels[sel], n_accepted, mode_tag),
            color=title_col, fontsize=10, pad=4)
        self.fig.patch.set_facecolor('#1a1a1a')
        self.canvas.draw_idle()

    # ── Click handler ─────────────────────────────────────────────────────────

    def _soma_contours(self, s):
        """Outline(s) of a soma mask in image coordinates, computed on the
        mask's bounding box and cached on the soma dict (a full-frame
        find_contours per soma per redraw is too slow with Cellpose's
        several hundred cells per channel)."""
        c = s.get('_contours')
        if c is None:
            m = s['mask2d']
            rr, cc = np.nonzero(m)
            if len(rr) == 0:
                c = []
            else:
                r0, r1 = max(0, rr.min() - 2), min(m.shape[0], rr.max() + 3)
                c0, c1 = max(0, cc.min() - 2), min(m.shape[1], cc.max() + 3)
                c = [cnt + np.array([r0, c0])
                     for cnt in measure.find_contours(
                         m[r0:r1, c0:c1].astype(np.float32), 0.5)]
            s['_contours'] = c
        return c

    def _on_click(self, event):
        if event.inaxes != self.ax or self.stack is None:
            return

        if event.button == 2:                       # middle button → pan
            self._pan_last = (event.x, event.y)
            return

        force = self._force_detect_mode or self._f_held
        if event.dblclick and not self._drawing_roi and not self._edit_mode and not force:
            self._reset_view()
            self._redraw()
            return

        if self._drawing_roi:
            if event.button != 1:
                return
            if event.dblclick:
                self._close_roi()
            else:
                self._roi_verts.append((event.xdata, event.ydata))
                self.lbl_roi.config(
                    text='Click to add vertices... ({} so far)'.format(
                        len(self._roi_verts)), fg='#b71c1c')
                self._redraw()
            return

        ex, ey = event.xdata, event.ydata

        if self._edit_mode:
            if self._edit_verts:
                disp = self.ax.transData.transform(np.array(self._edit_verts, float))
                d = np.hypot(disp[:, 0] - event.x, disp[:, 1] - event.y)
                i = int(np.argmin(d))
                if d[i] <= self.HANDLE_GRAB_PX:          # nearest handle, in screen pixels
                    if event.button == 3:                # right-click removes it
                        if len(self._edit_verts) > 4:
                            del self._edit_verts[i]
                            self._apply_edit_to_mask()
                        return
                    self._drag_idx = i
                    return
                if event.button == 1:                    # click on the outline inserts a handle
                    seg = self._nearest_segment(disp, event.x, event.y)
                    if seg is not None and seg[1] <= self.EDGE_GRAB_PX:
                        self._edit_verts.insert(seg[0] + 1, [ex, ey])
                        self._drag_idx = seg[0] + 1
                        self._redraw()
                        return
            if event.button != 1:
                return
            for patch in self._patches:
                contains, _ = patch.contains(event)
                if contains:
                    self._select_soma_for_edit(patch._soma_key)
                    return
            self._edit_key   = None
            self._edit_verts = []
            self._drag_idx   = -1
            self.lbl_edit.config(text='Click an outline to select it.')
            self._redraw()
            return

        if event.button != 1:
            return

        if force:
            self.lbl_force.config(text='Segmenting…')
            self.root.update()
            ch_idx, cand = self._force_detect_at(ex, ey)
            if cand is None:
                self.lbl_force.config(
                    text='No cell found. Try clicking closer to the soma centre.')
                return
            for si, s in enumerate(self.soma_by_ch.get(ch_idx, [])):
                if (s.get('source') == 'manual' and
                        np.hypot(s['cx'] - cand['cx'], s['cy'] - cand['cy']) < NMS_XY_RADIUS_UM / self.px_um):
                    self.lbl_force.config(text='Duplicate of an existing manual soma.')
                    return
            self.soma_by_ch[ch_idx].append(cand)
            nm = self._n_manual()
            self.lbl_force.config(
                text='Added manual soma #{} (Ch{}, z={}).'.format(
                    nm, ch_idx + 1, cand['z'] + 1))
            self._update_counts()
            self._redraw()
            return

        # Several outlines can contain the click (a force-detected soma drawn
        # over the rejected auto outline it replaced, or overlapping cells).
        # Prefer the one that belongs on THIS plane, then an accepted one, then
        # the smallest — never simply the first drawn, which is always the old
        # auto soma and made a force-detected soma impossible to toggle.
        z_now = self.current_z

        def _priority(key):
            ch, si = key
            s = self.soma_by_ch[ch][si]
            return (0 if (z_now < 0 or s['z'] == z_now) else 1,
                    1 if key in self.rejected else 0,
                    s.get('area_um2', 0.0))

        hits = {patch._soma_key for patch in self._patches
                if patch.contains(event)[0]
                and patch._soma_key[1] < len(self.soma_by_ch.get(patch._soma_key[0], []))}
        if hits:
            key = min(hits, key=_priority)
            if key in self.rejected:
                self.rejected.discard(key)
            else:
                self.rejected.add(key)
            self._update_counts()
            self._redraw()
            return

        best_key, best_rank = None, (30.0,)
        for patch in self._patches:
            ch, si = patch._soma_key
            soma_list = self.soma_by_ch.get(ch, [])
            if si >= len(soma_list):
                continue
            s = soma_list[si]
            d = np.hypot(s['cx'] - ex, s['cy'] - ey)
            rank = (round(d, 1),) + _priority((ch, si))
            if d < 30.0 and rank < best_rank:
                best_rank = rank
                best_key  = patch._soma_key
        if best_key is not None:
            if best_key in self.rejected:
                self.rejected.discard(best_key)
            else:
                self.rejected.add(best_key)
            self._update_counts()
            self._redraw()

    @staticmethod
    def _nearest_segment(disp, x, y):
        """(segment index, distance in screen px) of the closed-polygon edge
        nearest to the screen point (x, y)."""
        n = len(disp)
        if n < 2:
            return None
        p = np.array([x, y], float)
        best = None
        for j in range(n):
            a, b = disp[j], disp[(j + 1) % n]
            ab = b - a
            L2 = float(ab @ ab)
            t = 0.0 if L2 == 0 else float(np.clip((p - a) @ ab / L2, 0.0, 1.0))
            dist = float(np.hypot(*(a + t * ab - p)))
            if best is None or dist < best[1]:
                best = (j, dist)
        return best

    # ── Mouse move / release / scroll / view ───────────────────────────────────

    def _on_mouse_move(self, event):
        if self._pan_last is not None and event.x is not None:
            dx = event.x - self._pan_last[0]
            dy = event.y - self._pan_last[1]
            self._pan_last = (event.x, event.y)
            bbox = self.ax.get_window_extent()
            if bbox.width > 0 and bbox.height > 0:
                x0, x1 = self.ax.get_xlim()
                y0, y1 = self.ax.get_ylim()
                xs = (x1 - x0) / bbox.width
                ys = (y1 - y0) / bbox.height
                self._set_view(x0 - dx * xs, x1 - dx * xs, y0 - dy * ys, y1 - dy * ys)
            return

        if event.inaxes != self.ax:
            return
        if self._edit_mode and self._drag_idx >= 0 and event.xdata is not None:
            self._edit_verts[self._drag_idx] = [event.xdata, event.ydata]
            self._redraw()
            return
        if self._drawing_roi and len(self._roi_verts) > 0:
            self._redraw(preview_pt=(event.xdata, event.ydata))

    def _on_release(self, event):
        if event.button == 2:
            self._pan_last = None
            return
        if self._edit_mode and self._drag_idx >= 0:
            self._drag_idx = -1
            self._apply_edit_to_mask()

    def _set_view(self, x0, x1, y0, y1):
        """Apply a view in data coordinates (y0 > y1: images are drawn
        top-down), clamped inside the image.  A view as large as the image
        resets to the full, centred image."""
        ny, nx = self.stack.shape[2], self.stack.shape[3]
        w, h = abs(x1 - x0), abs(y0 - y1)
        if w >= nx - 1 or h >= ny - 1:
            self._reset_view()
            return
        lo_x, hi_x = min(x0, x1), max(x0, x1)
        if lo_x < -0.5:
            hi_x += -0.5 - lo_x; lo_x = -0.5
        if hi_x > nx - 0.5:
            lo_x -= hi_x - (nx - 0.5); hi_x = nx - 0.5
        lo_y, hi_y = min(y0, y1), max(y0, y1)
        if lo_y < -0.5:
            hi_y += -0.5 - lo_y; lo_y = -0.5
        if hi_y > ny - 0.5:
            lo_y -= hi_y - (ny - 0.5); hi_y = ny - 0.5
        self._zoom_xlim = [lo_x, hi_x]
        self._zoom_ylim = [hi_y, lo_y]
        self.ax.set_xlim(self._zoom_xlim)
        self.ax.set_ylim(self._zoom_ylim)
        self.canvas.draw_idle()

    def _reset_view(self):
        self._zoom_xlim = None
        self._zoom_ylim = None
        if self.stack is not None:
            ny, nx = self.stack.shape[2], self.stack.shape[3]
            self.ax.set_xlim(-0.5, nx - 0.5)
            self.ax.set_ylim(ny - 0.5, -0.5)
        self.canvas.draw_idle()

    def _on_scroll(self, event):
        """Zoom about the cursor; zooming out past the image re-centres it."""
        if event.inaxes != self.ax or event.xdata is None or self.stack is None:
            return
        f = 0.8 if event.step > 0 else 1.25
        x0, x1 = self.ax.get_xlim()
        y0, y1 = self.ax.get_ylim()
        xd, yd = event.xdata, event.ydata
        self._set_view(xd - (xd - x0) * f, xd + (x1 - xd) * f,
                       yd - (yd - y0) * f, yd + (y1 - yd) * f)

    # ── ROI polygon ───────────────────────────────────────────────────────────

    def _toggle_draw_roi(self):
        if self._drawing_roi:
            self._drawing_roi = False
            self._roi_verts   = []
            self.btn_draw_roi.config(text='Draw ROI', bg='#E65100')
            self._redraw()
        else:
            self._drawing_roi = True
            self._roi_verts   = []
            self._roi_polygon = None
            self._roi_mask    = None
            self.btn_draw_roi.config(text='Cancel ROI', bg='#b71c1c')
            self.lbl_roi.config(text='Click to add vertices...', fg='#b71c1c')
            self._redraw()

    def _clear_roi(self):
        self._drawing_roi = False
        self._roi_verts   = []
        self._roi_polygon = None
        self._roi_mask    = None
        self.btn_draw_roi.config(text='Draw ROI', bg='#E65100')
        self.lbl_roi.config(text='No ROI set (whole image)', fg='#555')
        self._redraw()

    def _close_roi(self):
        if len(self._roi_verts) < 3:
            self.lbl_roi.config(text='Need at least 3 points', fg='#c00')
            return
        verts = np.array(self._roi_verts)
        self._roi_polygon = verts
        if self.stack is not None:
            ny, nx = self.stack.shape[2], self.stack.shape[3]
            rr, cc = sk_polygon(verts[:, 1], verts[:, 0], shape=(ny, nx))
            mask   = np.zeros((ny, nx), bool)
            mask[rr, cc] = True
            self._roi_mask = mask
        self._drawing_roi = False
        self._roi_verts   = []
        self.btn_draw_roi.config(text='Draw ROI', bg='#E65100')
        self.lbl_roi.config(text='ROI set ({} vertices)'.format(len(verts)), fg='#1a7a1a')
        self._redraw()

    # ── Shape-edit helpers ────────────────────────────────────────────────────

    def _toggle_edit_mode(self):
        self._edit_mode = not self._edit_mode
        if self._edit_mode:
            if self._force_detect_mode:
                self._toggle_force_detect()
            self.btn_edit_mode.config(text='Edit outline: ON', bg='#6A1B9A')
            self.lbl_edit.config(text='Click an outline to select it, then drag handles.')
        else:
            self._edit_key   = None
            self._edit_verts = []
            self._drag_idx   = -1
            self.btn_edit_mode.config(text='Edit outline: OFF', bg='#4A148C')
            self.lbl_edit.config(text='')
        self._redraw()

    def _select_soma_for_edit(self, key):
        ch, si = key
        soma_list = self.soma_by_ch.get(ch, [])
        if si >= len(soma_list):
            return
        s    = soma_list[si]
        mask = s['mask2d']
        cnts = measure.find_contours(mask.astype(np.float32), 0.5)
        if not cnts:
            return

        # Snapshot pre-edit shape features on first selection only, so
        # re-selecting after a drag doesn't overwrite the original values.
        if 'area_pre_edit_um2' not in s:
            rp = measure.regionprops(mask.astype(np.uint8))
            if rp:
                rp = rp[0]
                hull_area = rp.convex_area if rp.convex_area > 0 else rp.area
                s['area_pre_edit_um2']        = round(s['area_um2'], 3)
                s['solidity_pre_edit']        = round(rp.area / hull_area, 4)
                s['compactness_pre_edit']     = round(
                    (4 * np.pi * rp.area) / (rp.perimeter ** 2)
                    if rp.perimeter > 0 else 0.0, 4)
                s['minor_axis_pre_edit_um']   = round(
                    rp.minor_axis_length * self.px_um, 3)

        # A handle every HANDLE_SPACING_UM of outline (HANDLE_MIN–HANDLE_MAX per
        # cell), evenly spaced; a click on the outline inserts more where needed.
        cnt   = max(cnts, key=len)
        seg   = np.hypot(np.diff(cnt[:, 0]), np.diff(cnt[:, 1]))
        perim = float(seg.sum())
        n_h   = int(np.clip(round(perim * self.px_um / self.HANDLE_SPACING_UM),
                            self.HANDLE_MIN, self.HANDLE_MAX))
        cum   = np.concatenate([[0.0], np.cumsum(seg)])
        idx   = np.searchsorted(cum, np.linspace(0, perim, n_h, endpoint=False))
        pts   = cnt[np.clip(idx, 0, len(cnt) - 1)]
        self._edit_key   = key
        self._edit_verts = [[pt[1], pt[0]] for pt in pts]
        self._drag_idx   = -1
        self.lbl_edit.config(
            text='Editing Ch{} soma {}  —  drag handles'.format(ch + 1, si + 1))
        self._redraw()

    def _apply_edit_to_mask(self):
        if self._edit_key is None or len(self._edit_verts) < 3:
            return
        ch, si = self._edit_key
        soma_list = self.soma_by_ch.get(ch, [])
        if si >= len(soma_list):
            return
        ny, nx = self.stack.shape[2], self.stack.shape[3]
        verts  = np.array(self._edit_verts)
        rr, cc = sk_polygon(verts[:, 1], verts[:, 0], shape=(ny, nx))
        new_mask      = np.zeros((ny, nx), bool)
        new_mask[rr, cc] = True
        s                  = soma_list[si]
        s['mask2d']        = new_mask
        s.pop('_contours', None)          # outline cache is stale
        s['area_um2']      = float(new_mask.sum()) * self.px_um ** 2
        s['shape_edited']  = True          # flag: boundary was manually adjusted
        rr2, cc2           = np.where(new_mask)
        if len(rr2):
            s['cy'] = float(rr2.mean())
            s['cx'] = float(cc2.mean())
        self._redraw()

    # ── Pixel statistics & scaling helpers ───────────────────────────────────

    def _dtype_max(self):
        if self.stack is None:
            return 65535.0
        if np.issubdtype(self.stack.dtype, np.integer):
            return float(np.iinfo(self.stack.dtype).max)
        return float(self.stack.max()) or 65535.0

    def _compute_minmax_stats(self, img):
        n = float(img.size)
        if n == 0:
            return 0.0, 0.0
        dmax    = self._dtype_max()
        pct_min = 100.0 * float((img == 0).sum()) / n
        pct_max = 100.0 * float((img >= dmax).sum()) / n
        return pct_min, pct_max

    def _apply_display_scaling(self, raw):
        if not self._intensity_scaling:
            return raw
        pct    = self._scale_pct_var.get()
        p_low  = float(np.percentile(raw, pct))
        p_high = float(np.percentile(raw, 100.0 - pct))
        if p_high <= p_low:
            return raw
        scaled = np.clip(raw, p_low, p_high)
        scaled = (scaled - p_low) / (p_high - p_low) * self._dtype_max()
        return scaled.astype(np.float32)

    def _toggle_minmax_mask(self):
        self._show_minmax_mask = not self._show_minmax_mask
        if self._show_minmax_mask:
            self.btn_minmax_mask.config(text='Min/Max Mask: ON',  bg='#1A237E')
        else:
            self.btn_minmax_mask.config(text='Min/Max Mask: OFF', bg='#37474F')
        self._redraw()

    def _toggle_intensity_scaling(self):
        self._intensity_scaling = not self._intensity_scaling
        if self._intensity_scaling:
            pct = self._scale_pct_var.get()
            self.btn_scaling.config(text='Intensity Scaling: ON',  bg='#BF360C')
            self.lbl_scaling.config(
                text='Clipping top & bottom {:.2f}% per plane/channel.'.format(pct))
        else:
            self.btn_scaling.config(text='Intensity Scaling: OFF', bg='#4E342E')
            self.lbl_scaling.config(text='')
        self._redraw()

    def _on_scale_pct_change(self):
        if self._intensity_scaling:
            pct = self._scale_pct_var.get()
            self.lbl_scaling.config(
                text='Clipping top & bottom {:.2f}% per plane/channel.'.format(pct))
            self._redraw()

    # ── Shape filters toggle ──────────────────────────────────────────────────

    # ── Export: standard + scaled measurements ───────────────────────────────

    def _export(self):
        if self.stack is None:
            messagebox.showinfo('No data', 'Load a file and run detection first.')
            return
        total = sum(len(v) for v in self.soma_by_ch.values())
        if total == 0:
            messagebox.showinfo('No detections', 'Run detection first.')
            return

        out_dir    = os.path.dirname(self.czi_path)
        _pfx       = os.path.basename(out_dir.rstrip(os.sep)).replace(' ', '')
        out_path   = os.path.join(out_dir, f'{_pfx}_soma_measurements.xlsx')
        out_scaled = os.path.join(out_dir, f'{_pfx}_intensity_scaled_soma_measurements.xlsx')
        ch_names   = [self.ch_name_var[0].get(), self.ch_name_var[1].get(),
                      self.ch3_name_var.get()]
        czi_fname  = os.path.basename(self.czi_path)
        appending  = os.path.isfile(out_path)
        scale_pct  = self._scale_pct_var.get()

        out_refiner = os.path.join(out_dir, f'{_pfx}_refiner_training_data.csv')
        sens_by_ch  = {0: 2.5, 1: 2.5}        # diagnostics sigma for the refiner CSV (classical engine removed)
        geno_runs   = self._read_genotype_runs()
        animal_id   = self.animal_id_var.get().strip()

        try:
            # Raw (unscaled) export
            n_cells, _n_bg = measure_and_export(
                self.stack, self.px_um, czi_fname,
                self.soma_by_ch, self.rejected,
                ch_names, out_path, genotype_runs=geno_runs,
                animal_id=animal_id)

            # Intensity-scaled export (always produced; uses current clip %)
            measure_and_export(
                self.stack, self.px_um, czi_fname,
                self.soma_by_ch, self.rejected,
                ch_names, out_scaled, scale_pct=scale_pct,
                genotype_runs=geno_runs, animal_id=animal_id)

            # Refiner training data (auto-exported every time)
            export_refiner_data(
                self.stack, self.px_um, czi_fname,
                self.soma_by_ch, self.rejected,
                sens_by_ch, ch_names, out_refiner)

            roi_path = self._save_rois(out_dir)
            action   = 'Appended to' if appending else 'Created'
            messagebox.showinfo(
                'Exported',
                '{a} {p}_soma_measurements.xlsx\n'
                '{a} {p}_intensity_scaled_soma_measurements.xlsx  (clip {clip:.2f}%)\n'
                '{a} {p}_refiner_training_data.csv\n'
                '{n} cells  ({r} rejected)\n'
                'ROI set: {roi}'.format(
                    a=action, p=_pfx, clip=scale_pct, n=n_cells, r=len(self.rejected),
                    roi=os.path.basename(roi_path) if roi_path else 'n/a'))
        except Exception as e:
            messagebox.showerror('Export error', str(e))

    # ── Batch processing ──────────────────────────────────────────────────────

    def _batch_process(self):
        directory = filedialog.askdirectory(
            title='Select directory containing images and *_ROIs.npz files')
        if not directory:
            return

        npz_files = sorted(glob.glob(os.path.join(directory, '*_ROIs.npz')))
        if not npz_files:
            messagebox.showinfo(
                'No NPZ files',
                'No *_ROIs.npz files found in:\n{}\n\n'
                'Open each image, run detection, curate, and\n'
                'click "Export measurements" to save ROI sets first.'.format(directory))
            return

        # Match each NPZ to its image file
        pairs = []
        missing = []
        for npz_path in npz_files:
            base = os.path.basename(npz_path)[:-len('_ROIs.npz')]
            found = False
            for ext in ['.czi', '.lsm']:
                img_path = os.path.join(directory, base + ext)
                if os.path.isfile(img_path):
                    pairs.append((img_path, npz_path))
                    found = True
                    break
            if not found:
                missing.append(base)

        if not pairs:
            messagebox.showinfo(
                'No matching images',
                'Found {} NPZ files but no matching .czi/.lsm '
                'images in the same directory.'.format(len(npz_files)))
            return

        # ── File-list confirmation dialog ──────────────────────────────────────
        confirmed = [False]
        dlg = tk.Toplevel(self.root)
        dlg.title('Batch Process Directory')
        dlg.resizable(True, True)
        dlg.grab_set()

        tk.Label(dlg, text='Directory:  {}'.format(directory),
                 font=('Arial', 9), anchor='w').pack(fill=tk.X, padx=10, pady=(10, 2))
        tk.Label(dlg,
                 text='Clip %: {:.2f}     Files to process: {}{}'.format(
                     self._scale_pct_var.get(), len(pairs),
                     '     ({}  NPZ without image skipped)'.format(len(missing)) if missing else ''),
                 font=('Arial', 9), anchor='w').pack(fill=tk.X, padx=10, pady=(0, 6))

        # Scrollable file list
        list_frame = tk.Frame(dlg)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=2)
        sb = tk.Scrollbar(list_frame, orient=tk.VERTICAL)
        lb = tk.Listbox(list_frame, yscrollcommand=sb.set,
                        font=('Courier', 9), selectmode=tk.BROWSE,
                        width=70, height=min(20, len(pairs) + len(missing) + 1))
        sb.config(command=lb.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        lb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        for img_path, npz_path in pairs:
            lb.insert(tk.END,
                      '\u2713  {}   \u2192  {}'.format(
                          os.path.basename(img_path),
                          os.path.basename(npz_path)))
        for base in missing:
            lb.insert(tk.END, '\u2717  {}   (no image found \u2014 skipped)'.format(base))
            lb.itemconfig(tk.END, fg='#AA0000')

        btn_frame = tk.Frame(dlg)
        btn_frame.pack(fill=tk.X, padx=10, pady=10)

        def _proceed():
            confirmed[0] = True
            dlg.destroy()

        tk.Button(btn_frame, text='Proceed', width=12, bg='#2E7D32', fg='white',
                  font=('Arial', 10, 'bold'), command=_proceed).pack(side=tk.LEFT, padx=4)
        tk.Button(btn_frame, text='Cancel', width=10,
                  font=('Arial', 10), command=dlg.destroy).pack(side=tk.LEFT)

        dlg.update_idletasks()
        dlg.minsize(500, 200)
        self.root.wait_window(dlg)
        if not confirmed[0]:
            return

        scale_pct   = self._scale_pct_var.get()
        ch_names    = [self.ch_name_var[0].get(), self.ch_name_var[1].get(),
                       self.ch3_name_var.get()]
        sens_by_ch  = {0: 2.5, 1: 2.5}        # diagnostics sigma for the refiner CSV (classical engine removed)
        _pfx        = os.path.basename(directory.rstrip(os.sep)).replace(' ', '')
        out_meas    = os.path.join(directory, f'{_pfx}_soma_measurements.xlsx')
        out_scaled  = os.path.join(directory, f'{_pfx}_intensity_scaled_soma_measurements.xlsx')
        out_refiner = os.path.join(directory, f'{_pfx}_refiner_training_data.csv')

        # Shared column headers for the combined workbooks.  Each NPZ carries
        # its own channel names (GFP-AKT1-Cre in one animal, GFP-AKT3-Cre in
        # another); a channel whose name varies across the batch is headed by
        # all of its names ("GFP-AKT1-Cre / GFP-AKT3-Cre") while the per-cell
        # 'Channel detected' / 'Source channel' columns keep each file's own.
        header_ch_names = batch_header_channel_names(pairs, ch_names)

        # A batch always starts fresh: outputs of an earlier run in this
        # directory are set aside under a timestamped name rather than appended
        # to, so the workbook's headers and rows describe exactly this run.
        _stamp = time.strftime('%Y%m%d-%H%M')
        for _p in (out_meas, out_scaled, out_refiner):
            if os.path.isfile(_p):
                _root, _ext = os.path.splitext(_p)
                try:
                    os.replace(_p, '{} (previous {}){}'.format(_root, _stamp, _ext))
                except OSError as e:
                    messagebox.showerror(
                        'Batch Process',
                        'Cannot set aside the previous output\n{}\n\n'
                        'Close it in Excel and run the batch again.\n({})'
                        .format(os.path.basename(_p), e))
                    return

        self.btn_batch.config(state=tk.DISABLED)
        self.btn_export.config(state=tk.DISABLED)

        def worker():
            n_ok   = 0
            errors = []

            for i, (img_path, npz_path) in enumerate(pairs):
                short = os.path.basename(img_path)[:35]
                self.root.after(0, lambda s=short, i=i: [
                    self.lbl_status.config(
                        text='Batch {}/{}: {}'.format(i + 1, len(pairs), s)),
                    self.lbl_batch.config(
                        text='{}/{}'.format(i + 1, len(pairs))),
                ])
                try:
                    # Load image
                    ext = os.path.splitext(img_path)[1].lower()
                    if ext == '.lsm':
                        stack, px_um, ch_info_lsm = load_lsm(img_path)
                        # Auto-assign channels for LSM too
                        order = auto_assign_channels(ch_info_lsm) if ch_info_lsm else None
                        if order and max(order) < stack.shape[0]:
                            stack = stack[order]
                    else:
                        with czifile.CziFile(img_path) as czi:
                            stack   = czi.asarray()
                            px_um   = get_pixel_size(czi.metadata())
                            ch_info = read_czi_channel_info(czi.metadata())
                        # Apply the same channel reordering as interactive mode
                        order = auto_assign_channels(ch_info) if ch_info else None
                        if order and max(order) < stack.shape[0]:
                            stack = stack[order]

                    # Reconstruct soma_by_ch and rejected from NPZ
                    data       = np.load(npz_path, allow_pickle=False)
                    n_rois     = int(data['ch_idx'].shape[0])
                    soma_by_ch = {}
                    rejected   = set()
                    src_arr    = data['source'] if 'source' in data else None

                    # Per-file channel names + genotype + animal ID from the NPZ
                    # metadata, so a directory can hold many genotypes/animals.
                    # Older NPZs lack an animal ID → fall back to a filename guess.
                    file_ch_names = ch_names
                    file_geno     = None
                    file_animal   = guess_animal_id(img_path)
                    if 'meta_json' in data:
                        try:
                            _meta = json.loads(str(data['meta_json']))
                            _mn   = _meta.get('ch_names')
                            if _mn and len(_mn) >= 3:
                                file_ch_names = _mn
                            file_geno = _meta.get('genotype_runs')
                            if _meta.get('animal_id'):
                                file_animal = _meta['animal_id']
                        except Exception:
                            pass

                    for flat_idx in range(n_rois):
                        ch = int(data['ch_idx'][flat_idx])
                        if ch not in soma_by_ch:
                            soma_by_ch[ch] = []
                        si = len(soma_by_ch[ch])
                        soma_by_ch[ch].append({
                            'z':              int(data['z'][flat_idx]),
                            'z_min':          int(data['z_min'][flat_idx]),
                            'z_max':          int(data['z_max'][flat_idx]),
                            'cy':             float(data['cy'][flat_idx]),
                            'cx':             float(data['cx'][flat_idx]),
                            'area_um2':       float(data['area_um2'][flat_idx]),
                            'cell_mean':      float(data['cell_mean'][flat_idx]),
                            'peak_intensity': float(data['peak_intensity'][flat_idx]),
                            'mask2d':         data['mask_{}'.format(flat_idx)].astype(bool),
                            'source':         str(src_arr[flat_idx])
                                              if src_arr is not None else 'auto',
                        })
                        if data['rejected'][flat_idx]:
                            rejected.add((ch, si))

                    czi_fname = os.path.basename(img_path)

                    measure_and_export(stack, px_um, czi_fname, soma_by_ch, rejected,
                                       file_ch_names, out_meas,
                                       genotype_runs=file_geno, animal_id=file_animal,
                                       header_ch_names=header_ch_names)
                    measure_and_export(stack, px_um, czi_fname, soma_by_ch, rejected,
                                       file_ch_names, out_scaled, scale_pct=scale_pct,
                                       genotype_runs=file_geno, animal_id=file_animal,
                                       header_ch_names=header_ch_names)
                    export_refiner_data(stack, px_um, czi_fname, soma_by_ch, rejected,
                                        sens_by_ch, file_ch_names, out_refiner)
                    n_ok += 1

                except Exception as e:
                    errors.append('{}: {}'.format(os.path.basename(img_path), str(e)))

            def _done():
                self.btn_batch.config(state=tk.NORMAL)
                self.btn_export.config(state=tk.NORMAL)
                summary = 'Batch complete: {}/{} files processed.'.format(
                    n_ok, len(pairs))
                if errors:
                    summary += '\n\nErrors:\n' + '\n'.join(errors[:10])
                    if len(errors) > 10:
                        summary += '\n…and {} more.'.format(len(errors) - 10)
                self.lbl_status.config(
                    text='Batch done: {}/{}'.format(n_ok, len(pairs)))
                self.lbl_batch.config(text='Last run: {}/{} OK'.format(n_ok, len(pairs)))
                messagebox.showinfo('Batch complete', summary)

            self.root.after(0, _done)

        threading.Thread(target=worker, daemon=True).start()

    # ── ROI save / load ───────────────────────────────────────────────────────

    def _save_rois(self, out_dir):
        basename = os.path.splitext(os.path.basename(self.czi_path))[0]
        out_path = os.path.join(out_dir, basename + '_ROIs.npz')

        # Per-file metadata: user-defined channel names (virus/antibody) and the
        # formatted mouse genotype, so loading this ROI set re-defines them and
        # batch runs can carry different genotypes per file in one directory.
        meta = {
            'ch_names': [self.ch_name_var[0].get(),
                         self.ch_name_var[1].get(),
                         self.ch3_name_var.get()],
            'genotype_runs': self._read_genotype_runs(),
            'animal_id': self.animal_id_var.get().strip(),
            'engine': 'cellpose',
            'engine_params': {'model': cellpose_engine.model_name(),
                              'cellprob': self.cp_vars['cellprob'].get(),
                              'diam': self.cp_vars['diam'].get(),
                              'min_slices': self.cp_vars['min_slices'].get(),
                              'open_um': self.cp_vars['open_um'].get(),
                              'min_area_um2': self.cp_vars['min_area'].get(),
                              'max_area_um2': self.cp_vars['max_area'].get(),
                              'exclude_z_edges': bool(self.cp_vars['exclude_z_edges'].get()),
                              'surface_keep_um2': cellpose_engine.SELECT_DEFAULTS.get('surface_keep_um2'),
                              'rep_rule': self._rep_rule_key(),
                              'min_pick': self.cp_vars['min_pick'].get(),
                              'min_contrast': [v.get() for v in self.cp_contrast_var],
                              'breadth_um': [v.get() for v in self.breadth_var]},
        }
        arrays = {
            'czi_fname':   np.array(os.path.basename(self.czi_path)),
            'px_um':       np.array(self.px_um),
            'stack_shape': np.array(self.stack.shape),
            'meta_json':   np.array(json.dumps(meta)),
        }
        ch_list, z_list, z_min_list, z_max_list = [], [], [], []
        cy_list, cx_list, area_list, mean_list, peak_list, rej_list, src_list = \
            [], [], [], [], [], [], []

        flat_idx = 0
        for ch_idx, soma_list in self.soma_by_ch.items():
            for si, s in enumerate(soma_list):
                arrays['mask_{}'.format(flat_idx)] = s['mask2d']
                ch_list.append(ch_idx)
                z_list.append(s['z'])
                z_min_list.append(s.get('z_min', s['z']))
                z_max_list.append(s.get('z_max', s['z']))
                cy_list.append(s['cy'])
                cx_list.append(s['cx'])
                area_list.append(s['area_um2'])
                mean_list.append(s['cell_mean'])
                peak_list.append(s.get('peak_intensity', 0.0))
                rej_list.append((ch_idx, si) in self.rejected)
                src_list.append(s.get('source', 'auto'))
                flat_idx += 1

        arrays.update({
            'ch_idx':         np.array(ch_list,   dtype=np.int32),
            'z':              np.array(z_list,     dtype=np.int32),
            'z_min':          np.array(z_min_list, dtype=np.int32),
            'z_max':          np.array(z_max_list, dtype=np.int32),
            'cy':             np.array(cy_list),
            'cx':             np.array(cx_list),
            'area_um2':       np.array(area_list),
            'cell_mean':      np.array(mean_list),
            'peak_intensity': np.array(peak_list),
            'rejected':       np.array(rej_list,  dtype=bool),
            'source':         np.array(src_list),
        })
        np.savez_compressed(out_path, **arrays)
        return out_path

    def _load_rois(self):
        path = filedialog.askopenfilename(
            title='Select ROI set file',
            filetypes=[('ROI set', '*_ROIs.npz'), ('All npz', '*.npz')])
        if not path:
            return
        try:
            data = np.load(path, allow_pickle=False)
        except Exception as e:
            messagebox.showerror('Load error', str(e))
            return

        if self.stack is not None and 'stack_shape' in data:
            saved = tuple(data['stack_shape'].tolist())
            if saved != self.stack.shape:
                if not messagebox.askyesno(
                        'Shape mismatch',
                        'Saved shape {} ≠ current {}.\nLoad anyway?'.format(
                            saved, self.stack.shape)):
                    return

        n          = int(data['ch_idx'].shape[0])
        soma_by_ch = {}
        rejected   = set()
        src_arr    = data['source'] if 'source' in data else None

        for flat_idx in range(n):
            ch = int(data['ch_idx'][flat_idx])
            if ch not in soma_by_ch:
                soma_by_ch[ch] = []
            si = len(soma_by_ch[ch])
            soma_by_ch[ch].append({
                'z':              int(data['z'][flat_idx]),
                'z_min':          int(data['z_min'][flat_idx]),
                'z_max':          int(data['z_max'][flat_idx]),
                'cy':             float(data['cy'][flat_idx]),
                'cx':             float(data['cx'][flat_idx]),
                'area_um2':       float(data['area_um2'][flat_idx]),
                'cell_mean':      float(data['cell_mean'][flat_idx]),
                'peak_intensity': float(data['peak_intensity'][flat_idx]),
                'mask2d':         data['mask_{}'.format(flat_idx)].astype(bool),
                'source':         str(src_arr[flat_idx]) if src_arr is not None else 'auto',
            })
            if data['rejected'][flat_idx]:
                rejected.add((ch, si))

        self.soma_by_ch = soma_by_ch
        self.rejected   = rejected

        # Restore user-defined channel names and mouse genotype, if present.
        if 'meta_json' in data:
            try:
                meta  = json.loads(str(data['meta_json']))
                names = meta.get('ch_names') or []
                if len(names) >= 3:
                    self.ch_name_var[0].set(names[0])
                    self.ch_name_var[1].set(names[1])
                    self.ch3_name_var.set(names[2])
                self._set_genotype_runs(meta.get('genotype_runs') or [])
                if meta.get('animal_id'):
                    self.animal_id_var.set(meta['animal_id'])
            except Exception:
                pass

        self._update_counts()
        self._update_info_strip()
        self._redraw()

    def _on_close(self):
        try:
            if getattr(self, '_jump_server', None) is not None:
                self._jump_server.stop()
        except Exception:
            pass
        plt.close('all')
        self.root.destroy()

    # ── Click-to-jump IPC handlers ────────────────────────────────────────────

    def _bring_to_front(self):
        """Force the window to the foreground. Windows blocks a background app
        from stealing focus with plain lift()/focus_force(), so briefly set the
        window topmost (which does raise it) and then release topmost so it does
        not stay pinned above everything else."""
        try:
            r = self.root
            r.deiconify()                       # restore if minimized
            r.lift()
            r.attributes('-topmost', True)
            r.update_idletasks()
            r.focus_force()
            # Release the always-on-top flag shortly after; by then the window
            # is already the frontmost, so it stays visible for the blink.
            r.after(900, lambda: r.attributes('-topmost', False))
        except Exception:
            pass

    def _handle_jump_command(self, msg):
        """Dispatch an inbound jump command from the Group Analysis GUI.

        Runs on the Tk main thread (the JumpServer drains its queue via
        root.after).  Tries to keep error handling forgiving so the user
        sees a meaningful messagebox instead of a stack trace if something
        like the .czi file is missing.
        """
        try:
            action = msg.get('action')
            if action != 'open_cell':
                return  # unknown action — ignore silently

            target_file = msg.get('file')
            if not target_file:
                return

            # ── Resolve filename to a full path ──
            resolved = self._resolve_czi_path(
                target_file, hints=msg.get('search_dirs') or [])
            if resolved is None:
                messagebox.showwarning(
                    'Cell-jump',
                    'Could not find\n  {}\n\nin the cached search '
                    'directories.\nPlease pick the folder containing this '
                    '.czi file.'.format(target_file))
                pick = filedialog.askdirectory(
                    title='Locate folder containing {}'.format(target_file))
                if not pick:
                    return
                self._jump_search_cache.insert(0, pick)
                resolved = self._resolve_czi_path(target_file, hints=[pick])
                if resolved is None:
                    messagebox.showerror(
                        'Cell-jump',
                        '{} not found in {}'.format(target_file, pick))
                    return

            # Cache the resolved directory for future jumps in this session.
            resolved_dir = os.path.dirname(resolved)
            if resolved_dir and resolved_dir not in self._jump_search_cache:
                self._jump_search_cache.insert(0, resolved_dir)

            # Bring the soma GUI to the front so the user sees the jump.
            self._bring_to_front()

            # ── Load the file if not already loaded ──
            if self.czi_path != resolved:
                self._load_file(resolved, prompt=False)

            # ── Auto-load matching ROI file (best-effort) ──
            self._auto_load_rois_for(resolved)

            if self.stack is None:
                return
            nz = self.stack.shape[1]

            # ── Identify the cell ──
            # The exported 'Slice' column is 1-based (plane index + 1);
            # centroids are µm = pixel × px_um.
            target_slice = msg.get('slice')
            z0 = (int(target_slice) - 1
                  if isinstance(target_slice, (int, float)) else None)
            cx_um, cy_um = msg.get('centroid_x_um'), msg.get('centroid_y_um')
            soma_num = msg.get('soma_num')
            cond     = str(msg.get('condition') or '').strip()
            label    = ('Soma #{}'.format(soma_num) if soma_num is not None else 'jumped cell')
            if cond and cond.lower() not in ('nan', 'none', ''):
                label += '  ' + cond

            key = None
            if cx_um is not None and cy_um is not None and self.px_um > 0:
                key = self._find_soma_near(float(cx_um) / self.px_um,
                                           float(cy_um) / self.px_um, z0)
            if key is not None:
                ch, si = key
                s      = self.soma_by_ch[ch][si]
                cx, cy = float(s['cx']), float(s['cy'])
                z_show = int(s['z'])
                rej    = key in self.rejected
                self._jump_target = {'key': key, 'cx': cx, 'cy': cy,
                                     'label': label + (' (rejected)' if rej else '')}
                self.ch_var.set(ch)
                status = ('Jumped to {} — Ch{} {}, plane {} of {}{}. '
                          'Magenta outline + label mark it; Esc clears, double-click resets the view.'
                          .format(label, ch + 1, self.ch_name_var[ch].get(), z_show + 1, nz,
                                  ' — currently REJECTED' if rej else ''))
            else:
                cx = float(cx_um) / self.px_um if cx_um is not None else self.stack.shape[3] / 2
                cy = float(cy_um) / self.px_um if cy_um is not None else self.stack.shape[2] / 2
                z_show = z0 if (z0 is not None and 0 <= z0 < nz) else max(self.current_z, 0)
                self._jump_target = {'key': None, 'cx': cx, 'cy': cy,
                                     'label': label + '  (no matching outline)'}
                status = ('Jumped to {} — no outline within 20 µm of the exported centroid '
                          '(is the right ROI set loaded?). The crosshair marks the exported position.'
                          .format(label))

            # ── Navigate: the cell's own plane, zoomed to a JUMP_VIEW_UM window ──
            if 0 <= z_show < nz:
                self.slice_var.set(z_show + 1)
                self.current_z = z_show
                self.lbl_slice.config(text='z={}'.format(z_show + 1))
            half = (self.JUMP_VIEW_UM / 2.0) / self.px_um
            self._set_view(cx - half, cx + half, cy + half, cy - half)
            self._redraw()
            self.lbl_status.config(text=status)
            self._bring_to_front()
            if key is not None:
                self._flash_soma_outline(key[0], key[1])      # attention pulse on top
            else:
                self._flash_marker_at(cx, cy)

        except Exception as e:
            import traceback
            traceback.print_exc()
            try:
                messagebox.showerror('Cell-jump error', str(e))
            except Exception:
                pass

    def _resolve_czi_path(self, basename, hints=()):
        """Find the full path of a .czi/.lsm file by basename.

        Search order:
          1. The current cached directories (most recent first)
          2. Any directories passed in `hints` (from the group GUI)
          3. Recursive descent of each cached dir (one level deep)
        Returns the first hit or None.
        """
        base = os.path.basename(basename)

        candidates = []
        for d in list(self._jump_search_cache) + list(hints):
            if d and d not in candidates:
                candidates.append(d)

        # Direct hits first (cheap).
        for d in candidates:
            p = os.path.join(d, base)
            if os.path.isfile(p):
                return p

        # One-level recursive scan (covers the common "datafiles in a
        # subfolder of the excel" layout without traversing huge trees).
        for d in candidates:
            if not os.path.isdir(d):
                continue
            try:
                for sub in os.listdir(d):
                    full = os.path.join(d, sub)
                    if os.path.isdir(full):
                        p = os.path.join(full, base)
                        if os.path.isfile(p):
                            return p
            except OSError:
                continue
        return None

    def _auto_load_rois_for(self, czi_path):
        """Auto-load <basename>_ROIs.npz from the same directory, if present.

        Silently does nothing when no companion file exists; surfaces only
        catastrophic load failures via a console warning.
        """
        try:
            base, _ext = os.path.splitext(czi_path)
            roi_path = base + '_ROIs.npz'
            if not os.path.isfile(roi_path):
                return
            data = np.load(roi_path, allow_pickle=False)
            n = int(data['ch_idx'].shape[0])
            soma_by_ch = {}
            rejected   = set()
            src_arr    = data['source'] if 'source' in data else None
            for flat_idx in range(n):
                ch = int(data['ch_idx'][flat_idx])
                soma_by_ch.setdefault(ch, [])
                si = len(soma_by_ch[ch])
                soma_by_ch[ch].append({
                    'z':              int(data['z'][flat_idx]),
                    'z_min':          int(data['z_min'][flat_idx]),
                    'z_max':          int(data['z_max'][flat_idx]),
                    'cy':             float(data['cy'][flat_idx]),
                    'cx':             float(data['cx'][flat_idx]),
                    'area_um2':       float(data['area_um2'][flat_idx]),
                    'cell_mean':      float(data['cell_mean'][flat_idx]),
                    'peak_intensity': float(data['peak_intensity'][flat_idx]),
                    'mask2d':         data['mask_{}'.format(flat_idx)].astype(bool),
                    'source':         str(src_arr[flat_idx]) if src_arr is not None else 'auto',
                })
                if data['rejected'][flat_idx]:
                    rejected.add((ch, si))
            self.soma_by_ch = soma_by_ch
            self.rejected   = rejected
            self._update_counts()
        except Exception as e:
            sys.stderr.write('Auto-ROI load failed: {}\n'.format(e))

    JUMP_VIEW_UM = 100.0     # click-to-jump: width of the zoomed window around the cell

    def _find_soma_near(self, px_x, px_y, z0=None, tol_um=20.0):
        """(ch, si) of the loaded soma nearest to a pixel position (any channel),
        within tol_um and within ±2 planes of z0 when given; non-rejected somas
        are preferred over rejected ones at similar distance."""
        best = None
        for ch, lst in self.soma_by_ch.items():
            for si, s in enumerate(lst):
                if z0 is not None:
                    zmin, zmax = s.get('z_min', s['z']), s.get('z_max', s['z'])
                    if not (zmin - 2 <= z0 <= zmax + 2):
                        continue
                d = float(np.hypot(s['cx'] - px_x, s['cy'] - px_y)) * self.px_um
                if d > tol_um:
                    continue
                score = d + (5.0 if (ch, si) in self.rejected else 0.0)
                if best is None or score < best[0]:
                    best = (score, (ch, si))
        return best[1] if best else None

    def _draw_jump_highlight(self, z):
        """Persistent click-to-jump marker: thick magenta outline (dimmed when
        the current plane is outside the soma's z-range), or a crosshair when
        no outline matched, plus a label with the exported soma number."""
        jt = self._jump_target
        cx, cy = jt['cx'], jt['cy']
        if jt.get('key') is not None:
            ch, si = jt['key']
            lst = self.soma_by_ch.get(ch, [])
            if si >= len(lst):
                return
            s = lst[si]
            on_plane = (z == -1) or (s.get('z_min', s['z']) <= z <= s.get('z_max', s['z']))
            for cnt in self._soma_contours(s):
                self.ax.plot(cnt[:, 1], cnt[:, 0], color='#FF00FF',
                             lw=3.2 if on_plane else 1.4,
                             alpha=0.95 if on_plane else 0.45, zorder=90)
        else:
            self.ax.plot([cx], [cy], marker='+', color='#FF00FF',
                         ms=26, mew=2.5, zorder=90)
        self.ax.annotate(jt['label'], xy=(cx, cy), xytext=(30, -30),
                         textcoords='offset points', color='white',
                         fontsize=9, fontweight='bold', zorder=95,
                         bbox=dict(boxstyle='round,pad=0.3', fc='#FF00FF',
                                   ec='none', alpha=0.92),
                         arrowprops=dict(arrowstyle='-|>', color='#FF00FF', lw=2))

    def _clear_jump_highlight(self):
        if self._jump_target is not None:
            self._jump_target = None
            if self.stack is not None:
                self._redraw()

    def _flash_nearest_soma(self, cx_um, cy_um, z=None, channel_hint=None,
                            radius_um=8.0):
        """Locate the loaded soma nearest to the (μm, μm) target and flash it.

        Uses the soma centres recorded in pixel coordinates and converts
        the target via self.px_um.  Channel preference is taken when
        provided; otherwise both channels are searched.  Flashes by
        drawing a thick magenta outline that pulses for ~1.6 s.
        """
        if not self.soma_by_ch or self.px_um <= 0:
            return

        target_px_x = cx_um / self.px_um
        target_px_y = cy_um / self.px_um
        tol_px      = radius_um / self.px_um

        ch_order = ([int(channel_hint)] if channel_hint is not None
                    else list(self.soma_by_ch.keys()))

        best = None  # (distance, ch, si)
        for ch in ch_order:
            for si, s in enumerate(self.soma_by_ch.get(ch, [])):
                if (ch, si) in self.rejected:
                    continue
                if z is not None:
                    z_min = s.get('z_min', s['z'])
                    z_max = s.get('z_max', s['z'])
                    if not (z_min <= z <= z_max):
                        continue
                d = ((s['cx'] - target_px_x) ** 2 +
                     (s['cy'] - target_px_y) ** 2) ** 0.5
                if best is None or d < best[0]:
                    best = (d, ch, si)

        if best is None or best[0] > tol_px * 2.5:
            # No match — fall back to marking the bare coordinates.
            self._flash_marker_at(target_px_x, target_px_y)
            return

        _, ch, si = best
        self._flash_soma_outline(ch, si)

    def _flash_soma_outline(self, ch, si, n_pulses=4, period_ms=200):
        """Draw a thick magenta outline over an existing ROI, pulsing
        between two alphas so the user notices it on first jump."""
        soma = self.soma_by_ch.get(ch, [None])[si] \
            if si < len(self.soma_by_ch.get(ch, [])) else None
        if soma is None:
            return
        cnts = measure.find_contours(soma['mask2d'].astype(np.float32), 0.5)
        artists = []
        for cnt in cnts:
            patch = MplPolygon(cnt[:, ::-1], closed=True,
                               fill=False, edgecolor='#FF00FF',
                               linewidth=4.0, alpha=0.95, zorder=100)
            self.ax.add_patch(patch)
            artists.append(patch)

        # (The jump handler has already navigated and zoomed to the cell.)
        self.canvas.draw_idle()

        # Pulse alpha to draw the eye.
        def pulse(step):
            if not artists:
                return
            for a in artists:
                try:
                    a.set_alpha(0.95 if step % 2 == 0 else 0.20)
                except Exception:
                    pass
            try:
                self.canvas.draw_idle()
            except Exception:
                return
            if step < n_pulses * 2 - 1:
                self.root.after(period_ms, lambda: pulse(step + 1))
            else:
                # Cleanup: remove the overlay so subsequent redraws are clean.
                self.root.after(period_ms * 3, lambda: self._remove_flash(artists))
        pulse(0)

    def _flash_marker_at(self, px_x, px_y, n_pulses=4, period_ms=200):
        """Flash a crosshair marker when no matching ROI was found."""
        try:
            line1 = self.ax.axhline(px_y, color='#FF00FF', lw=2, zorder=100)
            line2 = self.ax.axvline(px_x, color='#FF00FF', lw=2, zorder=100)
            artists = [line1, line2]
        except Exception:
            return
        self.canvas.draw_idle()

        def pulse(step):
            for a in artists:
                try:
                    a.set_alpha(0.9 if step % 2 == 0 else 0.2)
                except Exception:
                    pass
            try:
                self.canvas.draw_idle()
            except Exception:
                return
            if step < n_pulses * 2 - 1:
                self.root.after(period_ms, lambda: pulse(step + 1))
            else:
                self.root.after(period_ms * 3, lambda: self._remove_flash(artists))
        pulse(0)

    def _remove_flash(self, artists):
        """Remove temporary highlight artists and redraw."""
        for a in artists:
            try:
                a.remove()
            except Exception:
                pass
        try:
            self.canvas.draw_idle()
        except Exception:
            pass


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    if sys.stdout is None or sys.stderr is None:          # pythonw: no console
        _log = open(os.path.join(tempfile.gettempdir(), 'soma_detector_3.0.log'), 'a', encoding='utf-8')
        sys.stdout = sys.stdout or _log
        sys.stderr = sys.stderr or _log
    czi_path = sys.argv[1] if len(sys.argv) > 1 else None
    root = tk.Tk()
    root.geometry('1400x1000')
    _app = SomaDetectorGUI(root, czi_path=czi_path)     # kept alive for the mainloop
    root.mainloop()


if __name__ == '__main__':
    main()
