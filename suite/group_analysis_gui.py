"""
Group Analysis GUI — Statistical tests on classified soma measurements

Test of record: the clustered-data model of Moen, Fricano-Kugler, Luikart &
O'Malley (2016, PLoS One 11:e0146721) as used in Prina, Goyette et al. (2026,
Brain) — a linear mixed-effects model with a random intercept for each animal
(cells nested within animals), fitted by maximum likelihood, standard errors
cluster-robust on animal, Scheffé-adjusted pairwise contrasts (mixed_stats.py;
validated against the published Stata output).  Its Scheffé p-values drive the
significance brackets on the plot.

Descriptive secondary tests (scipy.stats, Bonferroni-corrected across the
pairwise family):
  * p_animal — Welch's t-test on per-animal means
  * p_cell   — Welch's t-test on cell-level values (overstates significance)
  * p_MW     — Mann-Whitney U on cell-level values (non-parametric)
If the mixed model cannot be fitted, the animal-level Bonferroni p-value
drives the brackets (cell-level when a group has fewer than two animals,
with that caveat noted in the results popup).

Conditions are detected dynamically from the filenames in the loaded file.
Standard filename pattern:  CC{animal}_{date}_{#}_{GENOTYPE}_{reporter}_...

WT cells:   redRubi reporter OR Cre- classification
KO cells:   any *Cre reporter → label derived from genotype token
  e.g.  Pten / pRubiGCre       → Pten KO
        PtenAKT1 / pRubiGAKT1Cre → Pten;Akt1 KO
        PtenAKT3 / pRubiGAKT3Cre → Pten;Akt3 KO
"""

import os
import sys
import subprocess
import re
import json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import numpy as np
import pandas as pd
try:
    import mixed_stats                 # the Brain-paper model (same folder): random-intercept
except ImportError:                    # mixed model, cluster-robust SE, Scheffé contrasts
    mixed_stats = None
import openpyxl
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure   # not pyplot: a pyplot figure owns a hidden
                                       # Tk root that keeps the process alive after
                                       # the window is closed
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy import stats as _scstats
from itertools import combinations
import cell_jump_ipc


def _launch_soma_detector():
    """Launch the Soma Detector application so the click-to-jump feature works
    even when it isn't already open. Uses the packaged executable when frozen
    (installed next to this one) or the .py script during development. Returns
    True if a launch was attempted."""
    try:
        if getattr(sys, 'frozen', False):
            exe_dir = os.path.dirname(sys.executable)
            for name in ('SomaDetector.exe',
                         'soma_area_and_intensity_gui_2.1.exe'):
                cand = os.path.join(exe_dir, name)
                if os.path.isfile(cand):
                    subprocess.Popen([cand], cwd=exe_dir)
                    return True
            return False
        here   = os.path.dirname(os.path.abspath(__file__))
        script = None
        for name in ('soma_area_and_intensity_gui_3.0.py',      # Cellpose engine (v3.0)
                     'soma_area_and_intensity_gui_2.1.py'):
            cand = os.path.join(here, name)
            if os.path.isfile(cand):
                script = cand
                break
        if script is None:
            return False
        # Prefer the Cellpose environment so v3.0 starts with its AI engine:
        # the portable python312 shipped with the program, then .venv-cellpose,
        # otherwise the interpreter running this GUI.
        portable_pyw = os.path.join(here, 'python312', 'pythonw.exe')
        venv_pyw = os.path.join(here, '.venv-cellpose', 'Scripts', 'pythonw.exe')
        pyw = os.path.join(os.path.dirname(sys.executable), 'pythonw.exe')
        launcher = next((p for p in (portable_pyw, venv_pyw, pyw) if os.path.isfile(p)),
                        sys.executable)
        subprocess.Popen([launcher, script], cwd=here)
        return True
    except Exception:
        return False

# ── Default color palette ─────────────────────────────────────────────────────

_WT_COLOR = '#AAAAAA'
_KO_PALETTE = [
    '#2ca02c', '#d62728', '#9467bd', '#1f77b4',
    '#ff7f0e', '#17becf', '#8c564b', '#e377c2',
]

# Phenotype / Classification values that identify Cre-negative (WT) cells.
# Covers both the old 'Cre-' standard value and common user-defined counter-only names.
_WT_PHENOTYPE_VALS = frozenset([
    'cre-', 'wt', 'mcherry', 'mcherry only',
    'counter-only', 'counter only', 'gfp only', 'gfp-only',
])


# ── Filename parsing (shared with soma_analysis_gui) ──────────────────────────

_SKIP_RE = re.compile(r'^(CC\d+|\d{6}|\d{1,2}[FM]?|\d{2})$', re.IGNORECASE)
_INTERFIELD_RE = re.compile(r'^\d+[FM]?$|^[FM]$', re.IGNORECASE)


def _parse_filename_parts(filename):
    """
    Extract (genotype_token, reporter_token) from a CZI/LSM filename.
    Standard pattern:  CC{animal}_{date}_{#}_{GENOTYPE}_{reporter}_...
    """
    base  = os.path.splitext(os.path.basename(str(filename)))[0]
    parts = [p for p in base.replace('#', '_').split('_') if p]

    genotype = None
    for p in parts:
        if _SKIP_RE.match(p):
            continue
        if genotype is None:
            genotype = p
        else:
            if _INTERFIELD_RE.match(p):
                continue
            return genotype, p

    clean = parts
    if len(clean) > 4:
        return clean[3], clean[4]
    if len(clean) > 3:
        return clean[3], 'Unknown'
    return (genotype or 'Unknown'), 'Unknown'


def _genotype_label(gt_tok, fn=''):
    """Map a genotype token to a display label.

    Builds the label combinatorially from every floxed allele encoded in the
    genotype token, so compound genotypes like PtenRictorAkt1Akt2Akt3 do not
    collapse onto simpler ones. Only the genotype token is inspected for allele
    tokens — the full filename is used only as a fallback when the genotype
    token is missing or unparseable, because reporter strings like
    "pRubiGAKT1T2ACre" contain spurious AKT1 substrings.
    """
    gt = gt_tok.upper()

    if 'C57' in gt and 'PTEN' not in gt:
        return 'C57'

    parts = []
    if 'PTEN' in gt:
        parts.append('Pten')
    if 'RAPTOR' in gt:
        parts.append('Raptor')
    if 'RICTOR' in gt:
        parts.append('Rictor')
    if 'AKT1' in gt:
        parts.append('Akt1')
    if 'AKT2' in gt:
        parts.append('Akt2')
    if 'AKT3' in gt:
        parts.append('Akt3')

    if parts:
        return ';'.join(parts)

    # Fallback to filename when the genotype token itself yields nothing.
    fu = fn.upper()
    if 'RAPTOR' in fu:
        return 'Pten;Raptor'
    if 'RICTOR' in fu:
        return 'Pten;Rictor'
    if 'PTEN' in fu:
        return 'Pten'
    if 'C57' in fu:
        return 'C57'
    return gt_tok   # unknown — use raw token


def _condition_from_filename(fn, cls_value=None):
    """
    Derive experimental condition label from a filename.

    WT  : redRubi reporter  OR  Cre- / mCherry / WT classification value
    KO  : label from genotype token + reporter
    """
    gt_tok, rep_tok = _parse_filename_parts(fn)
    sv_rep = rep_tok.lower()

    # ── WT determination ──
    is_wt = ('redrubi' in sv_rep or sv_rep in ('rr', 'rb'))
    if not is_wt and cls_value is not None:
        sv_cls = str(cls_value).strip().lower()
        is_wt  = sv_cls in ('cre-', 'wt', 'mcherry', 'mcherry only')

    if is_wt:
        return 'WT'

    # ── KO — build label from genotype token and reporter ──
    gt_label = _genotype_label(gt_tok, fn)

    # Detect rescue constructs co-expressed with Cre.
    # Common reporter naming: pRubiGAKT1Cre, pRubiGAKT1T2ACre, pRubiGAKT3T2ACre …
    if 'akt3cre' in sv_rep or 'akt3t2acre' in sv_rep:
        suffix = ' + AKT3'
    elif 'akt1cre' in sv_rep or 'akt1t2acre' in sv_rep:
        suffix = ' + AKT1'
    elif 'cre' in sv_rep:
        suffix = ''
    else:
        suffix = f' ({rep_tok})'

    return f'{gt_label} KO{suffix}'


# ── Animal ID auto-detection ──────────────────────────────────────────────────

def _detect_animal(fn):
    """Derive an animal ID from a CZI/LSM filename.

    Preferred pattern: CC{cage}_{date}_{animal_#}_{genotype}_...
    Returns "CC{cage}_{animal_#}" so the cage prefix (which encodes the
    floxed-allele combination) is preserved and two animals from different
    cages that share a date and number are not collapsed into one ID.

    Older / unusual filenames fall back to the previous heuristics.
    """
    s    = str(fn)
    base = os.path.splitext(os.path.basename(s))[0]

    # Preferred: CC{cage}_{date}_{animal_#}[F|M]_... (e.g. CC4425_101725_3_…)
    m = re.match(r'^(CC\d+)_\d{6,8}_(\d+[A-Za-z]?)(?:_|$)', base)
    if m:
        return f'{m.group(1)}_{m.group(2)}'

    # Fallbacks (older filename styles)
    stripped = re.sub(r'^CC\d+_', '', base)
    m = re.match(r'(\d{6,8}#\d+)', stripped)
    if m:
        return m.group(1)
    m = re.match(r'(\d{6,8}_\w+?)_', stripped)
    if m:
        return m.group(1)
    m = re.match(r'^([^_]+)', base)
    return m.group(1) if m else base


# ── Statistics helpers ────────────────────────────────────────────────────────

def _sig_label(p):
    if   p < 0.0001: return '****'
    elif p < 0.001:  return '***'
    elif p < 0.01:   return '**'
    elif p < 0.05:   return '*'
    else:            return 'ns'


def _load_excel(path):
    """Try two-row header (formatted export) then fall back to one-row."""
    for hdr in (1, 0):
        try:
            df = pd.read_excel(path, header=hdr)
            df.columns = [str(c).strip() for c in df.columns]
            if any(c.lower() == 'file' for c in df.columns):
                return df
        except Exception:
            continue
    raise ValueError('Could not read Excel file.')


def _is_control(cond):
    """True for the wild-type / Cre-negative control condition (old data used
    the literal 'WT'; genotype-labelled data uses 'Cre-')."""
    return str(cond).strip() in ('WT', 'Cre-')


def _genotype_runs_from_cell(v):
    """Extract [(text, italic, superscript), ...] runs from a Genotype cell
    (openpyxl rich text or a plain string)."""
    if v is None:
        return []
    if isinstance(v, str):
        return [(v, False, False)] if v.strip() else []
    runs = []
    try:
        for b in v:                       # CellRichText → TextBlocks / strings
            if isinstance(b, str):
                if b:
                    runs.append((b, False, False))
                continue
            font   = getattr(b, 'font', None)
            italic = bool(getattr(font, 'i', False)) if font else False
            sup    = (getattr(font, 'vertAlign', None) == 'superscript') if font else False
            t = str(getattr(b, 'text', ''))
            if t:
                runs.append((t, italic, sup))
    except TypeError:
        return [(str(v), False, False)]
    return runs


def _mathtext_escape(s):
    """Escape text for use inside a matplotlib mathtext ($...$) expression."""
    s = s.replace('\\', r'\backslash ')
    for ch in ('{', '}', '_', '^', '#', '%', '&', '$', '~'):
        s = s.replace(ch, '\\' + ch)
    return s.replace(' ', r'\ ')


def _genotype_forms(v):
    """From a Genotype cell return (plain, mathtext).

    plain    : caret notation, e.g. 'Pten^flx/flx;Rictor^flx/flx'
    mathtext : matplotlib string, e.g. r'$\\mathit{Pten}^{flx/flx}$' — renders
               gene names italic and allele designations as superscript.
    """
    runs = _genotype_runs_from_cell(v)
    if not runs:
        return '', ''
    plain_parts, math_parts = [], []
    for text, italic, sup in runs:
        plain_parts.append(('^' + text) if sup else text)
        seg = _mathtext_escape(text)
        if italic:
            seg = r'\mathit{' + seg + '}'
        if sup:
            seg = '^{' + seg + '}'
        math_parts.append(seg)
    return ''.join(plain_parts).strip(), '$' + ''.join(math_parts) + '$'


def _read_genotype_labels(path):
    """Return {file: (plain, mathtext)} from a workbook's Genotype column,
    reading rich text so the formatting survives (pandas flattens it)."""
    out = {}
    try:
        wb = openpyxl.load_workbook(path, rich_text=True)
        ws = wb['Cells'] if 'Cells' in wb.sheetnames else wb.active
        hdr_row, hdr = None, {}
        for hr in (2, 1):
            row = {str(ws.cell(hr, c).value).strip(): c
                   for c in range(1, ws.max_column + 1)
                   if ws.cell(hr, c).value is not None}
            if 'File' in row and 'Genotype' in row:
                hdr_row, hdr = hr, row
                break
        if hdr_row is None:
            return out
        fcol, gcol = hdr['File'], hdr['Genotype']
        for r in range(hdr_row + 1, ws.max_row + 1):
            fv = ws.cell(r, fcol).value
            if fv is None:
                continue
            fv = str(fv)
            if fv not in out:
                out[fv] = _genotype_forms(ws.cell(r, gcol).value)
    except Exception:
        pass
    return out


def _draw_brackets(ax, pairs_results, x_pos, y_data_max, y_range, base_fs=11):
    if not pairs_results:
        return

    tick_h = y_range * 0.018
    step   = y_range * 0.075
    gap    = y_range * 0.04

    pairs_results = sorted(pairs_results,
                           key=lambda r: abs(x_pos[r[1]] - x_pos[r[0]]))

    heights = []
    for i, (c1, c2, p_adj) in enumerate(pairs_results):
        x1 = x_pos[c1]
        x2 = x_pos[c2]
        if x1 > x2:
            x1, x2 = x2, x1
        candidate = y_data_max + gap + step * i
        for prev_c1, prev_c2, _, prev_y in heights:
            px1, px2 = sorted([x_pos[prev_c1], x_pos[prev_c2]])
            if not (x2 < px1 or x1 > px2):
                candidate = max(candidate, prev_y + step)
        heights.append((c1, c2, p_adj, candidate))

    for c1, c2, p_adj, y in heights:
        x1  = x_pos[c1]
        x2  = x_pos[c2]
        lbl = _sig_label(p_adj)
        ax.plot([x1, x2], [y, y], 'k-', lw=1.0, clip_on=False)
        ax.plot([x1, x1], [y - tick_h, y], 'k-', lw=1.0, clip_on=False)
        ax.plot([x2, x2], [y - tick_h, y], 'k-', lw=1.0, clip_on=False)
        fs = base_fs if lbl != 'ns' else max(7, base_fs - 3)
        fw = 'bold' if lbl != 'ns' else 'normal'
        ax.text((x1 + x2) / 2, y + tick_h * 0.3, lbl,
                ha='center', va='bottom', fontsize=fs, fontweight=fw)


# ── Main application ──────────────────────────────────────────────────────────

class GroupAnalysisApp:

    def __init__(self, root):
        self.root = root
        self.root.title('Group Analysis — Statistical tests')
        self.root.geometry('1460x860')
        self.root.minsize(1050, 700)

        self.df             = None
        self.out_dir        = None
        self.file_col       = None
        self.cls_col        = None
        self._outcome_cols  = []
        self._pair_results  = {}
        self._animal_map    = {}
        self._raw_df        = None
        self._loaded_paths  = []   # all currently loaded file paths
        self._zrel_enable   = tk.BooleanVar(value=False)
        self._zrel_cutoff   = tk.DoubleVar(value=0.50)
        self._cre_window    = tk.IntVar(value=4)
        self._excluded_animals   = set()
        self._highlighted_animal = None
        self._show_all_bright    = False
        self._font_size          = tk.IntVar(value=11)
        self._model_data         = None
        self._results_win        = None
        # Click-to-identify: list of (PathCollection, [df_indices]) pairs
        # rebuilt every redraw so pick events can map a point back to its row.
        self._picker_data        = []
        self._cell_info_win      = None

        # Dynamic per-dataset condition info (rebuilt in _rebuild_condition_ui)
        self._conditions_order  = []   # ordered list of condition strings
        self._user_cond_order   = []   # display order chosen with the ▲/▼ buttons (kept across loads)
        self._pair_state        = {}   # frozenset({c1, c2}) -> bracket checkbox state (kept across rebuilds)
        self._condition_colors  = {}   # condition → hex color
        self._condition_include = {}   # condition → BooleanVar (include in plot/stats)
        self._all_pairs         = []   # list of (c1,c2) tuples
        self._pair_vars         = {}   # (c1,c2) → BooleanVar
        self.ref_var            = tk.StringVar(value='WT')

        self._grouping_cfg      = None   # manual grouping map (set by Grouping setup)

        self._build_ui()

    # ── Build UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        top = tk.Frame(self.root, pady=4)
        top.pack(fill=tk.X, padx=8)

        tk.Button(top, text='Open file(s)',
                  command=self._open_file,
                  bg='#1565C0', fg='white',
                  font=('Arial', 10, 'bold'), width=14).pack(side=tk.LEFT)

        tk.Button(top, text='Add file(s)',
                  command=self._add_file,
                  bg='#1976D2', fg='white',
                  font=('Arial', 10, 'bold'), width=12).pack(side=tk.LEFT, padx=2)

        tk.Button(top, text='Clear',
                  command=self._clear_files,
                  bg='#757575', fg='white',
                  font=('Arial', 9, 'bold'), width=7).pack(side=tk.LEFT, padx=2)

        tk.Button(top, text='Edit Animal IDs',
                  command=self._show_animal_editor,
                  bg='#E65100', fg='white',
                  font=('Arial', 9, 'bold'), width=14).pack(side=tk.LEFT, padx=4)

        tk.Button(top, text='Grouping setup',
                  command=self._open_grouping_setup,
                  bg='#6A1B9A', fg='white',
                  font=('Arial', 9, 'bold'), width=14).pack(side=tk.LEFT, padx=2)

        self.lbl_file = tk.Label(top, text='No file loaded',
                                  fg='#555', font=('Arial', 9), wraplength=600,
                                  justify='left')
        self.lbl_file.pack(side=tk.LEFT, padx=8)

        self.lbl_n = tk.Label(top, text='', fg='#333',
                               font=('Arial', 9, 'bold'))
        self.lbl_n.pack(side=tk.RIGHT, padx=12)

        main = tk.Frame(self.root)
        main.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        fig_frame = tk.Frame(main)
        fig_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.fig = Figure(figsize=(9, 6.5))
        self.ax  = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.13, right=0.96, bottom=0.28, top=0.90)

        self.canvas = FigureCanvasTkAgg(self.fig, fig_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, fig_frame).update()
        # Click-to-identify: connect pick events on scatter points.
        self.canvas.mpl_connect('pick_event', self._on_pick)

        # ── control panel (scrollable) ────────────────────────────────────────
        ctrl_outer = tk.Frame(main, width=385, relief=tk.RIDGE, bd=1)
        ctrl_outer.pack(side=tk.RIGHT, fill=tk.Y, padx=(8, 0))
        ctrl_outer.pack_propagate(False)

        ctrl_cv = tk.Canvas(ctrl_outer, highlightthickness=0,
                            bg=ctrl_outer.cget('bg'))
        ctrl_sb = ttk.Scrollbar(ctrl_outer, orient='vertical',
                                command=ctrl_cv.yview)
        ctrl_cv.configure(yscrollcommand=ctrl_sb.set)
        ctrl_sb.pack(side=tk.RIGHT, fill=tk.Y)
        ctrl_cv.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        ctrl = tk.Frame(ctrl_cv, padx=12, pady=10)
        _ctrl_win = ctrl_cv.create_window((0, 0), window=ctrl, anchor='nw')
        ctrl_cv.bind('<Configure>',
                     lambda e: ctrl_cv.itemconfig(_ctrl_win, width=e.width))
        ctrl.bind('<Configure>',
                  lambda e: ctrl_cv.configure(
                      scrollregion=ctrl_cv.bbox('all')))
        ctrl_cv.bind('<Enter>', lambda _: ctrl_cv.bind_all(
            '<MouseWheel>',
            lambda e: ctrl_cv.yview_scroll(int(-e.delta / 120), 'units')))
        ctrl_cv.bind('<Leave>', lambda _: ctrl_cv.unbind_all('<MouseWheel>'))

        tk.Label(ctrl, text='Analysis Options',
                 font=('Arial', 12, 'bold')).pack(anchor='w', pady=(0, 6))

        # outcome variable
        tk.Label(ctrl, text='Outcome variable:',
                 font=('Arial', 9, 'bold')).pack(anchor='w')
        self.outcome_var = tk.StringVar()
        self.outcome_cb  = ttk.Combobox(ctrl, textvariable=self.outcome_var,
                                         state='readonly', width=34)
        self.outcome_cb.pack(anchor='w', pady=(2, 8))
        self.outcome_cb.bind('<<ComboboxSelected>>', lambda e: self._update_plot())

        # reference condition — rebuilt dynamically by _rebuild_condition_ui
        tk.Label(ctrl, text='Reference condition (model base):',
                 font=('Arial', 9, 'bold')).pack(anchor='w')
        self._ref_frame = tk.Frame(ctrl)
        self._ref_frame.pack(anchor='w', padx=10, pady=(2, 6))
        # placeholder text until a file is loaded
        self._ref_placeholder = tk.Label(
            self._ref_frame, text='(load a file)',
            font=('Arial', 8), fg='#999')
        self._ref_placeholder.pack(anchor='w')

        # bracket comparisons — rebuilt dynamically
        self._brk_lf = tk.LabelFrame(ctrl, text='Show brackets for:',
                                      font=('Arial', 9), padx=6, pady=4)
        self._brk_lf.pack(fill=tk.X, pady=4)
        self._brk_placeholder = tk.Label(
            self._brk_lf, text='(load a file)',
            font=('Arial', 8), fg='#999')
        self._brk_placeholder.pack(anchor='w')

        tk.Frame(ctrl, height=6).pack()

        # ── depth filter ──────────────────────────────────────────────────────
        depth_frame = tk.LabelFrame(ctrl, text='Slice depth filter',
                                    font=('Arial', 9), padx=6, pady=4)
        depth_frame.pack(fill=tk.X, pady=4)
        tk.Label(depth_frame,
                 text='Antibody penetration decreases with depth.\n'
                      'Keep only cells in the surface-most fraction\n'
                      'of each Z-stack (0 = surface, 1 = deepest).',
                 font=('Arial', 8), fg='#555', justify='left').pack(anchor='w')
        zrel_row = tk.Frame(depth_frame)
        zrel_row.pack(anchor='w', pady=(4, 0))
        tk.Checkbutton(zrel_row, text='Limit to z_rel ≤',
                       variable=self._zrel_enable,
                       font=('Arial', 9),
                       command=self._on_zrel_change).pack(side=tk.LEFT)
        tk.Spinbox(zrel_row, from_=0.10, to=1.00, increment=0.05,
                   textvariable=self._zrel_cutoff,
                   width=6, font=('Arial', 9), format='%.2f',
                   command=self._on_zrel_change).pack(side=tk.LEFT, padx=4)
        tk.Label(zrel_row, text='(surface fraction)',
                 font=('Arial', 8), fg='#777').pack(side=tk.LEFT)
        self.lbl_zrel_n = tk.Label(depth_frame, text='',
                                    font=('Arial', 8), fg='#1565C0',
                                    wraplength=255, justify='left')
        self.lbl_zrel_n.pack(anchor='w')

        tk.Frame(ctrl, height=6).pack()

        # ── Cre- normalisation ────────────────────────────────────────────────
        cre_frame = tk.LabelFrame(ctrl, text='Cre⁻ normalisation window',
                                  font=('Arial', 9), padx=6, pady=4)
        cre_frame.pack(fill=tk.X, pady=4)
        tk.Label(cre_frame,
                 text='F/F0 = cell signal / mean signal of Cre⁻ cells\n'
                      'within ±N slices in the same file.\n'
                      'Applies to the measured channel — the third channel\n'
                      '(e.g. CaV2.1, pS6, pMARCKS).',
                 font=('Arial', 8), fg='#555', justify='left').pack(anchor='w')
        cre_row = tk.Frame(cre_frame)
        cre_row.pack(anchor='w', pady=(4, 0))
        tk.Label(cre_row, text='Window ±', font=('Arial', 9)).pack(side=tk.LEFT)
        tk.Spinbox(cre_row, from_=1, to=50, increment=1,
                   textvariable=self._cre_window,
                   width=4, font=('Arial', 9),
                   command=self._on_cre_window_change).pack(side=tk.LEFT, padx=2)
        tk.Label(cre_row, text='slices', font=('Arial', 8),
                 fg='#777').pack(side=tk.LEFT, padx=2)
        self.lbl_cre_cov = tk.Label(cre_frame, text='',
                                     font=('Arial', 8), fg='#1565C0',
                                     wraplength=255, justify='left')
        self.lbl_cre_cov.pack(anchor='w')

        tk.Frame(ctrl, height=6).pack()

        # ── Animal explorer ───────────────────────────────────────────────────
        anim_frame = tk.LabelFrame(ctrl, text='Animal explorer',
                                   font=('Arial', 9), padx=6, pady=4)
        anim_frame.pack(fill=tk.X, pady=4)
        tk.Label(anim_frame,
                 text='Click to highlight.  ✕ = excluded from analysis.',
                 font=('Arial', 8), fg='#555').pack(anchor='w')
        lb_frame = tk.Frame(anim_frame)
        lb_frame.pack(fill=tk.X, pady=(2, 0))
        sb_a = ttk.Scrollbar(lb_frame, orient='vertical')
        self._animal_lb = tk.Listbox(lb_frame, height=7,
                                      font=('Courier', 8),
                                      yscrollcommand=sb_a.set,
                                      selectmode='single',
                                      activestyle='none',
                                      exportselection=False)
        sb_a.config(command=self._animal_lb.yview)
        self._animal_lb.pack(side=tk.LEFT, fill=tk.X, expand=True)
        sb_a.pack(side=tk.RIGHT, fill=tk.Y)
        self._animal_lb.bind('<<ListboxSelect>>', self._on_animal_select)
        btn_row_a = tk.Frame(anim_frame)
        btn_row_a.pack(anchor='w', pady=(3, 0))
        tk.Button(btn_row_a, text='Exclude / Include',
                  font=('Arial', 8), width=16,
                  command=self._toggle_animal_exclude).pack(side=tk.LEFT, padx=(0, 4))
        tk.Button(btn_row_a, text='Clear all',
                  font=('Arial', 8), width=8,
                  command=self._clear_exclusions).pack(side=tk.LEFT, padx=(0, 4))
        tk.Button(btn_row_a, text='Show all',
                  font=('Arial', 8), width=8,
                  command=self._show_all_animals).pack(side=tk.LEFT)
        self.lbl_excluded = tk.Label(anim_frame, text='',
                                      font=('Arial', 8), fg='#C62828')
        self.lbl_excluded.pack(anchor='w')

        tk.Frame(ctrl, height=6).pack()

        # ── Figure style ──────────────────────────────────────────────────────
        fig_frame2 = tk.LabelFrame(ctrl, text='Figure style',
                                   font=('Arial', 9), padx=6, pady=4)
        fig_frame2.pack(fill=tk.X, pady=4)
        font_row = tk.Frame(fig_frame2)
        font_row.pack(anchor='w')
        tk.Label(font_row, text='Font size:', font=('Arial', 9)).pack(side=tk.LEFT)
        tk.Spinbox(font_row, from_=7, to=22, increment=1,
                   textvariable=self._font_size, width=4, font=('Arial', 9),
                   command=self._on_font_change).pack(side=tk.LEFT, padx=4)
        tk.Button(fig_frame2, text='Export figure (300 dpi)',
                  font=('Arial', 8), command=self._export_figure
                  ).pack(anchor='w', pady=(3, 0))

        tk.Frame(ctrl, height=4).pack()

        # ── Run + results ─────────────────────────────────────────────────────
        tk.Button(ctrl, text='Run statistical tests',
                  command=self._run_model,
                  bg='#1B5E20', fg='white',
                  font=('Arial', 10, 'bold'),
                  width=28, height=2).pack(pady=4)
        res_btn_row = tk.Frame(ctrl)
        res_btn_row.pack(fill=tk.X)
        tk.Button(res_btn_row, text='View full results',
                  font=('Arial', 8), width=16,
                  command=self._show_results_window).pack(side=tk.LEFT, padx=(0, 4))
        tk.Button(res_btn_row, text='Export stats xlsx',
                  font=('Arial', 8), width=16,
                  command=self._export_stats).pack(side=tk.LEFT)
        self.lbl_model_summary = tk.Label(ctrl, text='',
                                           font=('Arial', 8), fg='#333',
                                           justify='left', wraplength=340)
        self.lbl_model_summary.pack(anchor='w', pady=(4, 0))

    # ── Dynamic condition UI rebuild ──────────────────────────────────────────

    def _rebuild_condition_ui(self):
        """
        Rebuild the condition list and bracket checkboxes whenever a new file
        is loaded.  Each condition row has:
          [☑ include]  [○ reference]  Condition label
        Unchecking "include" removes the condition from the plot and stats.
        """
        # ── Condition rows ────────────────────────────────────────────────────
        for w in self._ref_frame.winfo_children():
            w.destroy()

        if not self._conditions_order:
            tk.Label(self._ref_frame, text='(no conditions found)',
                     font=('Arial', 8), fg='#999').pack(anchor='w')
        else:
            # Default reference = the control (WT / Cre-) if present, else first
            if self.ref_var.get() not in self._conditions_order:
                ctrl = next((c for c in self._conditions_order if _is_control(c)),
                            self._conditions_order[0])
                self.ref_var.set(ctrl)

            # Column headers
            hdr = tk.Frame(self._ref_frame)
            hdr.pack(fill=tk.X, pady=(0, 2))
            tk.Label(hdr, text='Incl', font=('Arial', 7), fg='#888',
                     width=4).pack(side=tk.LEFT)
            tk.Label(hdr, text='Ref', font=('Arial', 7), fg='#888',
                     width=4).pack(side=tk.LEFT)
            tk.Label(hdr, text='Condition', font=('Arial', 7),
                     fg='#888').pack(side=tk.LEFT)
            tk.Label(hdr, text='order', font=('Arial', 7),
                     fg='#888').pack(side=tk.RIGHT, padx=(0, 4))

            ttk.Separator(self._ref_frame,
                          orient='horizontal').pack(fill=tk.X, pady=(0, 3))

            n_conds = len(self._conditions_order)
            for pos, cond in enumerate(self._conditions_order):
                color = self._condition_colors.get(cond, '#000')

                # Preserve existing include state across UI rebuilds
                if cond not in self._condition_include:
                    self._condition_include[cond] = tk.BooleanVar(value=True)
                inc_var = self._condition_include[cond]

                row = tk.Frame(self._ref_frame)
                row.pack(fill=tk.X, pady=1)

                # Include checkbox (no text — header labels it)
                tk.Checkbutton(
                    row, variable=inc_var,
                    command=self._on_include_change,
                    width=2
                ).pack(side=tk.LEFT)

                # Reference radio (no text — label follows)
                tk.Radiobutton(
                    row, variable=self.ref_var, value=cond,
                    width=1
                ).pack(side=tk.LEFT)

                # ▲/▼ move the group left/right on the plot (display order)
                tk.Button(row, text='▼', font=('Arial', 7), width=2, relief='flat',
                          state=tk.NORMAL if pos < n_conds - 1 else tk.DISABLED,
                          command=lambda c=cond: self._move_condition(c, +1)
                          ).pack(side=tk.RIGHT)
                tk.Button(row, text='▲', font=('Arial', 7), width=2, relief='flat',
                          state=tk.NORMAL if pos > 0 else tk.DISABLED,
                          command=lambda c=cond: self._move_condition(c, -1)
                          ).pack(side=tk.RIGHT)

                # Colored condition name label
                tk.Label(
                    row, text=cond, font=('Arial', 9), fg=color, anchor='w'
                ).pack(side=tk.LEFT, fill=tk.X, expand=True)

        # ── Bracket checkboxes ────────────────────────────────────────────────
        for w in self._brk_lf.winfo_children():
            w.destroy()

        self._all_pairs = list(combinations(self._conditions_order, 2))
        self._pair_vars = {}

        if not self._all_pairs:
            tk.Label(self._brk_lf, text='(load a file)',
                     font=('Arial', 8), fg='#999').pack(anchor='w')
        else:
            default_on = {(c1, c2)
                          for c1, c2 in self._all_pairs
                          if _is_control(c1) or _is_control(c2)}
            for c1, c2 in self._all_pairs:
                key = frozenset((c1, c2))
                var = tk.BooleanVar(value=self._pair_state.get(key, (c1, c2) in default_on))
                var.trace_add('write', lambda *a, k=key, v=var: self._pair_state.__setitem__(k, bool(v.get())))
                self._pair_vars[(c1, c2)] = var
                tk.Checkbutton(
                    self._brk_lf,
                    text=f'{c1}  ↔  {c2}',
                    variable=var,
                    font=('Arial', 8),
                    command=self._redraw_brackets
                ).pack(anchor='w')

    def _move_condition(self, cond, delta):
        """▲/▼: move a condition one place left/right in the display order.
        The order is remembered for the session (also across Add file(s))."""
        order = list(self._conditions_order)
        if cond not in order:
            return
        i = order.index(cond)
        j = i + delta
        if not (0 <= j < len(order)):
            return
        order[i], order[j] = order[j], order[i]
        self._conditions_order = order
        self._user_cond_order  = list(order)
        self._pair_results = {}          # pair keys are ordered tuples — recompute
        self._rebuild_condition_ui()
        self._update_plot()

    def _on_include_change(self):
        """Called when any include checkbox is toggled."""
        # If the current reference condition just got excluded, auto-pick a new one
        ref = self.ref_var.get()
        if ref in self._condition_include and not self._condition_include[ref].get():
            # Pick first included condition as new reference
            for c in self._conditions_order:
                if self._condition_include.get(c, tk.BooleanVar(value=True)).get():
                    self.ref_var.set(c)
                    break
        self._pair_results = {}   # invalidate model results
        self._update_plot()

    def _included_conditions(self):
        """Return the ordered list of conditions whose include checkbox is checked."""
        return [c for c in self._conditions_order
                if self._condition_include.get(c, tk.BooleanVar(value=True)).get()]

    # ── Animal ID editor ──────────────────────────────────────────────────────

    def _pipeline_animal_by_file(self):
        """{file: animal ID} as the pipeline sees it — the AnimalID column when
        present, else a filename guess. Ignores manual _animal_map overrides
        (those are shown/edited on top of this)."""
        df = self.df if self.df is not None else self._raw_df
        if df is None:
            return {}
        fcol = self.file_col or next(
            (c for c in df.columns if c.lower() == 'file'), df.columns[0])
        aid_col = next((c for c in df.columns
                        if str(c).strip().lower().replace(' ', '') == 'animalid'),
                       None)
        out = {}
        for fn, sub in df.groupby(fcol):
            aid = ''
            if aid_col is not None:
                v = sub[aid_col].iloc[0]
                if pd.notna(v) and str(v).strip():
                    aid = str(v).strip()
            out[str(fn)] = aid or _detect_animal(str(fn))
        return out

    def _pipeline_kocond_by_file(self):
        """{file: KO condition} from the prepared Condition column — the first
        non-control (Cre+) condition present in each file."""
        df = self.df
        if df is None or 'Condition' not in df.columns:
            return {}
        fcol = self.file_col
        out = {}
        for fn, sub in df.groupby(fcol):
            ko = sub.loc[~sub['Condition'].map(_is_control), 'Condition']
            out[str(fn)] = ko.iloc[0] if len(ko) else ''
        return out

    def _show_animal_editor(self):
        if self.df is None:
            messagebox.showinfo('No data', 'Load a file first.')
            return
        file_to_detected = self._pipeline_animal_by_file()
        from collections import defaultdict
        groups = defaultdict(list)
        for fn, aid in file_to_detected.items():
            groups[aid].append(fn)

        win = tk.Toplevel(self.root)
        win.title('Animal ID Editor')
        win.geometry('820x520')
        win.grab_set()

        tk.Label(win,
                 text='Animal IDs come from the AnimalID column defined in the '
                      'Soma Detector (filename guess if absent). Type an override '
                      'to correct any that are wrong, then click Apply.',
                 font=('Arial', 9), wraplength=800,
                 justify='left').pack(anchor='w', padx=10, pady=(8, 4))

        hdr = tk.Frame(win, bg='#CCCCCC')
        hdr.pack(fill=tk.X, padx=10)
        for txt, w in [('Animal ID (from data)', 22),
                       ('Sample Filename',     46),
                       ('Override (leave blank = keep)', 26)]:
            tk.Label(hdr, text=txt, font=('Arial', 9, 'bold'),
                     bg='#CCCCCC', width=w, anchor='w').pack(side=tk.LEFT, padx=2)

        container = tk.Frame(win)
        container.pack(fill=tk.BOTH, expand=True, padx=10, pady=4)
        canvas = tk.Canvas(container, highlightthickness=0)
        sb = ttk.Scrollbar(container, orient='vertical', command=canvas.yview)
        canvas.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        inner = tk.Frame(canvas)
        canvas.create_window((0, 0), window=inner, anchor='nw')
        inner.bind('<Configure>',
                   lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.bind('<MouseWheel>',
                    lambda e: canvas.yview_scroll(int(-1*(e.delta/120)), 'units'))

        override_vars = {}
        for detected_id in sorted(groups.keys()):
            sample = sorted(groups[detected_id])[0]
            current_override = self._animal_map.get(sample, '')
            row = tk.Frame(inner, bg='white' if len(override_vars) % 2 == 0
                           else '#F5F5F5')
            row.pack(fill=tk.X)
            tk.Label(row, text=detected_id, font=('Courier', 9),
                     width=22, anchor='w',
                     bg=row.cget('bg')).pack(side=tk.LEFT, padx=2)
            tk.Label(row, text=sample[:58], font=('Arial', 8),
                     width=46, anchor='w',
                     bg=row.cget('bg')).pack(side=tk.LEFT, padx=2)
            var = tk.StringVar(value=current_override)
            tk.Entry(row, textvariable=var,
                     font=('Arial', 9), width=26).pack(side=tk.LEFT, padx=2)
            override_vars[detected_id] = (var, groups[detected_id])

        def _apply():
            self._animal_map = {}
            for detected_id, (var, filenames) in override_vars.items():
                override = var.get().strip()
                if override and override != detected_id:
                    for fn in filenames:
                        self._animal_map[fn] = override
            self._prepare_data(self._raw_df)
            self._pair_results = {}
            self._update_plot()
            win.destroy()

        def _reset():
            for var, _ in override_vars.values():
                var.set('')

        btn_row = tk.Frame(win)
        btn_row.pack(pady=6)
        tk.Button(btn_row, text='Apply', command=_apply,
                  bg='#1B5E20', fg='white',
                  font=('Arial', 10, 'bold'), width=12).pack(side=tk.LEFT, padx=6)
        tk.Button(btn_row, text='Reset to Auto', command=_reset,
                  font=('Arial', 9), width=14).pack(side=tk.LEFT, padx=6)
        tk.Button(btn_row, text='Close', command=win.destroy,
                  font=('Arial', 9), width=10).pack(side=tk.LEFT, padx=6)

    # ── File loading ──────────────────────────────────────────────────────────

    def _open_file(self):
        """Open one or more files, replacing any currently loaded data."""
        paths = filedialog.askopenfilenames(
            title='Open classified soma measurements',
            filetypes=[('Excel files', '*.xlsx'), ('All files', '*.*')])
        if paths:
            self._load_files(list(paths), append=False)

    def _add_file(self):
        """Open one or more files and merge with currently loaded data."""
        paths = filedialog.askopenfilenames(
            title='Add classified soma measurements',
            filetypes=[('Excel files', '*.xlsx'), ('All files', '*.*')])
        if paths:
            self._load_files(list(paths), append=True)

    def _clear_files(self):
        self._loaded_paths  = []
        self._raw_df        = None
        self.df             = None
        self._animal_map    = {}
        self._excluded_animals   = set()
        self._highlighted_animal = None
        self._show_all_bright    = False
        self._pair_results  = {}
        self.lbl_file.config(text='No file loaded')
        self.lbl_n.config(text='')
        self.ax.cla()
        self.canvas.draw_idle()

    def _load_files(self, paths, append=False):
        try:
            new_dfs = [_load_excel(p) for p in paths]

            # Per-file mouse genotype (plain, mathtext) from the Genotype column
            # (rich text), used to label conditions instead of parsing filenames.
            geno_map = {}
            for p in paths:
                geno_map.update(_read_genotype_labels(p))

            if append and self._raw_df is not None:
                combined = pd.concat([self._raw_df] + new_dfs, ignore_index=True)
                self._loaded_paths.extend(paths)
                self._geno_forms_by_file = {**getattr(self, '_geno_forms_by_file', {}),
                                            **geno_map}
            else:
                combined = (pd.concat(new_dfs, ignore_index=True)
                            if len(new_dfs) > 1 else new_dfs[0])
                self._loaded_paths   = list(paths)
                self._animal_map     = {}
                self._excluded_animals = set()
                self._pair_results   = {}
                self._geno_forms_by_file = geno_map

            self.out_dir = os.path.dirname(self._loaded_paths[0])
            self._update_file_label()
            self._prepare_data(combined)
            self._update_plot()
        except Exception as e:
            messagebox.showerror('Error loading file(s)', str(e))

    def _update_file_label(self):
        n = len(self._loaded_paths)
        if n == 0:
            self.lbl_file.config(text='No file loaded')
        elif n == 1:
            self.lbl_file.config(text=os.path.basename(self._loaded_paths[0]))
        else:
            names = '  |  '.join(os.path.basename(p) for p in self._loaded_paths)
            self.lbl_file.config(text=f'{n} files: {names}')

    def _on_zrel_change(self, *_):
        self._pair_results = {}
        self._update_zrel_label()
        self._update_plot()

    def _on_cre_window_change(self, *_):
        if self._raw_df is None:
            return
        current_outcome = self.outcome_var.get()
        self._prepare_data(self._raw_df)
        if current_outcome in self._outcome_cols:
            self.outcome_var.set(current_outcome)
        self._pair_results = {}
        self._update_plot()

    # ── Grouping (auto + manual table override) ─────────────────────────────────
    #
    # _apply_grouping(df) sets df['Animal'] and df['Condition'] in place.
    # Default behaviour is the automatic filename-driven logic (_auto_grouping).
    # The "Grouping setup" table lets the user override Animal and the Cre+
    # Condition per file; Cre-negative cells always remain 'WT'.

    def _auto_grouping(self, df):
        """Automatic Animal/Condition derivation.

        Animal ID: prefer the explicit AnimalID column (user-defined in the soma
        detector); fall back to parsing the filename only when it is absent/blank.
        A manual "Edit Animal IDs" entry always overrides both.
        """
        aid_col = next((c for c in df.columns
                        if str(c).strip().lower().replace(' ', '') == 'animalid'),
                       None)

        def _animal(r):
            fn = str(r[self.file_col])
            if fn in self._animal_map:
                return self._animal_map[fn]
            if aid_col is not None:
                v = r[aid_col]
                if pd.notna(v) and str(v).strip():
                    return str(v).strip()
            return _detect_animal(fn)

        df['Animal'] = df.apply(_animal, axis=1)

        grp_col = next((c for c in df.columns if c.lower() == 'group'), None)
        _CRE_STATUS_ONLY = {
            'mcherry only', 'mcherry+only', 'mcherry-only',
            'gfp-cre only', 'gfp-cre+only', 'gfp-cre-only',
            'gfp only', 'gfp+only', 'gfp-only',
            'mcherry + gfp-cre', 'mcherry+gfp-cre',
            'gfp-cre + mcherry', 'gfp-cre+mcherry',
            'both positive', 'both negative',
        }
        group_is_cre_status_only = (
            grp_col is not None
            and df[grp_col].dropna().astype(str).str.strip().str.lower()
                  .isin(_CRE_STATUS_ONLY).all()
        )

        geno_map = getattr(self, '_geno_forms_by_file', {}) or {}

        def _cls_for_row(r):
            cls_val = r.get(self.cls_col) if self.cls_col else None
            if ((cls_val is None
                 or (isinstance(cls_val, float) and pd.isna(cls_val))
                 or str(cls_val).strip() == '') and grp_col):
                cls_val = r.get(grp_col)
            return cls_val

        # Display map: plain Condition → matplotlib label (mathtext for the
        # genotype so italics/superscript render on the plot).
        self._cond_display = {}

        # The Cre classifier exports each file's own Cre-channel name ('Cre
        # channel' column).  When a dataset mixes Cre constructs (e.g.
        # GFP-AKT1-Cre and GFP-AKT3-Cre images of the same genotype, each
        # over-expressing a different AKT) the construct is part of the
        # condition; with a single construct the plain 'Cre+' label is kept.
        cre_ch_col = next((c for c in df.columns
                           if str(c).strip().lower() == 'cre channel'), None)
        cre_ch_multi = (cre_ch_col is not None and
                        df[cre_ch_col].dropna().astype(str).str.strip()
                          .replace('', np.nan).dropna().nunique() > 1)

        def _cre_tag(r):
            if cre_ch_multi:
                cc = r.get(cre_ch_col)
                cc = '' if cc is None or (isinstance(cc, float) and pd.isna(cc)) else str(cc).strip()
                if cc:
                    return cc + '+'
            return 'Cre+'

        if geno_map and self.cls_col:
            # Condition = Cre status + mouse genotype:
            #   Cre-negative control → 'Cre-'
            #   Cre-positive KO      → '<genotype> Cre+'  (formatted on plots)
            #                       or '<genotype> <Cre channel>+' with mixed constructs
            # Genotypes are typed per file in the detector; spellings that differ
            # only in case or spacing ('AKT3^flx/flx' vs 'Akt3^flx/flx',
            # 'Akt1^flx/flx; Akt3^flx/flx') are the same condition — merge them
            # under the first spelling met.
            canon = {}

            def _cond_geno(r):
                cls_val = str(_cls_for_row(r) or '').strip().lower()
                if cls_val in _WT_PHENOTYPE_VALS:
                    self._cond_display['Cre-'] = 'Cre-'
                    return 'Cre-'
                tag = _cre_tag(r)
                plain, math = geno_map.get(str(r[self.file_col]), ('', ''))
                if plain:
                    cond = plain + ' ' + tag
                    key = re.sub(r'\s+', '', cond).casefold()
                    cond = canon.setdefault(key, cond)
                    self._cond_display.setdefault(cond, (math + ' ' + tag) if math else cond)
                    return cond
                self._cond_display[tag] = tag
                return tag
            df['Condition'] = df.apply(_cond_geno, axis=1)
        elif grp_col and not group_is_cre_status_only:
            df['Condition'] = df[grp_col].apply(lambda g: str(g).strip())
        elif self.cls_col:
            df['Condition'] = df.apply(
                lambda r: _condition_from_filename(r[self.file_col],
                                                   _cls_for_row(r)), axis=1)
        elif grp_col:
            df['Condition'] = df.apply(
                lambda r: _condition_from_filename(
                    r[self.file_col], r.get(grp_col)), axis=1)
        else:
            df['Condition'] = df[self.file_col].apply(_condition_from_filename)

    def _apply_grouping(self, df):
        # Auto first; manual entries then override on top.
        self._auto_grouping(df)
        mapping = (self._grouping_cfg or {}).get('map', {})
        if not mapping:
            return

        def animal_row(r):
            e = mapping.get(str(r[self.file_col]))
            return e['animal'] if (e and e.get('animal')) else r['Animal']

        def cond_row(r):
            e = mapping.get(str(r[self.file_col]))
            # Override only KO cells; the control (WT / Cre-) group stays put.
            if e and e.get('condition') and not _is_control(r['Condition']):
                return e['condition']
            return r['Condition']

        df['Animal']    = df.apply(animal_row, axis=1)
        df['Condition'] = df.apply(cond_row, axis=1)

    def _reapply_grouping(self):
        """Re-run grouping + downstream prep after a config change."""
        if self._raw_df is None:
            return
        current_outcome = self.outcome_var.get()
        self._prepare_data(self._raw_df)
        if current_outcome in self._outcome_cols:
            self.outcome_var.set(current_outcome)
        self._pair_results = {}
        self._update_plot()

    def _open_grouping_setup(self):
        if self._raw_df is None:
            messagebox.showinfo('Grouping setup', 'Load a file first.')
            return
        file_col = self.file_col or next(
            (c for c in self._raw_df.columns if c.lower() == 'file'),
            self._raw_df.columns[0])
        files = list(dict.fromkeys(self._raw_df[file_col].astype(str).tolist()))

        # Prefill animal + KO-condition from the pipeline (AnimalID column and
        # the genotype/Cre-derived Condition), falling back to a filename guess.
        existing    = (self._grouping_cfg or {}).get('map', {})
        pipe_animal = self._pipeline_animal_by_file()
        pipe_cond   = self._pipeline_kocond_by_file()
        prefill = {}
        for fn in files:
            e = existing.get(fn, {})
            animal = e.get('animal') or pipe_animal.get(fn) or _detect_animal(fn)
            cond   = e.get('condition') or pipe_cond.get(fn)
            if not cond:
                c = _condition_from_filename(fn, None)
                cond = c if not _is_control(c) else _genotype_label(
                    _parse_filename_parts(fn)[0], fn)
            prefill[fn] = {'animal': animal, 'condition': cond}

        win = tk.Toplevel(self.root)
        win.title('Manual Grouping Table')
        win.geometry('1040x700')
        win.grab_set()

        tk.Label(win, text='Animal ID and Condition are prefilled from the data '
                 '(AnimalID column + genotype/Cre status). Edit to override. '
                 'The Condition applies to Cre+ cells; Cre− cells stay the '
                 'control group. The File cell is selectable — highlight any part '
                 'and Ctrl+C, then paste into "Files containing" to bulk-assign.',
                 font=('Arial', 9), wraplength=1010, justify='left'
                 ).pack(anchor='w', padx=10, pady=(8, 4))

        # Column-width sliders
        widths = tk.LabelFrame(win, text='Column widths', font=('Arial', 9),
                               padx=6, pady=2)
        widths.pack(fill=tk.X, padx=10, pady=2)
        file_w   = tk.IntVar(value=56)
        animal_w = tk.IntVar(value=28)
        cond_w   = tk.IntVar(value=28)
        file_entries, animal_entries, cond_combos = [], [], []
        hdr_labels = {}

        def update_widths(*_):
            for e in file_entries:
                e.config(width=file_w.get())
            for e in animal_entries:
                e.config(width=animal_w.get())
            for c in cond_combos:
                c.config(width=cond_w.get())
            if hdr_labels:
                hdr_labels['file'].config(width=file_w.get())
                hdr_labels['animal'].config(width=animal_w.get())
                hdr_labels['cond'].config(width=cond_w.get())

        for lbl, var, lo, hi in [('File', file_w, 20, 110),
                                 ('Animal ID', animal_w, 10, 50),
                                 ('Condition', cond_w, 10, 50)]:
            cell = tk.Frame(widths); cell.pack(side=tk.LEFT, padx=8)
            tk.Label(cell, text=lbl, font=('Arial', 8)).pack(side=tk.LEFT)
            tk.Scale(cell, from_=lo, to=hi, orient=tk.HORIZONTAL, variable=var,
                     length=150, command=update_widths).pack(side=tk.LEFT)

        # Bulk-assign bar
        bulk = tk.LabelFrame(win, text='Bulk assign', font=('Arial', 9),
                             padx=6, pady=4)
        bulk.pack(fill=tk.X, padx=10, pady=4)
        tk.Label(bulk, text='Files containing:', font=('Arial', 9)).pack(side=tk.LEFT)
        sub_var = tk.StringVar()
        tk.Entry(bulk, textvariable=sub_var, width=28).pack(side=tk.LEFT, padx=3)
        tk.Label(bulk, text='→ set Condition:', font=('Arial', 9)).pack(side=tk.LEFT)
        bulk_cond = tk.StringVar()
        tk.Entry(bulk, textvariable=bulk_cond, width=24).pack(side=tk.LEFT, padx=3)

        row_vars = {}   # fn -> (animal_var, cond_var)

        def do_bulk():
            s = sub_var.get().strip().lower()
            val = bulk_cond.get().strip()
            if not s or not val:
                return
            for fn, (_av, cv) in row_vars.items():
                if s in os.path.basename(fn).lower():
                    cv.set(val)
        tk.Button(bulk, text='Apply to matching', command=do_bulk,
                  font=('Arial', 9), bg='#1976D2', fg='white').pack(side=tk.LEFT, padx=6)

        # Scrollable table
        table_wrap = tk.Frame(win)
        table_wrap.pack(fill=tk.BOTH, expand=True, padx=10, pady=4)
        canvas = tk.Canvas(table_wrap, highlightthickness=0)
        vsb    = ttk.Scrollbar(table_wrap, orient='vertical', command=canvas.yview)
        hsb    = ttk.Scrollbar(table_wrap, orient='horizontal', command=canvas.xview)
        inner  = tk.Frame(canvas)
        inner.bind('<Configure>',
                   lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=inner, anchor='nw')
        canvas.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        hsb.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.bind_all('<MouseWheel>',
                        lambda e: canvas.yview_scroll(int(-e.delta / 120), 'units'))

        # Header
        hdr = tk.Frame(inner, bg='#CCCCCC'); hdr.grid(row=0, column=0, sticky='w')
        hdr_labels['file'] = tk.Label(hdr, text='File', font=('Arial', 9, 'bold'),
                                      bg='#CCCCCC', width=file_w.get(), anchor='w')
        hdr_labels['file'].pack(side=tk.LEFT, padx=1)
        hdr_labels['animal'] = tk.Label(hdr, text='Animal ID', font=('Arial', 9, 'bold'),
                                        bg='#CCCCCC', width=animal_w.get(), anchor='w')
        hdr_labels['animal'].pack(side=tk.LEFT, padx=1)
        hdr_labels['cond'] = tk.Label(hdr, text='Condition (Cre+)',
                                      font=('Arial', 9, 'bold'),
                                      bg='#CCCCCC', width=cond_w.get(), anchor='w')
        hdr_labels['cond'].pack(side=tk.LEFT, padx=1)

        cond_values = sorted({prefill[f]['condition'] for f in files})
        for i, fn in enumerate(files, 1):
            rowf = tk.Frame(inner); rowf.grid(row=i, column=0, sticky='w')
            # File cell: read-only Entry so the text can be selected and copied
            fe = tk.Entry(rowf, font=('Consolas', 8), width=file_w.get())
            fe.insert(0, os.path.basename(fn))
            fe.config(state='readonly')
            fe.pack(side=tk.LEFT, padx=1)
            file_entries.append(fe)
            av = tk.StringVar(value=prefill[fn]['animal'])
            cv = tk.StringVar(value=prefill[fn]['condition'])
            ae = tk.Entry(rowf, textvariable=av, width=animal_w.get())
            ae.pack(side=tk.LEFT, padx=1)
            animal_entries.append(ae)
            cc = ttk.Combobox(rowf, textvariable=cv, values=cond_values,
                              width=cond_w.get())
            cc.pack(side=tk.LEFT, padx=1)
            cond_combos.append(cc)
            row_vars[fn] = (av, cv)

        # Buttons
        btns = tk.Frame(win); btns.pack(fill=tk.X, padx=10, pady=8)

        def collect():
            return {'map': {fn: {'animal': av.get().strip(),
                                 'condition': cv.get().strip()}
                            for fn, (av, cv) in row_vars.items()}}

        def do_apply():
            self._grouping_cfg = collect()
            self._reapply_grouping()
            win.destroy()

        def do_save():
            p = filedialog.asksaveasfilename(
                title='Save grouping map', defaultextension='.json',
                filetypes=[('JSON', '*.json')])
            if p:
                with open(p, 'w', encoding='utf-8') as f:
                    json.dump(collect(), f, indent=2)

        def do_load():
            p = filedialog.askopenfilename(
                title='Load grouping map', filetypes=[('JSON', '*.json')])
            if p:
                with open(p, encoding='utf-8') as f:
                    loaded = json.load(f)
                m = loaded.get('map', {})
                for fn, (av, cv) in row_vars.items():
                    if fn in m:
                        av.set(m[fn].get('animal', av.get()))
                        cv.set(m[fn].get('condition', cv.get()))

        tk.Button(btns, text='Apply', command=do_apply, bg='#2E7D32',
                  fg='white', font=('Arial', 9, 'bold'), width=10).pack(side=tk.LEFT)
        tk.Button(btns, text='Save map…', command=do_save,
                  font=('Arial', 9), width=12).pack(side=tk.LEFT, padx=4)
        tk.Button(btns, text='Load map…', command=do_load,
                  font=('Arial', 9), width=12).pack(side=tk.LEFT, padx=2)
        tk.Button(btns, text='Close', command=win.destroy,
                  font=('Arial', 9), width=8).pack(side=tk.RIGHT)

    # ── Animal explorer ───────────────────────────────────────────────────────

    def _populate_animal_list(self):
        if self.df is None:
            return
        self._animal_lb.delete(0, tk.END)
        cond_order = {c: i for i, c in enumerate(self._conditions_order)}
        rows = (self.df.groupby(['Animal', 'Condition'])
                .size().reset_index(name='n'))
        rows['_ord'] = rows['Condition'].map(cond_order).fillna(99)
        rows = rows.sort_values(['_ord', 'Animal'])
        for _, r in rows.iterrows():
            animal = r['Animal']
            prefix = '✕ ' if animal in self._excluded_animals else '● '
            label  = f"{prefix}{animal:<18}  {r['Condition']:<24}  {int(r['n'])}c"
            self._animal_lb.insert(tk.END, label)
            if animal in self._excluded_animals:
                self._animal_lb.itemconfig(tk.END, fg='#C62828')
            else:
                col = self._condition_colors.get(r['Condition'], '#333')
                self._animal_lb.itemconfig(tk.END, fg=col)
        self._update_excluded_label()

    def _selected_animal(self):
        sel = self._animal_lb.curselection()
        if not sel:
            return None
        text = self._animal_lb.get(sel[0])
        return text[2:].split('  ')[0].strip()

    def _on_animal_select(self, *_):
        self._highlighted_animal = self._selected_animal()
        self._show_all_bright    = False
        self._draw_plot()

    def _toggle_animal_exclude(self):
        animal = self._selected_animal()
        if animal is None:
            return
        if animal in self._excluded_animals:
            self._excluded_animals.discard(animal)
        else:
            self._excluded_animals.add(animal)
        self._pair_results = {}
        self._populate_animal_list()
        self._update_zrel_label()
        self._draw_plot()

    def _clear_exclusions(self):
        self._excluded_animals.clear()
        self._pair_results = {}
        self._populate_animal_list()
        self._update_zrel_label()
        self._draw_plot()

    def _show_all_animals(self):
        self._highlighted_animal = None
        self._show_all_bright    = True
        self._animal_lb.selection_clear(0, tk.END)
        self._draw_plot()

    def _update_excluded_label(self):
        n = len(self._excluded_animals)
        self.lbl_excluded.config(
            text=f'{n} animal{"s" if n != 1 else ""} excluded' if n else '')

    def _update_cre_label(self):
        """Coverage of the Cre⁻ normalisation, counted on the SAME cells the
        plots use (depth filter + animal exclusions applied).  F0 itself is
        computed at load from Cre⁻ cells at matching depth in the same file,
        regardless of the depth filter — a kept cell near the cutoff is still
        referenced against Cre⁻ cells at its own depth."""
        if self.df is None:
            self.lbl_cre_cov.config(text='')
            return
        norm_col = next((c for c in self.df.columns if c.endswith('/ Cre-')), None)
        if norm_col is None:
            self.lbl_cre_cov.config(text='No Cre⁻ column found')
            return
        df_f    = self._filtered_df()
        exp     = df_f[~df_f['Condition'].map(_is_control)]
        total   = len(exp)
        covered = int(exp[norm_col].notna().sum())
        n_ref   = int(self.df[self.cls_col].apply(
            lambda v: str(v).strip().lower() in _WT_PHENOTYPE_VALS).sum()) \
            if self.cls_col else 0
        try:
            win = int(self._cre_window.get())
        except Exception:
            win = 0
        if total > 0:
            self.lbl_cre_cov.config(
                text=f'Normalisation: {covered} of {total} kept Cre⁺ cells have a '
                     f'Cre⁻ reference within ±{win} slices ({100*covered/total:.0f}%)\n'
                     f'F0 from {n_ref} Cre⁻ cells at matching depth, all depths')
        else:
            self.lbl_cre_cov.config(text='Normalisation: no Cre⁺ cells in the kept set')

    def _update_zrel_label(self):
        if self.df is None:
            return
        df_f   = self._filtered_df()
        counts = df_f['Condition'].value_counts()
        animals = df_f.groupby('Condition')['Animal'].nunique()
        parts = [f"{c}: {animals.get(c,0)}a/{counts.get(c,0)}c"
                 for c in self._conditions_order if c in counts]
        self.lbl_n.config(text='   '.join(parts))
        if self._zrel_enable.get():
            total_kept = len(df_f)
            total_all  = len(self.df)
            self.lbl_zrel_n.config(
                text=f'Depth filter: keeping {total_kept} of {total_all} cells '
                     f'({100*total_kept/max(total_all,1):.0f}%), all groups')
        else:
            self.lbl_zrel_n.config(text='')
        # The normalisation coverage is reported on the same kept set.
        self._update_cre_label()

    def _filtered_df(self):
        df = self.df
        if df is None:
            return df
        if self._zrel_enable.get() and 'z_rel' in df.columns:
            try:
                cutoff = float(self._zrel_cutoff.get())
            except Exception:
                cutoff = 1.0
            df = df[df['z_rel'] <= cutoff]
        if self._excluded_animals:
            df = df[~df['Animal'].isin(self._excluded_animals)]
        return df

    def _prepare_data(self, raw):
        df = raw.copy()
        self._raw_df = raw

        self.file_col = next(
            (c for c in df.columns if c.lower() == 'file'), df.columns[0])
        self.cls_col  = next(
            (c for c in df.columns
             if 'class' in c.lower() or c.lower() == 'phenotype'), None)

        # Derive Animal and Condition. The manual grouping map (if configured
        # via "Grouping setup") overrides the automatic filename-driven labels.
        self._apply_grouping(df)

        # Relative Z-position
        slice_col = next((c for c in df.columns if c.lower() == 'slice'), None)
        if slice_col:
            file_max = df.groupby(self.file_col)[slice_col].transform('max')
            file_min = df.groupby(self.file_col)[slice_col].transform('min')
            span     = (file_max - file_min).clip(lower=1)
            df['z_rel'] = (df[slice_col] - file_min) / span
        else:
            df['z_rel'] = np.nan

        # Per-slice Cre- normalisation of the measured channel: each cell's
        # signal is divided by the mean of the Cre-negative cells within ±N
        # slices in the same file.
        _cre_neg_mask = (
            df[self.cls_col].apply(
                lambda v: str(v).strip().lower() in _WT_PHENOTYPE_VALS)
            if self.cls_col else pd.Series(False, index=df.index))
        if self.cls_col and _cre_neg_mask.any():
            # Normalize the MEASURED channel to the per-slice Cre-negative mean.
            # The measured channel is the last channel group in the export
            # (channel order is always detect-1, detect-2, measure). Keying off
            # the channel position rather than a 'p<marker>' name convention
            # fixes two failures: non-phospho antibodies (e.g. CaV2.1) were
            # missed, and capitalized detection-channel names (e.g.
            # PKCKW-GFP-Cre) were wrongly matched under case-insensitive 'p'.
            mean_cols    = [c for c in df.columns if str(c).startswith('Mean ')]
            measure_name = (str(mean_cols[-1])[len('Mean '):].strip()
                            if mean_cols else None)
            norm_cols = []
            if measure_name:
                for _pref in ('Mean', 'F/F0'):
                    _col = '{} {}'.format(_pref, measure_name)
                    if _col in df.columns and pd.api.types.is_numeric_dtype(df[_col]):
                        norm_cols.append(_col)

            slice_col = next((c for c in df.columns
                              if c.lower() == 'slice'), None)
            if norm_cols and slice_col:
                win     = int(self._cre_window.get())
                cre_neg = (df[_cre_neg_mask]
                           .rename(columns={slice_col: '__cre_slice'}))
                cell_keys = df[[self.file_col, slice_col]].drop_duplicates()
                windowed  = cell_keys.merge(
                    cre_neg[[self.file_col, '__cre_slice'] + norm_cols],
                    on=self.file_col)
                windowed = windowed[
                    (windowed['__cre_slice'] - windowed[slice_col]).abs() <= win]
                f0_win = (windowed.groupby([self.file_col, slice_col])[norm_cols]
                          .mean().add_suffix('__f0'))
                df = df.merge(f0_win.reset_index(),
                              on=[self.file_col, slice_col], how='left')
                for col in norm_cols:
                    df[f'{col} / Cre-'] = df[col] / df[f'{col}__f0']
                df = df.drop(columns=[c for c in df.columns
                                      if c.endswith('__f0')])

        self.df = df

        # ── Build dynamic conditions order and colors ──────────────────────────
        unique_conds = df['Condition'].dropna().unique().tolist()
        # Control (WT / Cre-) first, then KO conditions sorted alphabetically
        wt_conds  = sorted([c for c in unique_conds if _is_control(c)])
        ko_conds  = sorted([c for c in unique_conds
                            if not _is_control(c) and c != 'Unknown'])
        unk_conds = [c for c in unique_conds if c == 'Unknown']
        default_order = wt_conds + ko_conds + unk_conds
        # Keep the order the user chose with ▲/▼ (new conditions go to the end).
        self._conditions_order = ([c for c in self._user_cond_order if c in default_order]
                                  + [c for c in default_order if c not in self._user_cond_order])

        self._condition_colors = {}
        ko_idx = 0
        for cond in self._conditions_order:
            if _is_control(cond):
                self._condition_colors[cond] = _WT_COLOR
            elif cond == 'Unknown':
                self._condition_colors[cond] = '#CCCCCC'
            else:
                self._condition_colors[cond] = _KO_PALETTE[ko_idx % len(_KO_PALETTE)]
                ko_idx += 1

        # Rebuild condition-dependent UI
        self._rebuild_condition_ui()

        # Outcome columns
        skip = ['soma', 'slice', 'file', 'channel', 'source',
                'centroid', 'class', 'group', 'animal', 'phenotype',
                'condition', 'bg ', 'z_rel']
        self._outcome_cols = [
            c for c in df.columns
            if pd.api.types.is_numeric_dtype(df[c])
            and not any(p in c.lower() for p in skip)
        ]
        self.outcome_cb['values'] = self._outcome_cols
        default = next((c for c in self._outcome_cols if 'Area' in c),
                       self._outcome_cols[0] if self._outcome_cols else '')
        self.outcome_var.set(default)

        self._update_zrel_label()
        self._update_cre_label()
        self._populate_animal_list()

    # ── Plotting ──────────────────────────────────────────────────────────────

    def _update_plot(self, *_):
        if self.df is None:
            return
        self._draw_plot()

    def _redraw_brackets(self, *_):
        if self.df is None or not self._pair_results:
            return
        self._draw_plot()

    def _draw_plot(self):
        df      = self._filtered_df()
        out_col = self.outcome_var.get()
        if out_col not in df.columns:
            return
        # Reset click-to-identify map; it is rebuilt as scatter calls run.
        self._picker_data = []

        included = set(self._included_conditions())
        present  = [c for c in self._conditions_order
                    if c in df['Condition'].values and c in included]
        x_pos    = {c: i for i, c in enumerate(present)}

        self.ax.cla()
        rng   = np.random.default_rng(42)
        y_all = []
        hi_animal = self._highlighted_animal

        for cond in present:
            xi  = x_pos[cond]
            col = self._condition_colors.get(cond, '#888')
            sub = df[df['Condition'] == cond]

            if self._show_all_bright:
                sub_hi = sub
                sub_lo = sub.iloc[0:0]
            elif hi_animal:
                sub_hi = sub[sub['Animal'] == hi_animal]
                sub_lo = sub[sub['Animal'] != hi_animal]
            else:
                sub_hi = sub.iloc[0:0]
                sub_lo = sub

            vals_lo = pd.to_numeric(sub_lo[out_col], errors='coerce').dropna()
            vals_hi = pd.to_numeric(sub_hi[out_col], errors='coerce').dropna()

            if len(vals_lo) == 0 and len(vals_hi) == 0:
                continue

            y_all.extend(vals_lo.values)
            y_all.extend(vals_hi.values)

            if len(vals_lo):
                alpha_lo  = 0.18 if hi_animal else 0.55
                jitter_lo = rng.uniform(-0.30, 0.30, len(vals_lo))
                coll_lo = self.ax.scatter(xi + jitter_lo, vals_lo.values,
                                c=col, s=22, alpha=alpha_lo,
                                linewidths=0, zorder=2)
                # Click-to-identify: 6-pixel pick tolerance; remember the
                # DataFrame row indices in the same order as the plotted
                # points so _on_pick can resolve a click back to a cell.
                coll_lo.set_picker(6)
                self._picker_data.append((coll_lo, list(vals_lo.index)))

            if len(vals_hi):
                jitter_hi = rng.uniform(-0.30, 0.30, len(vals_hi))
                coll_hi = self.ax.scatter(xi + jitter_hi, vals_hi.values,
                                c=col, s=55, alpha=0.95,
                                edgecolors='black', linewidths=0.8, zorder=5)
                coll_hi.set_picker(6)
                self._picker_data.append((coll_hi, list(vals_hi.index)))

            all_vals = pd.to_numeric(sub[out_col], errors='coerce').dropna()
            if len(all_vals):
                med = all_vals.median()
                self.ax.plot([xi - 0.28, xi + 0.28], [med, med],
                             color='black', lw=2.0, solid_capstyle='round',
                             zorder=4)

        if self._pair_results:
            selected = []
            for (c1, c2), var in self._pair_vars.items():
                if var.get():
                    key = (c1, c2) if (c1, c2) in self._pair_results else (c2, c1)
                    if key in self._pair_results and c1 in x_pos and c2 in x_pos:
                        selected.append((c1, c2, self._pair_results[key]))
            if selected and y_all:
                y_max   = max(y_all)
                y_range = y_max - min(y_all)
                _draw_brackets(self.ax, selected, x_pos, y_max, y_range,
                               base_fs=int(self._font_size.get()))

        animals_per = (df[df['Condition'].isin(present)]
                       .groupby('Condition')['Animal'].nunique())
        # n per condition; genotype rendered via mathtext (italic/superscript).
        disp = getattr(self, '_cond_display', {})
        xlabels  = [f'{disp.get(c, c)}\n(n={animals_per.get(c, 0)})'
                    for c in present]
        fs       = int(self._font_size.get())
        tick_fs  = max(7, fs - 1)
        # Rotate if any label is long or there are many groups
        max_len  = max((len(c) for c in present), default=0)
        rotate   = max_len > 12 or len(present) > 5
        self.ax.set_xticks(range(len(present)))
        if rotate:
            self.ax.set_xticklabels(xlabels, fontsize=tick_fs,
                                    rotation=35, ha='right',
                                    rotation_mode='anchor')
        else:
            self.ax.set_xticklabels(xlabels, fontsize=tick_fs)
        self.ax.set_ylabel(out_col, fontsize=fs)
        self.ax.set_title(out_col, fontsize=fs + 1, fontweight='bold')
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.tick_params(axis='y', labelsize=fs)
        # Tighten bottom margin dynamically so rotated labels don't clip
        self.fig.subplots_adjust(
            bottom=0.30 if rotate else 0.14,
            left=0.13, right=0.96, top=0.90)

        if y_all and self._pair_results:
            cur_top = self.ax.get_ylim()[1]
            self.ax.set_ylim(bottom=self.ax.get_ylim()[0], top=cur_top * 1.02)

        self.canvas.draw_idle()

    # ── Click-to-identify ─────────────────────────────────────────────────────

    def _on_pick(self, event):
        """Resolve a pick event on a scatter point back to a DataFrame row
        and show the cell's identifying info in a popup."""
        if self.df is None or not getattr(event, 'ind', None):
            return
        # Find which scatter collection was hit and look up its row indices.
        df_indices = None
        for coll, indices in self._picker_data:
            if coll is event.artist:
                df_indices = indices
                break
        if df_indices is None:
            return
        idx_in_coll = int(event.ind[0])
        if idx_in_coll < 0 or idx_in_coll >= len(df_indices):
            return
        df_row_idx = df_indices[idx_in_coll]
        if df_row_idx not in self.df.index:
            return
        self._show_cell_info(self.df.loc[df_row_idx])

    def _show_cell_info(self, row):
        """Show a popup with cell identifiers and intensities, with
        copy-to-clipboard buttons for cross-referencing in other GUIs."""
        # Re-use a single popup window — open a second click replaces the
        # contents rather than spawning a stack of identical windows.
        if self._cell_info_win and self._cell_info_win.winfo_exists():
            self._cell_info_win.destroy()

        win = tk.Toplevel(self.root)
        win.title('Cell Info  —  click any point on the plot')
        win.geometry('640x460')
        self._cell_info_win = win

        txt = scrolledtext.ScrolledText(win, font=('Courier', 9),
                                         bg='#F5F5F5', wrap='none')
        txt.pack(fill=tk.BOTH, expand=True, padx=8, pady=(8, 4))

        # Build a structured, copy-friendly summary of the cell.
        lines = []
        lines.append('Identifiers')
        lines.append('-' * 56)
        for label, key in [
            ('File',       self.file_col),
            ('Animal',     'Animal'),
            ('Condition',  'Condition'),
            ('Class',      self.cls_col),
            ('Soma #',     'Soma #'),
            ('Slice',      'Slice'),
        ]:
            if key and key in row.index:
                lines.append(f'  {label:<11}: {row[key]}')

        lines.append('')
        lines.append('Geometry')
        lines.append('-' * 56)
        for c in row.index:
            cl = str(c).lower()
            if 'area' in cl or 'centroid' in cl or cl == 'z_rel':
                v = row[c]
                if isinstance(v, float):
                    lines.append(f'  {c:<22}: {v:.3f}')
                else:
                    lines.append(f'  {c:<22}: {v}')

        lines.append('')
        lines.append('Intensities')
        lines.append('-' * 56)
        for c in row.index:
            sc = str(c)
            if any(sc.startswith(p) for p in ('Mean ', 'BG ', 'F/F0 ', 'Range ')) \
                    or sc.endswith('/ Cre-'):
                v = row[c]
                if isinstance(v, float):
                    lines.append(f'  {sc:<26}: {v:.4f}')
                else:
                    lines.append(f'  {sc:<26}: {v}')

        text = '\n'.join(lines)
        txt.insert('1.0', text)
        txt.configure(state='disabled')

        # ── Buttons ───────────────────────────────────────────────────────────
        btn_row = tk.Frame(win)
        btn_row.pack(pady=(0, 8))

        file_val  = row.get(self.file_col, '')
        soma_val  = row.get('Soma #', '')
        slice_val = row.get('Slice', '')

        # Pull centroid (in µm) for the IPC payload — column header uses
        # the µ character which may be encoded different ways across files,
        # so match by 'centroid x' / 'centroid y' substring instead.
        cx_um = cy_um = None
        for c in row.index:
            cl = str(c).lower()
            if 'centroid x' in cl and cx_um is None:
                try: cx_um = float(row[c])
                except Exception: pass
            elif 'centroid y' in cl and cy_um is None:
                try: cy_um = float(row[c])
                except Exception: pass

        # Short identifier (best for pasting into the soma analysis GUI).
        ident_short = f'{file_val}\tSoma {soma_val}\tSlice {slice_val}'

        def _copy(payload):
            win.clipboard_clear()
            win.clipboard_append(payload)
            win.update()  # ensure clipboard persists after window closes

        def _open_in_soma_gui():
            """Send a jump command to the soma detector. Fall back to
            clipboard copy when no server is listening."""
            try:
                slice_val_int = int(slice_val)
            except Exception:
                slice_val_int = None
            search_dirs = []
            if self.out_dir:
                search_dirs.append(self.out_dir)
            for p in (self._loaded_paths or []):
                d = os.path.dirname(p)
                if d and d not in search_dirs:
                    search_dirs.append(d)
            payload = {
                'action':        'open_cell',
                'file':          str(file_val),
                'search_dirs':   search_dirs,
                'slice':         slice_val_int,
                'soma_num':      int(soma_val) if pd.notna(soma_val) else None,
                'centroid_x_um': cx_um,
                'centroid_y_um': cy_um,
                'condition':     str(row.get('Condition', '')),
            }
            if cell_jump_ipc.send_jump(payload):
                # Success — close the popup so focus shifts to the soma GUI.
                win.destroy()
                return

            # Not running — try to launch it, then poll until its IPC server
            # is up (a fresh launch takes a few seconds), and jump once it is.
            if not _launch_soma_detector():
                _copy(ident_short)
                messagebox.showinfo(
                    'Soma detector not found',
                    'Could not reach or launch the Soma Detector.\n\n'
                    'The identifier was copied to the clipboard. Open the Soma '
                    'Detector manually and click "Open in soma detector" again.')
                return

            def _retry(n=0):
                if not (win.winfo_exists()):
                    return
                if cell_jump_ipc.send_jump(payload):
                    win.destroy()
                    return
                if n < 40:                      # ~40 x 300 ms = 12 s
                    self.root.after(300, lambda: _retry(n + 1))
                else:
                    _copy(ident_short)
                    messagebox.showinfo(
                        'Soma detector starting',
                        'The Soma Detector is opening but did not respond in '
                        'time. Once it has finished loading, click "Open in '
                        'soma detector" again.\n\nThe identifier was copied to '
                        'the clipboard as a backup.')
            self.root.after(800, _retry)   # give the app a moment to start

        tk.Button(btn_row, text='Open in soma detector',
                  font=('Arial', 9, 'bold'),
                  bg='#1B5E20', fg='white', width=22,
                  command=_open_in_soma_gui
                  ).pack(side=tk.LEFT, padx=4)
        tk.Button(btn_row, text='Copy identifier',
                  font=('Arial', 9), width=14,
                  command=lambda: _copy(ident_short)
                  ).pack(side=tk.LEFT, padx=4)
        tk.Button(btn_row, text='Copy all info',
                  font=('Arial', 9), width=12,
                  command=lambda: _copy(text)
                  ).pack(side=tk.LEFT, padx=4)
        tk.Button(btn_row, text='Close',
                  font=('Arial', 9), width=8,
                  command=win.destroy
                  ).pack(side=tk.LEFT, padx=4)

    # ── Statistical tests ─────────────────────────────────────────────────────

    def _run_model(self):
        """Run the group comparisons on the current dataset.

        Test of record — the clustered-data model of Moen, Fricano-Kugler,
        Luikart & O'Malley (2016, PLoS One 11:e0146721), as used in Prina,
        Goyette et al. (2026, Brain): a linear mixed-effects model with a
        random intercept for each animal (cells nested within animals),
        fitted by maximum likelihood, with standard errors cluster-robust
        on animal and Scheffé-adjusted pairwise contrasts (Stata: xtmixed
        … || animalid:, vce(robust); pwcompare, mcompare(scheffe)).  It is
        implemented in mixed_stats.py and drives the significance brackets.

        Descriptive secondary tests, Bonferroni-corrected across the
        comparisons:
          * p_animal — Welch's t-test on per-animal means
          * p_cell   — Welch's t-test on cell-level values (ignores the
            nesting of cells within animals; over-states significance)
          * p_MW     — Mann-Whitney U on cell-level values
        Only if the mixed model cannot be fitted does the animal-level
        adjusted p-value drive the brackets (cell-level when a group has
        fewer than two animals).
        """
        if self.df is None:
            messagebox.showinfo('No data', 'Load a file first.')
            return

        df       = self._filtered_df().copy()
        out_col  = self.outcome_var.get()
        ref_cond = self.ref_var.get()

        df['_y'] = pd.to_numeric(df[out_col], errors='coerce')
        df = df.dropna(subset=['_y', 'Condition', 'Animal'])
        df = df[df['Condition'] != 'Unknown']

        included = set(self._included_conditions())
        df       = df[df['Condition'].isin(included)]
        present  = [c for c in self._conditions_order
                    if c in df['Condition'].values and c in included]
        if len(present) < 2:
            messagebox.showwarning('Warning', 'Need at least 2 conditions.')
            return
        if ref_cond not in present:
            ref_cond = present[0]

        # Aggregate cells to one mean per animal for the animal-level test.
        animal_means = (df.groupby(['Condition', 'Animal'])['_y']
                          .mean().reset_index())

        lines = []
        sep   = '=' * 78
        lines += [
            sep,
            'PAIRWISE STATISTICAL TESTS  (test of record: the mixed-effects model, below)',
            '',
            f'  Outcome   : {out_col}',
            f'  Reference : {ref_cond}',
            f'  N cells   : {len(df):,}',
            f'  N animals : {df["Animal"].nunique()}',
            sep,
            '',
            'Condition summary:',
            f'{"Condition":<32} {"n_anim":>7} {"n_cells":>8} '
            f'{"mean(cells)":>12} {"SD(cells)":>10}',
            '-' * 78,
        ]
        for cond in present:
            sub_cells   = df[df['Condition'] == cond]['_y']
            sub_animals = animal_means[animal_means['Condition'] == cond]['_y']
            lines.append(
                f'{cond:<32} {len(sub_animals):>7d} {len(sub_cells):>8d} '
                f'{sub_cells.mean():>12.3f} {sub_cells.std():>10.3f}')

        # ── Pairwise tests ────────────────────────────────────────────────────
        pair_stats = []
        for g1, g2 in combinations(present, 2):
            va = df[df['Condition'] == g1]['_y'].dropna()
            vb = df[df['Condition'] == g2]['_y'].dropna()
            ma = animal_means[animal_means['Condition'] == g1]['_y'].dropna()
            mb = animal_means[animal_means['Condition'] == g2]['_y'].dropna()

            delta = vb.mean() - va.mean()

            # Cell-level Welch t-test (always computable when n>=2 in both).
            try:
                _, p_cell = _scstats.ttest_ind(va, vb, equal_var=False)
            except Exception:
                p_cell = float('nan')

            # Animal-level Welch t-test (the inferential one of record).
            if len(ma) >= 2 and len(mb) >= 2:
                try:
                    _, p_animal = _scstats.ttest_ind(ma, mb, equal_var=False)
                except Exception:
                    p_animal = float('nan')
            else:
                p_animal = float('nan')

            # Mann-Whitney U (cell-level, two-sided).
            try:
                _, p_mw = _scstats.mannwhitneyu(
                    va, vb, alternative='two-sided')
            except Exception:
                p_mw = float('nan')

            # Standard error of the cell-level mean difference (Welch).
            n_a, n_b = len(va), len(vb)
            if n_a > 1 and n_b > 1:
                se_delta = ((va.var(ddof=1) / n_a +
                             vb.var(ddof=1) / n_b)) ** 0.5
            else:
                se_delta = float('nan')

            pair_stats.append(
                (g1, g2, delta, se_delta, p_cell, p_animal, p_mw,
                 len(ma), len(mb), n_a, n_b))

        # Bonferroni correction across the family of pairwise tests.
        n_comp = len(pair_stats)

        def _bonf(p):
            if pd.isna(p):
                return p
            return min(1.0, p * n_comp)

        # ── The model of record: the one used in Prina et al., Brain 2026 ──────
        # Stata  xtmixed y i.genotype || animalid:, robust  +  pwcompare, mcompare(scheffe)
        # = linear mixed model with an animal random intercept (ML), cluster-robust
        # standard errors, Scheffé-adjusted contrasts.  Reproduces the published
        # tables exactly (mixed_stats.py, validated on the paper's datasheets).
        mm, mm_err = None, ''
        if mixed_stats is None:
            mm_err = 'mixed_stats.py not found next to this script'
        else:
            try:
                mm = mixed_stats.genotype_model(df['_y'].values, df['Animal'].values,
                                                df['Condition'].values, reference=ref_cond)
                if not mm['converged']:
                    mm_err = 'did not converge'
            except Exception as e:
                mm, mm_err = None, str(e)

        pairwise_rows = []
        self._pair_results = {}
        for (g1, g2, delta, se, p_cell, p_animal, p_mw,
             na_anim, nb_anim, na_cell, nb_cell) in pair_stats:
            p_cell_adj   = _bonf(p_cell)
            p_animal_adj = _bonf(p_animal)
            p_mw_adj     = _bonf(p_mw)
            cm = None
            if mm is not None:
                try:
                    cm = mixed_stats.contrast(mm, g1, g2)      # g2 − g1, like delta
                except Exception:
                    cm = None
            mm_p = cm['p_scheffe'] if cm and np.isfinite(cm['p_scheffe']) else float('nan')
            # Brackets show the mixed-model Scheffé p when the model fitted,
            # otherwise the animal-level Bonferroni p (cell-level when a group
            # has fewer than 2 animals).
            if not pd.isna(mm_p):
                bracket_p = mm_p
            else:
                bracket_p = p_animal_adj if not pd.isna(p_animal_adj) \
                                         else p_cell_adj
            self._pair_results[(g1, g2)] = bracket_p
            pairwise_rows.append({
                'g1': g1, 'g2': g2,
                'delta': delta, 'se': se,
                'mm_est': cm['estimate'] if cm else float('nan'),
                'mm_se':  cm['se_robust'] if cm else float('nan'),
                'mm_z':   cm['z'] if cm else float('nan'),
                'mm_p_raw': cm['p'] if cm else float('nan'),
                'mm_p_scheffe': mm_p,
                'p_cell_raw':   p_cell,   'p_cell_adj':   p_cell_adj,
                'p_animal_raw': p_animal, 'p_animal_adj': p_animal_adj,
                'p_mw_raw':     p_mw,     'p_mw_adj':     p_mw_adj,
                'bracket_p':    bracket_p,
                'n_animals_a':  na_anim,  'n_animals_b':  nb_anim,
                'n_cells_a':    na_cell,  'n_cells_b':    nb_cell,
            })

        # ── Build the popup text ──────────────────────────────────────────────
        def _fmt_p(p):
            if pd.isna(p):
                return '   n/a'
            if p < 0.0001:
                return '<0.0001'
            return f'{p:7.4f}'

        lines += ['', 'MIXED-EFFECTS MODEL  (test of record: Moen et al. 2016, PLoS One; as in Prina et al., Brain 2026)']
        if mm is not None:
            lines += [
                '  y ~ condition + (1 | animal), maximum likelihood; standard errors',
                '  cluster-robust on animal; p-values Scheffé-adjusted over the '
                f'{mm["k"]} conditions',
                f'  N cells {mm["n_obs"]:,}   N animals {mm["n_animals"]}   '
                f'var(animal) {mm["sigma2_animal"]:.4f}   var(residual) {mm["sigma2_residual"]:.4f}   '
                f'log-likelihood {mm["loglik"]:.3f}',
                '',
                f'{"Comparison (g2 − g1)":<60} {"contrast":>9} {"robust SE":>10} '
                f'{"z":>7} {"p Scheffé":>10}',
                '-' * 100,
            ]
            for r in pairwise_rows:
                label = f'{r["g1"][:28]} vs {r["g2"][:28]}'
                sig = _sig_label(r['mm_p_scheffe']) if not pd.isna(r['mm_p_scheffe']) else ''
                lines.append(
                    f'{label:<60} {r["mm_est"]:>+9.4f} {r["mm_se"]:>10.4f} '
                    f'{r["mm_z"]:>7.2f} {_fmt_p(r["mm_p_scheffe"]):>10}  {sig}')
        else:
            lines.append(f'  not available: {mm_err}')

        lines += [
            '',
            'DESCRIPTIVE PAIRWISE TESTS  (Bonferroni-adjusted; for reference only)',
            f'{"Comparison":<38} {"Δ mean":>8}   '
            f'{"p_animal":>8} {"p_cell":>8} {"p_MW":>8}',
            '-' * 78,
        ]
        for r in pairwise_rows:
            label = f'{r["g1"][:17]} vs {r["g2"][:17]}'
            lines.append(
                f'{label:<38} {r["delta"]:>+8.3f}   '
                f'{_fmt_p(r["p_animal_adj"]):>8} '
                f'{_fmt_p(r["p_cell_adj"]):>8} '
                f'{_fmt_p(r["p_mw_adj"]):>8}')

        lines += [
            '',
            'Test legend:',
            '  mixed model  contrast = difference of condition means from the model;',
            '               p Scheffé = P[chi2(k−1) ≥ z²] as Stata pwcompare, mcompare(scheffe).',
            '  p_animal     Welch t-test on per-animal means.',
            '  p_cell       Welch t-test on cell-level values (overstates significance).',
            '  p_MW         Mann-Whitney U on cell-level values.',
            '',
            'Brackets on the plot show the mixed-model Scheffé p-value; if the model',
            'could not be fitted, the animal-level Bonferroni p-value (cell-level when a',
            'group has fewer than 2 animals).',
            '',
            '* p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001',
            sep,
        ]

        self._model_data = {
            'outcome':    out_col,
            'reference':  ref_cond,
            'n_cells':    len(df),
            'n_animals':  df['Animal'].nunique(),
            # variance components of the mixed model (nan if it was not fitted)
            'animal_var': mm['sigma2_animal'] if mm else float('nan'),
            'residual':   mm['sigma2_residual'] if mm else float('nan'),
            'mixed':      {'fitted': mm is not None, 'error': mm_err,
                           'k': mm['k'] if mm else None,
                           'loglik': mm['loglik'] if mm else float('nan')},
            'pairwise':   pairwise_rows,
            'text':       '\n'.join(lines),
        }

        # ── Side-panel summary line ───────────────────────────────────────────
        ref_pairs = [r for r in pairwise_rows
                     if r['g1'] == ref_cond or r['g2'] == ref_cond]
        summary_parts = []
        for r in ref_pairs:
            other = r['g2'] if r['g1'] == ref_cond else r['g1']
            sig   = _sig_label(r['bracket_p']) if not pd.isna(r['bracket_p']) \
                                              else 'n/a'
            summary_parts.append(f'{other}: {sig}')
        self.lbl_model_summary.config(
            text='vs ' + ref_cond + '  |  ' + '   '.join(summary_parts)
                 + f'\nN={len(df):,} cells, {df["Animal"].nunique()} animals  '
                   '(animal-level Welch)')

        if self._results_win and self._results_win.winfo_exists():
            self._update_results_window()
        self._draw_plot()

    # ── Font / figure export ──────────────────────────────────────────────────

    def _on_font_change(self, *_):
        self._draw_plot()

    def _export_figure(self):
        if self.df is not None:
            path = filedialog.asksaveasfilename(
                title='Save figure',
                defaultextension='.png',
                filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'),
                           ('SVG', '*.svg'), ('TIFF', '*.tif')])
            if path:
                self.fig.savefig(path, dpi=300, bbox_inches='tight')
                messagebox.showinfo('Saved', os.path.basename(path))

    # ── Results popup ─────────────────────────────────────────────────────────

    def _show_results_window(self):
        if self._model_data is None:
            messagebox.showinfo('No results', 'Run the model first.')
            return
        if self._results_win and self._results_win.winfo_exists():
            self._results_win.lift()
            return
        win = tk.Toplevel(self.root)
        win.title('Model Results')
        win.geometry('640x540')
        self._results_win = win
        txt = scrolledtext.ScrolledText(win, font=('Courier', 9), bg='#F5F5F5')
        txt.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        self._results_text_widget = txt
        btn_row = tk.Frame(win)
        btn_row.pack(pady=(0, 8))
        tk.Button(btn_row, text='Export stats to xlsx',
                  font=('Arial', 9, 'bold'),
                  bg='#1565C0', fg='white', width=20,
                  command=self._export_stats).pack(side=tk.LEFT, padx=6)
        tk.Button(btn_row, text='Close',
                  font=('Arial', 9), width=10,
                  command=win.destroy).pack(side=tk.LEFT, padx=6)
        self._update_results_window()

    def _update_results_window(self):
        if not (self._results_win and self._results_win.winfo_exists()):
            return
        if not self._model_data:
            return
        self._results_text_widget.delete('1.0', tk.END)
        self._results_text_widget.insert('1.0', self._model_data.get('text', ''))

    # ── Stats xlsx export ─────────────────────────────────────────────────────

    def _export_stats(self):
        if not self._model_data or 'pairwise' not in self._model_data:
            messagebox.showinfo('No results', 'Run the model first.')
            return
        path = filedialog.asksaveasfilename(
            title='Save statistics',
            defaultextension='.xlsx',
            filetypes=[('Excel', '*.xlsx')])
        if not path:
            return

        from openpyxl.styles import Font, PatternFill, Alignment
        wb       = openpyxl.Workbook()
        hdr_fill = PatternFill('solid', start_color='1565C0')
        hdr_font = Font(bold=True, color='FFFFFF', name='Arial', size=10)
        dat_font = Font(name='Arial', size=10)

        ws = wb.active
        ws.title = 'Pairwise'
        headers  = ['Outcome', 'Group 1', 'Group 2',
                    'Δ mean', 'SE (cell-level)',
                    'n_animals_1', 'n_animals_2',
                    'n_cells_1', 'n_cells_2',
                    'Mixed model contrast (2 − 1)', 'Mixed model robust SE', 'Mixed model z',
                    'p_mixed_raw', 'p_mixed_Scheffe',
                    'p_animal_raw', 'p_animal_adj',
                    'p_cell_raw',   'p_cell_adj',
                    'p_MW_raw',     'p_MW_adj',
                    'Significance (mixed-model Scheffé; else animal-level Bonferroni)']
        for ci, h in enumerate(headers, 1):
            c = ws.cell(row=1, column=ci, value=h)
            c.font = hdr_font; c.fill = hdr_fill
            c.alignment = Alignment(horizontal='center')
        md = self._model_data

        def _r(x, n=6):
            try:
                v = float(x)
                return round(v, n) if np.isfinite(v) else ''
            except Exception:
                return ''

        for ri, r in enumerate(md['pairwise'], 2):
            bracket_p = r['bracket_p']
            row = [md['outcome'], r['g1'], r['g2'],
                   _r(r['delta'], 4), _r(r['se'], 4),
                   r['n_animals_a'], r['n_animals_b'],
                   r['n_cells_a'],   r['n_cells_b'],
                   _r(r.get('mm_est'), 4), _r(r.get('mm_se'), 4), _r(r.get('mm_z'), 2),
                   _r(r.get('mm_p_raw')), _r(r.get('mm_p_scheffe')),
                   _r(r['p_animal_raw']), _r(r['p_animal_adj']),
                   _r(r['p_cell_raw']),   _r(r['p_cell_adj']),
                   _r(r['p_mw_raw']),     _r(r['p_mw_adj']),
                   _sig_label(bracket_p) if not pd.isna(bracket_p) else 'n/a']
            for ci, v in enumerate(row, 1):
                ws.cell(row=ri, column=ci, value=v).font = dat_font
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = 16

        ws2  = wb.create_sheet('Model Info')
        info = [
            ('Outcome variable',        md['outcome']),
            ('Reference condition',     md['reference']),
            ('N cells',                 md['n_cells']),
            ('N animals',               md['n_animals']),
            ('Primary test',            'Linear mixed-effects model with a random intercept for each animal '
                                        '(cells nested within animals), fitted by maximum likelihood; '
                                        'standard errors cluster-robust on animal; pairwise contrasts '
                                        'Scheffé-adjusted. This is the clustered-data approach of Moen, '
                                        'Fricano-Kugler, Luikart & O\'Malley (2016) PLoS One 11:e0146721, '
                                        'as applied in Prina, Goyette et al. (2026) Brain '
                                        '(Stata: xtmixed … || animalid:, vce(robust); pwcompare, mcompare(scheffe))'),
            ('Implementation',          'mixed_stats.py (NumPy/SciPy) of the Soma Detector Suite v3.0; reproduces '
                                        'the published Stata contrasts, robust SEs, z and Scheffé p-values of '
                                        'Prina et al. 2026 to within 4e-7'),
            ('Mixed model fitted',      'yes' if md.get('mixed', {}).get('fitted') else
                                        'no: ' + str(md.get('mixed', {}).get('error', ''))),
            ('Conditions in model (k)', md.get('mixed', {}).get('k', '')),
            ('var(animal)',             _r(md.get('animal_var'))),
            ('var(residual)',           _r(md.get('residual'))),
            ('Log-likelihood',          _r(md.get('mixed', {}).get('loglik'), 4)),
            ('Secondary tests',         'Descriptive only: Welch t-test on per-animal means; Welch t-test '
                                        '(cell-level); Mann-Whitney U (cell-level); Bonferroni-adjusted'),
            ('Correction',              'Scheffé (mixed model, the test of record); Bonferroni (secondary tests)'),
        ]
        for ri, (k, v) in enumerate(info, 1):
            ws2.cell(row=ri, column=1, value=k).font = Font(bold=True, name='Arial', size=10)
            ws2.cell(row=ri, column=2, value=v).font = Font(name='Arial', size=10)
        ws2.column_dimensions['A'].width = 28
        ws2.column_dimensions['B'].width = 22

        ws3      = wb.create_sheet('Condition Summary')
        df_f     = self._filtered_df()
        out_col  = md['outcome']
        if out_col in df_f.columns:
            sum_hdrs = ['Condition', 'N animals', 'N cells', 'Median', 'Mean', 'SD']
            for ci, h in enumerate(sum_hdrs, 1):
                c = ws3.cell(row=1, column=ci, value=h)
                c.font = hdr_font; c.fill = hdr_fill
                c.alignment = Alignment(horizontal='center')
            for ri, cond in enumerate(self._conditions_order, 2):
                sub = pd.to_numeric(
                    df_f[df_f['Condition'] == cond][out_col], errors='coerce').dropna()
                if not len(sub):
                    continue
                n_an = df_f[df_f['Condition'] == cond]['Animal'].nunique()
                for ci, v in enumerate([cond, n_an, len(sub),
                                         round(sub.median(), 3),
                                         round(sub.mean(), 3),
                                         round(sub.std(), 3)], 1):
                    ws3.cell(row=ri, column=ci, value=v).font = dat_font
            for col in ws3.columns:
                ws3.column_dimensions[col[0].column_letter].width = 16

        # ── Individual Cells sheet (Stata-ready long format) ──────────────────
        ws4 = wb.create_sheet('Individual Cells')
        included = self._included_conditions()
        stata_hdrs = [out_col, 'animalid', 'conditionid']
        for ci, h in enumerate(stata_hdrs, 1):
            c = ws4.cell(row=1, column=ci, value=h)
            c.font = hdr_font; c.fill = hdr_fill
            c.alignment = Alignment(horizontal='center')
        if out_col in df_f.columns:
            stata_df = df_f[df_f['Condition'].isin(included)][
                [out_col, 'Animal', 'Condition']].copy()
            stata_df[out_col] = pd.to_numeric(stata_df[out_col], errors='coerce')
            stata_df = stata_df.dropna(subset=[out_col])
            # Sort by condition order then animal for clean Stata import
            cond_order_map = {c: i for i, c in enumerate(self._conditions_order)}
            stata_df['_ord'] = stata_df['Condition'].map(cond_order_map).fillna(99)
            stata_df = stata_df.sort_values(['_ord', 'Animal']).drop(columns='_ord')
            for ri, (_, row) in enumerate(stata_df.iterrows(), 2):
                ws4.cell(row=ri, column=1, value=round(float(row[out_col]), 4)).font = dat_font
                ws4.cell(row=ri, column=2, value=str(row['Animal'])).font = dat_font
                ws4.cell(row=ri, column=3, value=str(row['Condition'])).font = dat_font
        ws4.column_dimensions['A'].width = 28
        ws4.column_dimensions['B'].width = 18
        ws4.column_dimensions['C'].width = 28

        wb.save(path)
        messagebox.showinfo('Saved', os.path.basename(path))


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    root = tk.Tk()
    GroupAnalysisApp(root)
    root.mainloop()
