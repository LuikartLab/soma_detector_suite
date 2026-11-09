"""
Cre Classifier GUI  —  three linked 2D scatter plots
Open soma_measurements.xlsx, adjust Cre / Counter cutoffs via sliders,
classify somas by fluorescence phenotype, export results.

Phenotype groups (user-named):
  Cre-only      Cre channel >= cutoff  AND  Counter channel < cutoff
  Counter-only  Cre channel < cutoff   AND  Counter channel >= cutoff
  Both          both channels >= cutoff
  Unclassified  both below cutoff  (omitted from export)

Channel assignment:
  Auto-detected — any Mean column containing 'Cre' → Cre channel.
  Both channels can be reassigned via dropdowns.
"""

import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
import time
import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.cell.rich_text import CellRichText
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure   # not pyplot: a pyplot figure owns a hidden
                                       # Tk root that keeps the process alive after
                                       # the window is closed
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

COL_CRE  = '#2ca02c'   # green  – Cre channel  (matches GFP-Cre)
COL_CTR  = '#d62728'   # red    – Counter channel (matches mCherry)
COL_BOTH = '#9467bd'   # purple – Both positive

# Intensity metric the cutoffs apply to (per channel).  'Mean' is the raw soma
# mean (original behaviour).  'F/F0' divides it by that slice's background and
# 'Mean-BG' subtracts it: both remove the image-to-image background differences
# (green BG 13–32 across the Fig 1C iv stacks) that make one raw cutoff sit at
# a different place in each image's negative population.
METRICS = [('Mean', 'Mean (raw)'), ('F/F0', 'F/F0 (mean / background)'), ('Mean-BG', 'Mean - BG')]
METRIC_LABEL = dict(METRICS)
METRIC_DEFAULT_CUT = {'Mean': 25.0, 'F/F0': 1.8, 'Mean-BG': 10.0}

import json
import subprocess
try:
    import cell_jump_ipc                     # same folder; jump to a cell in the Soma Detector
except ImportError:                          # standalone copy of this script
    cell_jump_ipc = None


def _launch_soma_detector():
    """Start the Soma Detector next to this script (portable python312, the
    Cellpose venv, or the running interpreter).  True if a process was started."""
    here = os.path.dirname(os.path.abspath(__file__))
    script = next((os.path.join(here, n) for n in ('soma_area_and_intensity_gui_3.0.py',
                                                    'soma_area_and_intensity_gui_2.1.py')
                   if os.path.isfile(os.path.join(here, n))), None)
    if script is None:
        return False
    portable = os.path.join(here, 'python312', 'pythonw.exe')
    venv     = os.path.join(here, '.venv-cellpose', 'Scripts', 'pythonw.exe')
    pyw      = os.path.join(os.path.dirname(sys.executable), 'pythonw.exe')
    launcher = next((p for p in (portable, venv, pyw) if os.path.isfile(p)), sys.executable)
    try:
        subprocess.Popen([launcher, script], cwd=here)
        return True
    except OSError:
        return False
COL_UNDET = '#bbbbbb'  # grey   – Unclassified


def _detect_group(filename):
    """Derive experimental group from filename (Pten / Pten/Raptor / Pten/Rictor)."""
    fn = str(filename).upper()
    if 'RAPTOR' in fn:
        return 'Pten/Raptor'
    if 'RICTOR' in fn:
        return 'Pten/Rictor'
    return 'Pten'

# Internal phenotype keys
_CRE_ONLY = 'cre_only'
_CTR_ONLY = 'ctr_only'
_BOTH     = 'both'
_UNDET    = 'Unclassified'


class CreClassifierApp:

    def __init__(self, root):
        self.root = root
        self.root.title('Cre Classifier')
        self.root.geometry('1380x960')
        self.root.minsize(900, 800)

        self.df          = None
        self.raw_df      = None
        self.out_dir     = None
        self.soma_col    = None
        self.area_col    = None
        self.src_col     = None
        self.file_col    = None
        self.cre_col     = None
        self.counter_col = None
        self._mean_cols  = []
        # intensity metric per channel (see METRICS) and the last cutoff used
        # with each metric, so switching back restores the user's value
        self._metric_cre = tk.StringVar(value='Mean')
        self._metric_ctr = tk.StringVar(value='Mean')
        self._cut_memory = {}
        # bleed-through compensation: Cre positive when
        #   cre_metric >= cutoff + slope * counter_metric   (slope 0 = flat cutoff)
        self._slope_var = tk.DoubleVar(value=0.0)
        # per-cell manual overrides {cell key: 'Cre+' | 'Cre-' | 'Both positive'},
        # persisted beside the workbook; the selected cell (row index) for the
        # override buttons and the jump to the Soma Detector
        self._overrides  = {}
        self._selected   = None
        self._sel_rings  = []

        # one scatter per axis (updated in-place)
        self._scatters   = []
        self._lines_ch1  = []
        self._lines_ch2  = []

        self._build_ui()

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        top = tk.Frame(self.root, pady=4)
        top.pack(fill=tk.X, padx=8)

        tk.Button(top, text='Open soma_measurements.xlsx',
                  command=self._open_file,
                  bg='#1565C0', fg='white',
                  font=('Arial', 10, 'bold'), width=28).pack(side=tk.LEFT)

        self.lbl_file = tk.Label(top, text='No file loaded', fg='#555',
                                 font=('Arial', 9))
        self.lbl_file.pack(side=tk.LEFT, padx=8)

        self.lbl_counts = tk.Label(top, text='', fg='#222',
                                   font=('Arial', 10, 'bold'))
        self.lbl_counts.pack(side=tk.RIGHT, padx=12)

        main = tk.Frame(self.root)
        main.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        # ── figure panel ──────────────────────────────────────────────────────
        fig_frame = tk.Frame(main)
        fig_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.fig  = Figure(figsize=(13, 5.2))
        self.axes = self.fig.subplots(1, 3)
        self.fig.subplots_adjust(left=0.07, right=0.97,
                                  bottom=0.12, top=0.91, wspace=0.32)
        for ax in self.axes:
            ax.tick_params(labelsize=7)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        self.canvas = FigureCanvasTkAgg(self.fig, fig_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, fig_frame)
        toolbar.update()
        # click a point to select that cell (double-click: show it in the Soma Detector)
        self.canvas.mpl_connect('button_press_event', self._on_click)

        # ── control panel ─────────────────────────────────────────────────────
        ctrl = tk.Frame(main, width=290, padx=12, pady=10,
                        relief=tk.RIDGE, bd=1)
        ctrl.pack(side=tk.RIGHT, fill=tk.Y, padx=(8, 0))
        ctrl.pack_propagate(False)

        # Channel assignment
        ch_frame = tk.LabelFrame(ctrl, text='Channel Assignment',
                                  font=('Arial', 9), padx=6, pady=6)
        ch_frame.pack(fill=tk.X, pady=(0, 6))

        tk.Label(ch_frame, text='Cre channel:',
                 font=('Arial', 8, 'bold'), fg=COL_CRE).grid(
                     row=0, column=0, sticky='w', pady=2)
        self._cre_menu_var = tk.StringVar(value='(load file)')
        self._cre_menu = tk.OptionMenu(ch_frame, self._cre_menu_var,
                                        '(load file)')
        self._cre_menu.config(font=('Arial', 8), width=13, anchor='w')
        self._cre_menu.grid(row=0, column=1, sticky='ew', padx=(4, 0))

        tk.Label(ch_frame, text='Counter channel:',
                 font=('Arial', 8, 'bold'), fg=COL_CTR).grid(
                     row=1, column=0, sticky='w', pady=2)
        self._counter_menu_var = tk.StringVar(value='(load file)')
        self._counter_menu = tk.OptionMenu(ch_frame, self._counter_menu_var,
                                            '(load file)')
        self._counter_menu.config(font=('Arial', 8), width=13, anchor='w')
        self._counter_menu.grid(row=1, column=1, sticky='ew', padx=(4, 0))

        # Intensity metric the cutoffs apply to, per channel
        for row_i, (lbl, var, col) in enumerate(
                [('Cre metric:', self._metric_cre, COL_CRE),
                 ('Counter metric:', self._metric_ctr, COL_CTR)], start=2):
            tk.Label(ch_frame, text=lbl, font=('Arial', 8, 'bold'),
                     fg=col).grid(row=row_i, column=0, sticky='w', pady=2)
            disp = tk.StringVar(value=METRIC_LABEL[var.get()])
            om = tk.OptionMenu(ch_frame, disp, *[l for _, l in METRICS],
                               command=lambda label, v=var: self._set_metric(v, label))
            om.config(font=('Arial', 8), width=13, anchor='w')
            om.grid(row=row_i, column=1, sticky='ew', padx=(4, 0))

        ch_frame.columnconfigure(1, weight=1)

        # Cutoff sliders
        tk.Label(ctrl, text='Classification Cutoffs',
                 font=('Arial', 12, 'bold')).pack(anchor='w', pady=(0, 6))

        self.ch1_var = tk.DoubleVar(value=25.0)
        self.ch2_var = tk.DoubleVar(value=25.0)
        self._ch1_label_var = tk.StringVar(value='Cre channel minimum:')
        self._ch2_label_var = tk.StringVar(value='Counter channel minimum:')

        self._make_cutoff_widget(ctrl,
            label_var=self._ch1_label_var,
            var=self.ch1_var, slider_attr='sl_ch1', color=COL_CRE)

        # Bleed-through compensation: tilt the Cre cutoff with the counter signal
        srow = tk.Frame(ctrl)
        srow.pack(anchor='w', pady=(2, 0))
        tk.Label(srow, text='Bleed-through slope:', font=('Arial', 8),
                 fg=COL_CRE).pack(side=tk.LEFT)
        sp = tk.Spinbox(srow, from_=-1.0, to=1.0, increment=0.005,
                        textvariable=self._slope_var, width=7,
                        font=('Arial', 9), format='%.3f', command=self._on_slope_change)
        sp.pack(side=tk.LEFT, padx=(4, 2))
        sp.bind('<Return>',   lambda e: self._on_slope_change())
        sp.bind('<FocusOut>', lambda e: self._on_slope_change())
        tk.Button(srow, text='Fit', font=('Arial', 8), width=4,
                  command=self._fit_slope).pack(side=tk.LEFT, padx=2)
        tk.Button(srow, text='0', font=('Arial', 8), width=2,
                  command=lambda: (self._slope_var.set(0.0), self._on_slope_change())).pack(side=tk.LEFT)

        tk.Frame(ctrl, height=6).pack()

        self._make_cutoff_widget(ctrl,
            label_var=self._ch2_label_var,
            var=self.ch2_var, slider_attr='sl_ch2', color=COL_CTR)

        tk.Frame(ctrl, height=6).pack()

        # Genotype group names
        grp_frame = tk.LabelFrame(ctrl, text='Genotype Group Names  (used in Group column)',
                                   font=('Arial', 9), padx=6, pady=6)
        grp_frame.pack(fill=tk.X, pady=4)

        self._name_cre  = tk.StringVar(value='Cre-only')
        self._name_ctr  = tk.StringVar(value='Counter-only')
        self._name_both = tk.StringVar(value='Both positive')

        self._grp_rows = [
            (COL_CRE,  '● Cre-only:',      self._name_cre),
            (COL_CTR,  '● Counter-only:',   self._name_ctr),
            (COL_BOTH, '● Both positive:',  self._name_both),
        ]
        for row_i, (col, lbl_text, var) in enumerate(self._grp_rows):
            tk.Label(grp_frame, text=lbl_text, fg=col,
                     font=('Arial', 8, 'bold')).grid(
                         row=row_i, column=0, sticky='w', pady=2)
            e = tk.Entry(grp_frame, textvariable=var, width=16,
                         font=('Arial', 8))
            e.grid(row=row_i, column=1, sticky='ew', padx=(4, 0))

        grp_frame.columnconfigure(1, weight=1)

        # Colour legend
        leg = tk.LabelFrame(ctrl, text='Classification', font=('Arial', 9),
                             padx=6, pady=4)
        leg.pack(fill=tk.X, pady=4)
        for col, text in [
                (COL_CRE,  '● Cre-only    Cre ≥ cut, Counter < cut'),
                (COL_CTR,  '● Counter-only   Cre < cut, Counter ≥ cut'),
                (COL_BOTH, '● Both positive   both ≥ cutoff'),
                (COL_UNDET,'● Unclassified   both below')]:
            tk.Label(leg, text=text, fg=col,
                     font=('Arial', 8), anchor='w').pack(anchor='w')

        # Selected cell: manual override + jump to the Soma Detector
        sel = tk.LabelFrame(ctrl, text='Selected cell  (click a point; double-click opens it)',
                            font=('Arial', 9), padx=6, pady=4)
        sel.pack(fill=tk.X, pady=4)
        self.lbl_sel = tk.Label(sel, text='none', font=('Arial', 8), fg='#333',
                                justify='left', anchor='w', wraplength=250)
        self.lbl_sel.pack(anchor='w')
        brow = tk.Frame(sel)
        brow.pack(anchor='w', pady=(3, 0))
        for text, val, col, key in (('Cre+', 'Cre+', COL_CRE, '1'),
                                    ('Cre-', 'Cre-', COL_CTR, '2'),
                                    ('Both', 'Both positive', COL_BOTH, '3'),
                                    ('Auto', None, '#555', '0')):
            tk.Button(brow, text='{} ({})'.format(text, key), fg=col, width=7,
                      font=('Arial', 8, 'bold'),
                      command=lambda v=val: self._set_override(v)).pack(side=tk.LEFT, padx=1)
            self.root.bind(key, lambda e, v=val: self._set_override(v))
        tk.Button(sel, text='Show in Soma Detector', font=('Arial', 8),
                  command=self._show_in_detector).pack(anchor='w', pady=(4, 0))

        tk.Frame(ctrl, height=4).pack()

        tk.Button(ctrl,
                  text='Export classified workbook',
                  command=self._export,
                  bg='#2E7D32', fg='white',
                  font=('Arial', 10, 'bold'),
                  width=26, height=2).pack(pady=6)

        self.lbl_status = tk.Label(ctrl, text='', fg='#333',
                                   font=('Arial', 8),
                                   wraplength=250, justify='left')
        self.lbl_status.pack(anchor='w')

    def _make_cutoff_widget(self, parent, label_var, var, slider_attr, color):
        tk.Label(parent, textvariable=label_var,
                 font=('Arial', 9, 'bold'), fg=color).pack(anchor='w')

        row = tk.Frame(parent)
        row.pack(anchor='w', pady=(2, 0))

        spin = tk.Spinbox(row, from_=0, to=99999, increment=0.2,
                          textvariable=var, width=8,
                          font=('Arial', 10), format='%.1f',
                          command=self._update)
        spin.pack(side=tk.LEFT)
        spin.bind('<Return>',   lambda e: self._update())
        spin.bind('<FocusOut>', lambda e: self._update())
        setattr(self, slider_attr.replace('sl_', 'sp_'), spin)

        sl = tk.Scale(parent, from_=0, to=500, resolution=0.2,
                      orient=tk.HORIZONTAL, variable=var,
                      length=250, troughcolor=color,
                      command=lambda _: self._update())
        sl.pack(anchor='w')
        setattr(self, slider_attr, sl)

    # ── Channel assignment ─────────────────────────────────────────────────────

    def _short_name(self, col):
        return col[5:] if col.startswith('Mean ') else col

    def _populate_channel_menus(self):
        short_names = [self._short_name(c) for c in self._mean_cols]
        for menu_widget, var in [(self._cre_menu,     self._cre_menu_var),
                                  (self._counter_menu, self._counter_menu_var)]:
            menu = menu_widget['menu']
            menu.delete(0, 'end')
            for name in short_names:
                menu.add_command(
                    label=name,
                    command=lambda n=name, v=var: self._set_channel(v, n))
        self._cre_menu_var.set(self._short_name(self.cre_col))
        self._counter_menu_var.set(self._short_name(self.counter_col))

    def _set_channel(self, var, name):
        var.set(name)
        self._on_channel_change()

    def _on_channel_change(self):
        if self.df is None:
            return
        cre_short     = self._cre_menu_var.get()
        counter_short = self._counter_menu_var.get()
        cre_col     = next((c for c in self._mean_cols
                            if self._short_name(c) == cre_short), None)
        counter_col = next((c for c in self._mean_cols
                            if self._short_name(c) == counter_short), None)
        if not cre_col or not counter_col or cre_col == counter_col:
            return
        self.cre_col     = cre_col
        self.counter_col = counter_col
        self._update_dynamic_labels()
        self._update_group_name_defaults()
        self._prepare_data()
        self._configure_cutoff_widgets()
        self._init_plots()
        self._update()

    def _update_dynamic_labels(self):
        cre_name = self._short_name(self.cre_col)
        ctr_name = self._short_name(self.counter_col)
        self._ch1_label_var.set(f'Cre channel ({cre_name}) {self._metric_cre.get()} minimum:')
        self._ch2_label_var.set(f'Counter channel ({ctr_name}) {self._metric_ctr.get()} minimum:')

    def _update_group_name_defaults(self):
        """Set sensible default group names from the actual channel names."""
        cre_name = self._short_name(self.cre_col)
        ctr_name = self._short_name(self.counter_col)
        self._name_cre.set(f'{cre_name} only')
        self._name_ctr.set(f'{ctr_name} only')
        self._name_both.set(f'{ctr_name} + {cre_name}')

    # ── Data loading ──────────────────────────────────────────────────────────

    def _open_file(self):
        path = filedialog.askopenfilename(
            title='Open soma_measurements.xlsx',
            filetypes=[('Excel files', '*.xlsx'), ('All files', '*.*')])
        if not path:
            return
        try:
            self.raw_df   = pd.read_excel(path, sheet_name='Cells', header=1)
            self._load_genotype_map(path)
            self.out_dir  = os.path.dirname(path)
            self.in_fname = os.path.basename(path)
            self.lbl_file.config(text=self.in_fname)
            self._resolve_columns()
            self._prepare_data()
            self._load_overrides()
            self._configure_cutoff_widgets()
            self._init_plots()
            self._update()
            self.lbl_status.config(
                text='Loaded {:,} unique somas.\nCre: {}\nCounter: {}'.format(
                    len(self.df),
                    self._short_name(self.cre_col),
                    self._short_name(self.counter_col)))
        except Exception as e:
            messagebox.showerror('Error loading file', str(e))

    def _load_genotype_map(self, path):
        """Read the Genotype column per File directly from the workbook, so the
        rich-text (italic/superscript) formatting survives — pandas flattens it.
        Genotype is per-file, so one value per File is sufficient."""
        self._geno_rt_by_file    = {}
        self._geno_plain_by_file = {}
        try:
            wb = openpyxl.load_workbook(path, rich_text=True)
            ws = wb['Cells']
            hdr = {}
            for c in range(1, ws.max_column + 1):
                v = ws.cell(row=2, column=c).value
                if v is not None:
                    hdr[str(v).strip()] = c
            fcol, gcol = hdr.get('File'), hdr.get('Genotype')
            if not fcol or not gcol:
                return
            for r in range(3, ws.max_row + 1):
                fv = ws.cell(row=r, column=fcol).value
                if fv is None:
                    continue
                fv = str(fv)
                if fv in self._geno_rt_by_file:
                    continue
                gv = ws.cell(row=r, column=gcol).value
                self._geno_rt_by_file[fv]    = gv
                self._geno_plain_by_file[fv] = '' if gv is None else str(gv)
        except Exception:
            pass

    def _resolve_columns(self):
        cols = [str(c).strip() for c in self.raw_df.columns]
        self.raw_df.columns = cols

        self.soma_col = next(
            (c for c in cols if 'Soma' in c and '#' in c), cols[3])
        self.src_col  = next(
            (c for c in cols if 'Source' in c), None)
        self.file_col = next(
            (c for c in cols if c.lower() == 'file'), cols[0])
        self.area_col = next(
            (c for c in cols if 'Area' in c), None)

        self._mean_cols = [c for c in cols if c.startswith('Mean')]

        cre_cols = [c for c in self._mean_cols if 'cre' in c.lower()]
        non_cre  = [c for c in self._mean_cols if 'cre' not in c.lower()]

        if cre_cols:
            self.cre_col     = cre_cols[0]
            self.counter_col = (non_cre[0] if non_cre
                                else (self._mean_cols[1]
                                      if len(self._mean_cols) > 1
                                      else self._mean_cols[0]))
        else:
            self.cre_col     = self._mean_cols[0] if self._mean_cols else None
            self.counter_col = (self._mean_cols[1]
                                if len(self._mean_cols) > 1 else None)

        self.det_col = next(
            (c for c in cols if 'Channel' in c and 'detected' in c.lower()), None)
        self._cre_name_by_file = self._cre_channel_names_by_file()

        self._populate_channel_menus()
        self._update_dynamic_labels()
        self._update_group_name_defaults()

    # ── Intensity metric ──────────────────────────────────────────────────────

    def _metric_values(self, mean_col, metric, df=None):
        """Per-cell values of `metric` for the channel whose raw column is
        `mean_col`: the raw mean, F/F0 (mean ÷ slice background) or mean − BG.
        F/F0 and BG columns are exported by the detector next to each Mean
        column; when absent the raw mean is used."""
        df = self.df if df is None else df
        m = df[mean_col].values.astype(float)
        if metric == 'Mean' or mean_col is None:
            return m
        name = self._short_name(mean_col)
        if metric == 'F/F0':
            if 'F/F0 ' + name in df.columns:
                return df['F/F0 ' + name].values.astype(float)
            if 'BG ' + name in df.columns:
                bg = df['BG ' + name].values.astype(float)
                return np.where(bg > 0, m / np.where(bg > 0, bg, 1.0), np.nan)
        if metric == 'Mean-BG' and 'BG ' + name in df.columns:
            return m - df['BG ' + name].values.astype(float)
        return m

    def _set_metric(self, var, label):
        key = next((k for k, l in METRICS if l == label), 'Mean')
        if key == var.get():
            return
        # remember the cutoff used with the metric we are leaving
        cut_var = self.ch1_var if var is self._metric_cre else self.ch2_var
        self._cut_memory[(id(var), var.get())] = cut_var.get()
        var.set(key)
        # widen/narrow the slider BEFORE restoring the value: a Scale clamps its
        # variable to the range it currently has
        self._configure_cutoff_widgets()
        cut_var.set(self._cut_memory.get((id(var), key), METRIC_DEFAULT_CUT[key]))
        if self.df is None:
            return
        self._update_dynamic_labels()
        self._init_plots()
        self._update()

    def _configure_cutoff_widgets(self):
        """Slider range / resolution and spinbox step for each channel's metric
        (raw means run to ~255, F/F0 to ~15)."""
        if self.df is None:
            return
        for col, mvar, sl, sp in [(self.cre_col, self._metric_cre, self.sl_ch1, self.sp_ch1),
                                  (self.counter_col, self._metric_ctr, self.sl_ch2, self.sp_ch2)]:
            if col is None:
                continue
            vals = self._metric_values(col, mvar.get())
            vals = vals[np.isfinite(vals)]
            vmax = float(vals.max()) if len(vals) else 1.0
            vmin = float(vals.min()) if len(vals) else 0.0
            if mvar.get() == 'F/F0':
                sl.config(from_=0, to=float(np.ceil(vmax)) + 1, resolution=0.05)
                sp.config(from_=0, increment=0.05, format='%.2f')
            else:
                sl.config(from_=min(0, int(np.floor(vmin))), to=int(vmax) + 20, resolution=0.2)
                sp.config(from_=min(0, int(np.floor(vmin))), increment=0.2, format='%.1f')

    # ── Selected cell, manual overrides, jump to the detector ─────────────────

    def _cell_keys(self):
        """Stable per-cell keys (file | slice | centroid) — they survive a
        re-export of the workbook, which renumbers the somas."""
        df = self.df
        cx = next((c for c in df.columns if 'Centroid X' in c), None)
        cy = next((c for c in df.columns if 'Centroid Y' in c), None)
        sl = next((c for c in df.columns if c.strip().lower() == 'slice'), None)
        f  = df[self.file_col].astype(str).values
        z  = df[sl].values if sl else np.zeros(len(df))
        x  = df[cx].values.astype(float) if cx else np.zeros(len(df))
        y  = df[cy].values.astype(float) if cy else np.zeros(len(df))
        return ['{}|{}|{:.1f}|{:.1f}'.format(fi, zi, xi, yi) for fi, zi, xi, yi in zip(f, z, x, y)]

    def _overrides_path(self):
        return os.path.join(self.out_dir, self.in_fname + '.overrides.json')

    def _load_overrides(self):
        self._overrides = {}
        self._selected = None
        try:
            with open(self._overrides_path(), encoding='utf-8') as fh:
                data = json.load(fh)
            keys = set(self._cell_keys())
            self._overrides = {k: v for k, v in data.items() if k in keys
                               and v in ('Cre+', 'Cre-', 'Both positive')}
            if self._slope_var.get() == 0.0 and isinstance(data.get('__slope__'), (int, float)):
                self._slope_var.set(float(data['__slope__']))
        except (OSError, ValueError):
            pass

    def _save_overrides(self):
        try:
            data = dict(self._overrides)
            data['__slope__'] = float(self._slope_var.get())
            with open(self._overrides_path(), 'w', encoding='utf-8') as fh:
                json.dump(data, fh, indent=1)
        except OSError as e:
            messagebox.showwarning('Overrides not saved', str(e))

    def _on_click(self, event):
        if self.df is None or event.inaxes is None or event.button != 1:
            return
        try:
            ai = list(self.axes).index(event.inaxes)
        except ValueError:
            return
        # nearest point in screen pixels (axes have very different scales)
        xy = self._scatters[ai].get_offsets()
        pts = event.inaxes.transData.transform(np.asarray(xy))
        d = np.hypot(pts[:, 0] - event.x, pts[:, 1] - event.y)
        i = int(np.nanargmin(d)) if len(d) else None
        if i is None or d[i] > 12:
            return
        self._select(i)
        if event.dblclick:
            self._show_in_detector()

    def _select(self, i):
        self._selected = i
        row = self.df.iloc[i]
        lab = self._classify()[i]
        names = {_CRE_ONLY: 'Cre+', _CTR_ONLY: 'Cre-', _BOTH: 'Both positive', _UNDET: 'Unclassified'}
        key = self._cell_keys()[i]
        ov = self._overrides.get(key)
        cre_v = self._metric_values(self.cre_col, self._metric_cre.get())[i]
        ctr_v = self._metric_values(self.counter_col, self._metric_ctr.get())[i]
        sl = next((c for c in self.df.columns if c.strip().lower() == 'slice'), None)
        self.lbl_sel.config(text='{}\nSoma #{}  slice {}  area {:.0f} µm²\n{} {} = {:.2f}   {} {} = {:.2f}\nclass: {}{}'.format(
            str(row[self.file_col])[:36], row[self.soma_col], row[sl] if sl else '?',
            float(row[self.area_col]) if self.area_col else float('nan'),
            self._metric_cre.get(), self._short_name(self.cre_col), cre_v,
            self._metric_ctr.get(), self._short_name(self.counter_col), ctr_v,
            names.get(lab, '?'), '  (manual override)' if ov else ''))
        for ring, sc in zip(self._sel_rings, self._scatters):
            ring.set_offsets(np.asarray(sc.get_offsets())[[i]])
            ring.set_visible(True)
        self.canvas.draw_idle()

    def _set_override(self, value):
        if self.df is None or self._selected is None:
            return
        key = self._cell_keys()[self._selected]
        if value is None:
            self._overrides.pop(key, None)
        else:
            self._overrides[key] = value
        self._save_overrides()
        self._update()
        self._select(self._selected)

    def _show_in_detector(self):
        if self.df is None or self._selected is None:
            return
        if cell_jump_ipc is None:
            messagebox.showinfo('Soma Detector link', 'cell_jump_ipc.py is not next to this script.')
            return
        row = self.df.iloc[self._selected]
        sl = next((c for c in self.df.columns if c.strip().lower() == 'slice'), None)
        cx = next((c for c in self.df.columns if 'Centroid X' in c), None)
        cy = next((c for c in self.df.columns if 'Centroid Y' in c), None)
        names = {_CRE_ONLY: 'Cre+', _CTR_ONLY: 'Cre-', _BOTH: 'Both positive'}
        payload = {'action': 'open_cell', 'file': str(row[self.file_col]),
                   'search_dirs': [self.out_dir],
                   'slice': int(row[sl]) if sl and pd.notna(row[sl]) else None,
                   'soma_num': int(row[self.soma_col]) if pd.notna(row[self.soma_col]) else None,
                   'centroid_x_um': float(row[cx]) if cx else None,
                   'centroid_y_um': float(row[cy]) if cy else None,
                   'condition': names.get(self._classify()[self._selected], '')}
        if cell_jump_ipc.send_jump(payload):
            return
        if not _launch_soma_detector():
            messagebox.showinfo('Soma Detector not found',
                                'Could not reach or launch the Soma Detector.')
            return

        def _retry(n=0):
            if cell_jump_ipc.send_jump(payload):
                return
            if n < 40:                       # ~12 s for a fresh launch
                self.root.after(300, lambda: _retry(n + 1))
            else:
                messagebox.showinfo('Soma Detector starting',
                                    'The Soma Detector is opening but did not respond in time; '
                                    'click "Show in Soma Detector" again once it has loaded.')
        self.root.after(300, _retry)

    def _on_slope_change(self, *_):
        if self.df is not None:
            self._save_overrides()          # the slope is stored with the overrides
        self._update()

    def _fit_slope(self):
        """Bleed-through slope = robust (Theil–Sen) slope of the Cre metric
        against the counter metric among the cells that are Cre-negative at the
        current flat cutoff."""
        if self.df is None:
            return
        from scipy import stats
        cre = self._metric_values(self.cre_col, self._metric_cre.get())
        ctr = self._metric_values(self.counter_col, self._metric_ctr.get())
        neg = np.isfinite(cre) & np.isfinite(ctr) & (cre < self.ch1_var.get())
        if neg.sum() < 10:
            messagebox.showinfo('Fit', 'Fewer than 10 Cre-negative cells at the current cutoff.')
            return
        slope = float(stats.theilslopes(cre[neg], ctr[neg])[0])
        self._slope_var.set(round(slope, 3))
        self._save_overrides()
        self._update()

    def _cre_channel_names_by_file(self):
        """{file: Cre-channel name of that image}.

        Channel names are set per image in the detector (e.g. GFP-AKT1-Cre in
        the Akt1 animals, GFP-AKT3-Cre in the Akt3 animals) but the workbook's
        column headers only carry the first file's names.  The per-cell
        'Channel detected' / 'Source channel' values do carry each image's own
        names, so the Cre channel of a file is the detection-channel name that
        is not the counter channel.  Files where that cannot be recovered
        (no Cre-detected cell at all) fall back to the header name."""
        cre_default = self._short_name(self.cre_col) if self.cre_col else ''
        ctr_name    = (self._short_name(self.counter_col)
                       if self.counter_col else '').lower()
        cols = [c for c in (self.det_col, self.src_col) if c]
        out  = {}
        if not cols or self.file_col not in self.raw_df.columns:
            return out
        for fn, grp in self.raw_df.groupby(self.file_col):
            names = set()
            for c in cols:
                for v in grp[c].dropna().astype(str):
                    names.update(p.strip() for p in v.split('+'))
            names = {n for n in names if n and n.lower() != ctr_name}
            out[str(fn)] = names.pop() if len(names) == 1 else cre_default
        return out

    def _cre_name_for_file(self, fn):
        default = self._short_name(self.cre_col) if self.cre_col else ''
        return (getattr(self, '_cre_name_by_file', None) or {}).get(str(fn), default)

    def _prepare_data(self):
        """Deduplicate to one row per soma.

        For somas detected in only one channel: keep that row.
        For somas detected in both channels: randomly pick one row,
        using a fixed seed so results are reproducible across reloads.
        """
        df  = self.raw_df.copy()
        # A soma is identified by (File, Soma #): soma numbers repeat across
        # files, so deduplicating on the number alone would wrongly collapse
        # somas from different files into one.
        df['_uid'] = df[self.file_col].astype(str) + '||' + \
            df[self.soma_col].astype(str)
        src = self.src_col
        if src:
            cre_short  = self._short_name(self.cre_col)
            # A row comes from the Cre channel when its source channel is the
            # file's own Cre-channel name (per-image names differ, e.g.
            # GFP-AKT1-Cre vs GFP-AKT3-Cre); the header-name substring match
            # is kept for older workbooks without per-cell channel names.
            cre_mask   = (df[src].str.contains(cre_short, na=False, case=False)
                          | pd.Series([str(s).strip().lower()
                                       == self._cre_name_for_file(f).lower()
                                       for s, f in zip(df[src], df[self.file_col])],
                                      index=df.index))
            cre_uids   = set(df.loc[ cre_mask, '_uid'])
            ctr_uids   = set(df.loc[~cre_mask, '_uid'])
            dual_uids  = cre_uids & ctr_uids

            if dual_uids:
                # Assign each row an independent random priority score,
                # then keep the highest-scoring row per soma.  This gives
                # a true per-soma 50/50 coin flip regardless of dataset size.
                rng      = np.random.default_rng(42)
                dual_df  = df[df['_uid'].isin(dual_uids)].copy()
                dual_df['_p'] = rng.random(len(dual_df))
                dual   = (dual_df.sort_values('_p', ascending=False)
                                 .drop_duplicates(subset=['_uid'])
                                 .drop(columns=['_p']))
                single = df[~df['_uid'].isin(dual_uids)]
                dedup  = pd.concat([dual, single])
            else:
                # No dual-channel somas — keep Cre rows, fall back to counter rows
                other = df[~df['_uid'].isin(cre_uids)]
                dedup = pd.concat([df[cre_mask], other])
        else:
            dedup = df.copy()
        self.df = (dedup.drop_duplicates(subset=['_uid'])
                        .drop(columns=['_uid'])
                        .reset_index(drop=True))

    # ── Plot initialisation ────────────────────────────────────────────────────

    def _init_plots(self):
        self._scatters  = []
        self._lines_ch1 = []
        self._lines_ch2 = []

        cre_vals  = self._metric_values(self.cre_col, self._metric_cre.get())
        ctr_vals  = self._metric_values(self.counter_col, self._metric_ctr.get())
        area_vals = self.df[self.area_col].values.astype(float)

        cre_name = self._short_name(self.cre_col)
        ctr_name = self._short_name(self.counter_col)
        cre_m, ctr_m = self._metric_cre.get(), self._metric_ctr.get()

        for ax in self.axes:
            ax.cla()

        ax_data = [
            (cre_vals,  ctr_vals),
            (area_vals, cre_vals),
            (area_vals, ctr_vals),
        ]
        titles  = [f'{cre_name}  vs  {ctr_name}',
                   f'Area  vs  {cre_name}',
                   f'Area  vs  {ctr_name}']
        xlabels = [f'{cre_m} {cre_name}', 'Area (µm²)', 'Area (µm²)']
        ylabels = [f'{ctr_m} {ctr_name}', f'{cre_m} {cre_name}', f'{ctr_m} {ctr_name}']

        n = len(cre_vals)
        self._sel_rings = []
        self._selected  = None
        if hasattr(self, 'lbl_sel'):
            self.lbl_sel.config(text='none')
        for ai, ax in enumerate(self.axes):
            xd, yd = ax_data[ai]
            sc = ax.scatter(xd, yd,
                            c=[COL_UNDET] * n,
                            s=16, alpha=0.75, marker='o', linewidths=0)
            self._scatters.append(sc)
            ring = ax.scatter([], [], s=120, facecolors='none', edgecolors='black',
                              linewidths=1.6, zorder=5)
            ring.set_visible(False)
            self._sel_rings.append(ring)

        for ax, t, xl, yl in zip(self.axes, titles, xlabels, ylabels):
            ax.set_title(t, fontsize=9, fontweight='bold')
            ax.set_xlabel(xl, fontsize=8)
            ax.set_ylabel(yl, fontsize=8)
            ax.tick_params(labelsize=7)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        ch1 = self.ch1_var.get()
        ch2 = self.ch2_var.get()

        ax0 = self.axes[0]
        l1, = ax0.plot([ch1, ch1], ax0.get_ylim(),
                       color=COL_CRE, lw=1.2, ls='--', alpha=0.8)
        l2, = ax0.plot(ax0.get_xlim(), [ch2, ch2],
                       color=COL_CTR, lw=1.2, ls='--', alpha=0.8)
        self._lines_ch1.append(l1)
        self._lines_ch2.append(l2)

        ax1 = self.axes[1]
        l1, = ax1.plot(ax1.get_xlim(), [ch1, ch1],
                       color=COL_CRE, lw=1.2, ls='--', alpha=0.8)
        self._lines_ch1.append(l1)
        self._lines_ch2.append(None)

        ax2 = self.axes[2]
        l2, = ax2.plot(ax2.get_xlim(), [ch2, ch2],
                       color=COL_CTR, lw=1.2, ls='--', alpha=0.8)
        self._lines_ch1.append(None)
        self._lines_ch2.append(l2)

        self.canvas.draw_idle()

    # ── Classification ─────────────────────────────────────────────────────────

    def _classify(self):
        """Return array of _CRE_ONLY / _CTR_ONLY / _BOTH / _UNDET per cell."""
        cre = self._metric_values(self.cre_col, self._metric_cre.get())
        ctr = self._metric_values(self.counter_col, self._metric_ctr.get())
        ch1 = self.ch1_var.get()
        ch2 = self.ch2_var.get()
        k   = self._slope_var.get()
        # tilted Cre cutoff: bleed-through of the counter signal into the Cre
        # channel raises the negatives' Cre level in proportion to the counter
        cre_pos = cre >= ch1 + k * np.nan_to_num(ctr)
        ctr_pos = ctr >= ch2
        result = np.full(len(cre), _UNDET, dtype=object)
        result[ cre_pos & ~ctr_pos] = _CRE_ONLY
        result[~cre_pos &  ctr_pos] = _CTR_ONLY
        result[ cre_pos &  ctr_pos] = _BOTH
        if self._overrides:
            code = {'Cre+': _CRE_ONLY, 'Cre-': _CTR_ONLY, 'Both positive': _BOTH}
            for i, key in enumerate(self._cell_keys()):
                ov = self._overrides.get(key)
                if ov in code:
                    result[i] = code[ov]
        return result

    def _phenotype_color(self, labels):
        return np.where(labels == _CRE_ONLY, COL_CRE,
               np.where(labels == _CTR_ONLY, COL_CTR,
               np.where(labels == _BOTH,     COL_BOTH, COL_UNDET)))

    # ── Live update ────────────────────────────────────────────────────────────

    def _update(self, *_):
        if self.df is None or not self._scatters:
            return

        labels = self._classify()
        n_cre  = (labels == _CRE_ONLY).sum()
        n_ctr  = (labels == _CTR_ONLY).sum()
        n_both = (labels == _BOTH).sum()
        n_und  = (labels == _UNDET).sum()

        n_ov = len(self._overrides)
        self.lbl_counts.config(
            text='Cre+: {:,}   Cre-: {:,}   Both: {:,}   Unclassified: {:,}{}'.format(
                n_cre, n_ctr, n_both, n_und,
                '   |   manual overrides: {:,}'.format(n_ov) if n_ov else ''))

        colors = self._phenotype_color(labels)
        # overridden cells get a black rim so the manual decisions stay visible
        keys = self._cell_keys() if n_ov else None
        edges = (['black' if k in self._overrides else 'none' for k in keys]
                 if keys else 'none')
        widths = ([1.2 if k in self._overrides else 0 for k in keys] if keys else 0)
        for sc in self._scatters:
            sc.set_facecolors(colors)
            sc.set_edgecolors(edges)
            sc.set_linewidths(widths)

        ch1 = self.ch1_var.get()
        ch2 = self.ch2_var.get()
        k   = self._slope_var.get()

        ax0 = self.axes[0]            # Cre cutoff as a line x = ch1 + k*y
        y0, y1 = ax0.get_ylim()
        self._lines_ch1[0].set_xdata([ch1 + k * y0, ch1 + k * y1])
        self._lines_ch1[0].set_ydata([y0, y1])
        self._lines_ch2[0].set_xdata(ax0.get_xlim())
        self._lines_ch2[0].set_ydata([ch2, ch2])

        ax1 = self.axes[1]            # with a slope, show the cutoff at the median counter value
        ctr = self._metric_values(self.counter_col, self._metric_ctr.get())
        c1_mid = ch1 + k * float(np.nanmedian(ctr)) if k else ch1
        self._lines_ch1[1].set_xdata(ax1.get_xlim())
        self._lines_ch1[1].set_ydata([c1_mid, c1_mid])

        ax2 = self.axes[2]
        self._lines_ch2[2].set_xdata(ax2.get_xlim())
        self._lines_ch2[2].set_ydata([ch2, ch2])

        self.canvas.draw_idle()

    # ── Export ─────────────────────────────────────────────────────────────────

    def _export(self):
        if self.df is None:
            messagebox.showinfo('No data', 'Load a file first.')
            return

        labels   = self._classify()
        ch1      = self.ch1_var.get()
        ch2      = self.ch2_var.get()
        cre_name = self._short_name(self.cre_col)
        ctr_name = self._short_name(self.counter_col)

        # Standardised classification values (used by group_analysis_gui.py)
        cls_map = {_CRE_ONLY: 'Cre+', _CTR_ONLY: 'Cre-', _BOTH: 'Both positive'}
        # User-defined genotype group names
        grp_map = {
            _CRE_ONLY: self._name_cre.get()  or 'Cre+',
            _CTR_ONLY: self._name_ctr.get()  or 'Cre-',
            _BOTH:     self._name_both.get() or 'Both positive',
        }

        # Labels are positionally aligned to self.df rows — assign them directly
        # (soma numbers are not unique across files, so a soma#→label dict would
        # collide).
        cls_vals = [cls_map.get(l, _UNDET) for l in labels]
        # Each file's own Cre-channel name (GFP-AKT1-Cre vs GFP-AKT3-Cre ...):
        # the Group names are templates built from the header name, so swap in
        # the file's name, and export it as its own column so Group Analysis
        # can keep cells with different Cre constructs in separate groups.
        cre_files = [self._cre_name_for_file(f)
                     for f in self.df[self.file_col].astype(str)]

        def _grp(label, cre_f):
            g = grp_map.get(label, _UNDET)
            if label != _UNDET and cre_name and cre_f and cre_name in g:
                g = g.replace(cre_name, cre_f)
            return g
        grp_vals = [_grp(l, c) for l, c in zip(labels, cre_files)]

        out = self.df.copy()
        # Preserve incoming animal IDs (plain text) before we reposition columns.
        _aid_col = next((c for c in out.columns
                         if str(c).strip().lower().replace(' ', '') == 'animalid'),
                        None)
        animal_series = (out[_aid_col] if _aid_col is not None
                         else pd.Series('', index=out.index))
        drop = [c for c in out.columns
                if ('Channel' in c and 'detected' in c.lower())
                or ('Source'  in c and 'channel'  in c.lower())
                or str(c).strip() == 'Genotype'          # re-inserted below
                or (_aid_col is not None and c == _aid_col)]
        out = out.drop(columns=drop, errors='ignore')

        try:
            file_col = next(c for c in out.columns if c.lower() == 'file')
            file_idx = list(out.columns).index(file_col)
        except StopIteration:
            file_col = out.columns[0]
            file_idx = 0
        # Order: File | AnimalID | Genotype | Classification | Group | …
        geno_plain = out[file_col].map(
            lambda f: self._geno_plain_by_file.get(str(f), ''))
        out.insert(file_idx + 1, 'AnimalID',       animal_series.values)
        out.insert(file_idx + 2, 'Genotype',       geno_plain)
        out.insert(file_idx + 3, 'Classification', cls_vals)
        out.insert(file_idx + 4, 'Group',          grp_vals)
        out.insert(file_idx + 5, 'Cre channel',    cre_files)
        keys = self._cell_keys()
        out.insert(file_idx + 6, 'Override',
                   [self._overrides.get(k_, '') for k_ in keys])

        n_before = len(out)
        out      = out[out['Classification'] != _UNDET].copy()
        n_excl   = n_before - len(out)

        n_cre  = (out['Classification'] == 'Cre+').sum()
        n_ctr  = (out['Classification'] == 'Cre-').sum()
        n_both = (out['Classification'] == 'Both positive').sum()

        mean_cols = [c for c in out.columns if str(c).startswith('Mean ')]
        ch_names  = [c[5:] for c in mean_cols]

        out_path = os.path.join(self.out_dir, 'classified_' + self.in_fname)
        try:
            cre_m, ctr_m = self._metric_cre.get(), self._metric_ctr.get()
            fmt = lambda m, v: ('{:.2f}' if m == 'F/F0' else '{:.1f}').format(v)
            settings = [('Source workbook', self.in_fname),
                        ('Cre channel', cre_name), ('Cre metric', METRIC_LABEL[cre_m]),
                        ('Cre cutoff (>=)', ch1),
                        ('Counter channel', ctr_name), ('Counter metric', METRIC_LABEL[ctr_m]),
                        ('Counter cutoff (>=)', ch2),
                        ('Bleed-through slope (Cre cutoff += slope x counter)', float(self._slope_var.get())),
                        ('Manual overrides (Override column)', int(sum(1 for k_ in keys if k_ in self._overrides))),
                        ('Cells exported', len(out)), ('Excluded (below both cutoffs)', n_excl),
                        ('Exported', time.strftime('%Y-%m-%d %H:%M'))]
            self._write_formatted_xlsx(out, out_path, ch_names,
                                       geno_rt_by_file=self._geno_rt_by_file,
                                       file_col=file_col, settings=settings)

            self.lbl_status.config(
                text='Exported {:,} cells\n'
                     'Cre+ ({}): {:,}\n'
                     'Cre- ({}): {:,}\n'
                     'Both ({}): {:,}\n'
                     'Excluded: {:,}\n'
                     'Cre ({}) {} >= {}\n'
                     'Counter ({}) {} >= {}'.format(
                         len(out),
                         grp_map[_CRE_ONLY], n_cre,
                         grp_map[_CTR_ONLY], n_ctr,
                         grp_map[_BOTH],     n_both,
                         n_excl,
                         cre_name, cre_m, fmt(cre_m, ch1),
                         ctr_name, ctr_m, fmt(ctr_m, ch2)))
            messagebox.showinfo(
                'Exported',
                'classified_' + self.in_fname + '\n\n'
                'Cre+ ({}): {:,}\nCre- ({}): {:,}\nBoth ({}): {:,}\n'
                'Excluded (unclassified): {:,}\n\n'
                'Cre channel ({}): {} >= {}\n'
                'Counter channel ({}): {} >= {}\n\n'
                'Metric and cutoffs are recorded on the sheet '
                '"Classification settings".'.format(
                    grp_map[_CRE_ONLY], n_cre,
                    grp_map[_CTR_ONLY], n_ctr,
                    grp_map[_BOTH],     n_both,
                    n_excl,
                    cre_name, cre_m, fmt(cre_m, ch1),
                    ctr_name, ctr_m, fmt(ctr_m, ch2)))
        except PermissionError:
            messagebox.showerror(
                'Export error',
                f'Cannot write to:\n{out_path}\n\n'
                'The file is open in Excel. Please close it and try again.')
        except Exception as e:
            messagebox.showerror('Export error', str(e))

    def _write_formatted_xlsx(self, df, path, ch_names,
                              geno_rt_by_file=None, file_col='File', settings=None):
        geno_rt_by_file = geno_rt_by_file or {}
        hdr_font = Font(name='Arial', bold=True, color='FFFFFF', size=11)
        sub_font = Font(name='Arial', bold=True, size=10)
        dat_font = Font(name='Arial', size=10)

        fills = {
            'meta': PatternFill('solid', start_color='4A4A4A'),
            'ch1':  PatternFill('solid', start_color='1F7D1F'),
            'ch2':  PatternFill('solid', start_color='B22222'),
            'ch3':  PatternFill('solid', start_color='1A5276'),
            'sub':  PatternFill('solid', start_color='D9D9D9'),
        }
        ch_fkeys = ['ch1', 'ch2', 'ch3']

        ctr  = Alignment(horizontal='center', vertical='center', wrap_text=True)
        thin = Side(style='thin',   color='AAAAAA')
        med  = Side(style='medium', color='555555')
        brd  = Border(left=thin, right=thin, top=thin, bottom=thin)
        brdR = Border(left=thin, right=med,  top=thin, bottom=thin)

        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = 'Cells'

        cols   = [str(c).strip() for c in df.columns]
        n_chs  = len(ch_names)
        # Info columns = everything before the first channel (Mean ...) column.
        n_info = next((i for i, c in enumerate(cols) if c.startswith('Mean ')),
                      len(cols))
        name_to_col = {c: i for i, c in enumerate(cols, 1)}
        geno_col_idx = name_to_col.get('Genotype')
        file_col_idx = name_to_col.get(str(file_col).strip()) \
            or name_to_col.get('File') or 1

        ws.merge_cells('A1:{}1'.format(get_column_letter(n_info)))
        c = ws['A1']
        c.value = 'Cell Info'; c.font = hdr_font
        c.fill = fills['meta']; c.alignment = ctr

        for i, ch_name in enumerate(ch_names):
            col_s = n_info + 1 + i * 4
            col_e = col_s + 3
            span  = '{}1:{}1'.format(get_column_letter(col_s),
                                     get_column_letter(col_e))
            ws.merge_cells(span)
            c = ws[span.split(':')[0]]
            c.value = 'Ch{} - {}'.format(i + 1, ch_name)
            c.font  = hdr_font
            c.fill  = fills[ch_fkeys[i]] if i < 3 else fills['ch3']
            c.alignment = ctr

        for ci, lbl in enumerate(cols, 1):
            c = ws.cell(row=2, column=ci, value=lbl)
            c.font = sub_font; c.fill = fills['sub']
            c.alignment = ctr; c.border = brd

        def _info_width(name):
            n = name.lower()
            if n == 'file':                     return 34
            if n == 'animalid':                 return 20
            if n == 'genotype':                 return 22
            if n in ('classification', 'group'): return 14
            if n == 'cre channel':              return 16
            if 'soma' in n:                     return 8
            if 'area' in n:                     return 11
            if n == 'slice':                    return 7
            if 'centroid' in n:                 return 9
            return 12
        base_widths = [_info_width(c) for c in cols[:n_info]]
        ch_widths   = [11, 11, 9, 10] * n_chs
        for ci, w in enumerate(base_widths + ch_widths, 1):
            ws.column_dimensions[get_column_letter(ci)].width = w

        ws.row_dimensions[1].height = 22
        ws.row_dimensions[2].height = 34
        ws.freeze_panes = 'A3'

        # Special columns resolved by header name (layout-independent).
        _centroid_cols = {name_to_col[c] for c in cols if 'Centroid' in c}
        _brdR_cols     = set(_centroid_cols)
        if 'Slice' in name_to_col:    _brdR_cols.add(name_to_col['Slice'])
        if geno_col_idx:              _brdR_cols.add(geno_col_idx)
        _ff0_cols   = set()
        _range_cols = set()
        for i in range(n_chs):
            base = n_info + 1 + i * 4
            _ff0_cols.add(base + 2)
            _range_cols.add(base + 3)
            _brdR_cols.add(base + 3)

        for ri, row_vals in enumerate(df.itertuples(index=False), 3):
            row_file = str(row_vals[file_col_idx - 1])
            for ci, val in enumerate(row_vals, 1):
                if ci == geno_col_idx:
                    # Write the preserved rich-text genotype for this file.
                    # Copy the CellRichText so the same object isn't shared
                    # across cells.
                    rt = geno_rt_by_file.get(row_file)
                    if isinstance(rt, CellRichText):
                        cell_val = CellRichText(list(rt))
                    elif rt is not None:
                        cell_val = rt
                    else:
                        cell_val = val if isinstance(val, str) and val else None
                    c = ws.cell(row=ri, column=ci, value=cell_val)
                    c.font      = dat_font
                    c.alignment = Alignment(horizontal='center', vertical='center')
                    c.border    = brdR if ci in _brdR_cols else brd
                    continue
                if not isinstance(val, str) and pd.isna(val):
                    val = None
                c = ws.cell(row=ri, column=ci, value=val)
                c.font      = dat_font
                c.alignment = Alignment(horizontal='center', vertical='center')
                c.border    = brdR if ci in _brdR_cols else brd
                if isinstance(val, float):
                    if   ci in _ff0_cols:      c.number_format = '0.00'
                    elif ci in _range_cols:    c.number_format = '0.0'
                    elif ci in _centroid_cols: c.number_format = '0.00'
                    else:                      c.number_format = '0.000'

        if settings:
            # Which metric and cutoffs produced the Classification column.
            ws3 = wb.create_sheet('Classification settings')
            ws3.column_dimensions['A'].width = 30
            ws3.column_dimensions['B'].width = 44
            for ri, (k, v) in enumerate(settings, 1):
                ws3.cell(row=ri, column=1, value=k).font = sub_font
                ws3.cell(row=ri, column=2, value=v).font = dat_font
        wb.save(path)


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    root = tk.Tk()
    CreClassifierApp(root)
    root.mainloop()
