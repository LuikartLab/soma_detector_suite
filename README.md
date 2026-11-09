# Soma Detector Suite v3.0

Automated three-dimensional soma detection, interactive curation, fluorescence quantification and group-level
statistics for multi-channel confocal stacks (Zeiss CZI / LSM) of genetically manipulated neurons — built by a
neuroscience laboratory in collaboration with AI coding assistants, and validated against the laboratory's
published hand-traced data (Prina, Goyette et al., *Brain* 2026). The methods paper that describes and validates
the suite is published separately; this repository holds only what is needed to run the programs.

The detection engine is **Cellpose-SAM** (Pachitariu, Rariden & Stringer 2025), fine-tuned on the laboratory's
hand-traced somata; the suite adds z-linking of per-plane labels, plane selection, a learned ranking of
candidates for curation, a curation interface, Cre-phenotype classification, and a group-analysis tool that
fits the mixed-effects model of the published study.

## What is here

| folder | contents |
|---|---|
| [`suite/`](suite/) | the programs: `soma_area_and_intensity_gui_3.0.py` (detector + curation), `cellpose_engine.py`, `cre_classifier_gui.py`, `group_analysis_gui.py`, `mixed_stats.py` (the statistical model), helpers and launchers; `models/` (ranking model; see its README for the Cellpose weights) |
| [`docs/`](docs/) | the portable-bundle install notes and the v2.2 user guide |
| [`legacy/`](legacy/) | the classical (v2.1) detector, CPU-only, with its own launcher |

The repository can be placed in any folder. The launchers run from their own location, the programs find each
other and the model files relative to their own file, log to the system temp folder, and take image data and
workbooks through file dialogs; nothing refers to a fixed path.

## Requirements

* Windows 10/11 (the interfaces are Tkinter; the engine itself is portable).
* An NVIDIA GPU is strongly recommended: detection of a 30-plane two-channel stack takes ~18 s on an RTX 5080
  and ~50× longer on CPU. 8 GB of GPU memory suffices for detection; fine-tuning needs 16 GB.
* Python 3.12 with the packages in [`requirements.txt`](requirements.txt) (PyTorch CUDA 12.8 build, Cellpose 4.2).

## Install

**From source** — create the environment inside `suite/`, where the launchers look for it:
```bat
cd suite
python -m venv .venv-cellpose
.venv-cellpose\Scripts\activate
pip install -r ..\requirements.txt --extra-index-url https://download.pytorch.org/whl/cu128
```
Put the fine-tuned weights `cpsam_luikart_v2` into `suite/models/` (see [`suite/models/README.md`](suite/models/README.md));
Cellpose downloads the base weights on first use. Then double-click `suite\Soma_Detector_3.0.bat` — it starts the
detector, Group Analysis and the Cre classifier with `python312\` or `.venv-cellpose\` found next to itself, else
the `pythonw` on PATH — or run `python suite\soma_area_and_intensity_gui_3.0.py` with the environment active.

**Portable bundle (no Python install)** — the laboratory distributes a self-contained folder (private Python
with all packages and both models, ~5 GB zipped); see [`docs/INSTALL_portable_bundle.txt`](docs/INSTALL_portable_bundle.txt).

## Workflow

1. **Soma Detector** — open a stack, confirm channels / animal / genotype in the file-info dialog, detect, curate
   (click to reject, hold F and click to force-detect, edit outlines), save the ROI set (`*_ROIs.npz`),
   then *Batch Process* a directory of curated stacks to the measurement workbooks.
2. **Cre classifier** — gate Cre-positive / Cre-negative / double-positive cells on raw or background-normalised
   (F/F₀) intensity, with per-cell overrides; writes `classified_*.xlsx` with a settings sheet.
3. **Group Analysis** — conditions from genotype and Cre construct, per-image Cre-negative normalisation within a
   z-window, and the published statistical model (random-intercept mixed model, cluster-robust SE, Scheffé contrasts).

## Validation

The suite was validated cell by cell against the hand-traced data of the published study (63 confocal stacks,
1,120 somata: area concordance CCC 0.94, median area ratio 0.99, no size dependence; normalised reporter
intensity r = 0.85, 0.95 on the best-curated set), and `mixed_stats.py` reproduces the study's published
mixed-effects contrasts exactly from its supplementary datasheets. The analysis is reported in the methods paper
(see Citation).

## Citation

See [`CITATION.cff`](CITATION.cff). Please also cite Cellpose-SAM (Pachitariu et al., 2025) and Cellpose
(Stringer et al., 2021), and the validation dataset's study (Prina, Goyette et al., *Brain* 2026,
doi:10.1093/brain/awag098).

## License

MIT — see [`LICENSE`](LICENSE).
