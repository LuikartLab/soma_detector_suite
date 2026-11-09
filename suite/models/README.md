# Model files

| file | what | where to get it |
|---|---|---|
| `soma_pick_model.joblib` (0.3 MB, included) | the human-pick ranking model (scikit-learn HistGradientBoosting; 47,097 tracks from 76 curated stacks) | this folder |
| `cellpose/cpsam_v2` (1.2 GB) | the published Cellpose-SAM weights | Cellpose downloads them on first use into `~/.cellpose/models/`; the portable bundle ships a copy in `models/cellpose/` and `cellpose_engine.py` points Cellpose there when that folder exists |
| `cpsam_luikart_v2` (1.2 GB) | Cellpose-SAM fine-tuned on the laboratory's hand-traced somata (2026-09-23; see the paper) — the default model of v3.0 | **[release asset / Zenodo DOI — to be added]**; place it in this folder. Without it the engine falls back to `cpsam_v2` and the window title says so (and the outline-breadth default returns to −0.5 µm) |

`cellpose_engine.py` selects the weights with the `CUSTOM_MODEL` constant near its top.
