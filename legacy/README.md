# Legacy: the classical detector (v2.1)

`soma_area_and_intensity_gui_2.1.py` is the last version of the threshold-and-watershed detector that preceded
the Cellpose-SAM engine (described in the methods paper as "the classical detector"). It needs no GPU and writes
the same workbooks and ROI archives as v3.0, so the companion tools in `../suite` work with its output.
`Soma_Detector_2.1.bat` starts it together with the current Group Analysis and Cre classifier from `../suite`;
it runs from wherever this folder is placed and uses the same Python search order as the v3.0 launcher
(`python312\` or `.venv-cellpose\` inside `../suite`, else `pythonw` on PATH). `cell_jump_ipc.py` is a copy of
the suite's module, which the 2.1 GUI imports.
