@echo off
rem Soma Detector Suite v3.0 launcher.  Runs from wherever this folder lives.
rem Python search order: python312\ (portable, shipped with the program),
rem .venv-cellpose\ (PyTorch CUDA + Cellpose-SAM), then the system pythonw.
set "HERE=%~dp0"
if exist "%HERE%python312\pythonw.exe" (
    set "PY=%HERE%python312\pythonw.exe"
) else if exist "%HERE%.venv-cellpose\Scripts\pythonw.exe" (
    set "PY=%HERE%.venv-cellpose\Scripts\pythonw.exe"
) else (
    set "PY=pythonw"
)
start "Soma Detector 3.0" "%PY%" "%HERE%soma_area_and_intensity_gui_3.0.py"
timeout /t 2 /nobreak >nul
start "Group Analysis" "%PY%" "%HERE%group_analysis_gui.py"
start "Cre Classifier" "%PY%" "%HERE%cre_classifier_gui.py"
