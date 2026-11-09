@echo off
rem Classical (CPU-only) Soma Detector v2.1 with the current Group Analysis and Cre
rem classifier from ..\suite.  Runs from wherever this folder lives.  Python search
rem order (same as the v3.0 launcher): python312\ or .venv-cellpose\ inside ..\suite,
rem then the system pythonw.
set "HERE=%~dp0"
set "SUITE=%HERE%..\suite\"
if exist "%SUITE%python312\pythonw.exe" (
    set "PY=%SUITE%python312\pythonw.exe"
) else if exist "%SUITE%.venv-cellpose\Scripts\pythonw.exe" (
    set "PY=%SUITE%.venv-cellpose\Scripts\pythonw.exe"
) else (
    set "PY=pythonw"
)
start "Soma Detector 2.1" "%PY%" "%HERE%soma_area_and_intensity_gui_2.1.py"
timeout /t 2 /nobreak >nul
start "Group Analysis" "%PY%" "%SUITE%group_analysis_gui.py"
start "Cre Classifier" "%PY%" "%SUITE%cre_classifier_gui.py"
