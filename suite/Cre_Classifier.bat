@echo off
rem Cre classifier launcher (Soma Detector Suite v3.0).  Same Python search as
rem Soma_Detector_3.0.bat: python312\ (portable), .venv-cellpose\, system pythonw.
set "HERE=%~dp0"
if exist "%HERE%python312\pythonw.exe" (
    set "PY=%HERE%python312\pythonw.exe"
) else if exist "%HERE%.venv-cellpose\Scripts\pythonw.exe" (
    set "PY=%HERE%.venv-cellpose\Scripts\pythonw.exe"
) else (
    set "PY=pythonw"
)
start "Cre Classifier" "%PY%" "%HERE%cre_classifier_gui.py"
