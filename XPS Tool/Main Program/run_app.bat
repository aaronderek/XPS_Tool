@echo off
REM Quick start script for the 2DEG Visualization Tool on Windows
setlocal

REM Ensure the working directory is this script's folder
pushd "%~dp0"

echo =========================================
echo 2DEG Visualization Tool
echo =========================================
echo.
echo Starting Streamlit application...
echo The app will open in your browser at http://localhost:8501
echo.

py -3.11 -V >NUL 2>&1
if errorlevel 1 (
    echo [Error] Python 3.11 is not available on this system.
    echo Install Python 3.11 and rerun this script.
    goto :pause_and_exit
)

py -3.11 -c "import streamlit" >NUL 2>&1
if errorlevel 1 (
    echo [Error] Streamlit is not installed for Python 3.11.
    echo Run:  py -3.11 -m pip install -r requirements.txt
    echo and try again.
    goto :pause_and_exit
)

py -3.11 -m streamlit run app.py
if errorlevel 1 (
    echo.
    echo [Error] Failed to start the Streamlit application.
) else (
    echo.
    echo Streamlit has exited.
)

:pause_and_exit
echo.
pause
popd
endlocal
