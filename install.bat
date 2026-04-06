@echo off
REM DPT2 Setup Script for Windows
REM Usage: install.bat

echo ==========================================
echo DPT2 Installation Script (Windows)
echo ==========================================
echo.

REM Check Python
echo Checking Python installation...
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found
    echo Please install Python 3.9 or higher from https://www.python.org/
    pause
    exit /b 1
)

for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo [OK] Found Python %PYTHON_VERSION%

REM Create virtual environment
echo.
echo Creating virtual environment...
if not exist "venv" (
    python -m venv venv
    echo [OK] Virtual environment created
) else (
    echo [OK] Virtual environment already exists
)

REM Activate virtual environment
echo Activating virtual environment...
call venv\Scripts\activate.bat
if errorlevel 1 (
    echo [ERROR] Failed to activate virtual environment
    pause
    exit /b 1
)
echo [OK] Virtual environment activated

REM Install DPT2
echo.
echo Installing DPT2 and dependencies...
python -m pip install -e .
if errorlevel 1 (
    echo [ERROR] Installation failed
    pause
    exit /b 1
)
echo [OK] Installation successful

REM Check for BLAST
echo.
echo Checking for BLAST+ ^(optional^)...
where blastn >nul 2>&1
if errorlevel 1 (
    echo [WARNING] BLAST+ not found ^(optional - needed only for --blast-screen^)
    echo Download from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
) else (
    for /f "tokens=*" %%i in ('blastn -version 2^>^&1 ^| findstr /N ".*"') do (
        set BLAST_INFO=%%i
        goto :found_blast
    )
    :found_blast
    echo [OK] BLAST+ found
)

REM Summary
echo.
echo ==========================================
echo Setup Complete!
echo ==========================================
echo.
echo To activate the environment in future sessions:
echo   venv\Scripts\activate.bat
echo.
echo To design primers:
echo   python scripts\design_and_score_primers.py --genes RPP30 MRGPRX1
echo.
echo For more help:
echo   type INSTALL.md
echo.
pause
