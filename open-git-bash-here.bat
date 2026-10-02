@echo off
setlocal
REM Open Git Bash in THIS folder.
REM This machine has no standard Git for Windows, so there is no
REM "Git Bash Here" right-click entry. This launcher opens the portable
REM Git Bash that ships with WorkBuddy, already cd'd into this folder.

cd /d "%~dp0"

set "GITROOT=C:\Users\docta\.workbuddy\binaries\PortableGit\versions\1.2.0"
set "BASH=%GITROOT%\bin\bash.exe"
if not exist "%BASH%" set "BASH=%GITROOT%\usr\bin\bash.exe"
if not exist "%BASH%" set "BASH=%GITROOT%\git-bash.exe"
if not exist "%BASH%" (
  echo.
  echo ERROR: bash.exe not found under %GITROOT%
  echo.
  pause
  exit /b 1
)

"%BASH%" --login -i
endlocal
