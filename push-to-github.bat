@echo off
setlocal
REM Push this local repository to GitHub.
REM
REM BEFORE running: create an empty PUBLIC repository on github.com
REM (do NOT tick "Add a README", "Add .gitignore" or "Choose a license").
REM Then run this file and paste the repository URL when asked.
REM
REM NOTE: this machine has no standard Git for Windows, only the WorkBuddy
REM PortableGit build, which is NOT on the PATH of an Explorer-launched
REM cmd.exe. We therefore call git.exe by ABSOLUTE PATH and never rely on PATH.
REM
REM If GitHub asks for a password, it will NOT accept your account password.
REM Use a Personal Access Token instead:
REM   GitHub -> Settings -> Developer settings -> Personal access tokens
REM   -> Tokens (classic) -> Generate new token -> tick "repo"

cd /d "%~dp0"

set "GITROOT=C:\Users\docta\.workbuddy\binaries\PortableGit\versions\1.2.0"
set "GIT=%GITROOT%\cmd\git.exe"
if not exist "%GIT%" set "GIT=%GITROOT%\mingw64\bin\git.exe"
if not exist "%GIT%" (
  echo.
  echo ERROR: git.exe not found under %GITROOT%
  echo Expected one of:
  echo   %GITROOT%\cmd\git.exe
  echo   %GITROOT%\mingw64\bin\git.exe
  echo.
  pause
  exit /b 1
)

echo Using git: %GIT%
"%GIT%" --version
if errorlevel 1 (
  echo Git was found but failed to run.
  pause
  exit /b 1
)
echo.

set "REPO=%~1"
if "%REPO%"=="" set /p REPO=Repository URL (example: https://github.com/USER/hearing-loss-cross-endpoint-scDRS.git):
if "%REPO%"=="" (
  echo No URL entered. Nothing was pushed.
  pause
  exit /b 1
)

echo Remote: %REPO%
echo.

"%GIT%" remote remove origin 2>nul
"%GIT%" remote add origin "%REPO%"
if errorlevel 1 (
  echo.
  echo Failed to add remote. Check the URL.
  pause
  exit /b 1
)

"%GIT%" push -u origin main
if errorlevel 1 (
  echo.
  echo Push failed. Common causes:
  echo   1. Repository does not exist yet, or URL typo
  echo   2. Authentication rejected - use a Personal Access Token, not your password
  echo   3. Repository is private - Zenodo needs it to be PUBLIC
  echo.
  echo If a browser or GitHub login window opened, finish it, then run this file again.
) else (
  echo.
  echo Push succeeded. Open the repository page on github.com to confirm the files are there.
)

pause
endlocal
