@echo off
REM Push this local repository to GitHub.
REM
REM BEFORE running: create an empty PUBLIC repository on github.com
REM (do NOT tick "Add a README", "Add .gitignore" or "Choose a license").
REM Then run this file and paste the repository URL when asked.
REM
REM If GitHub asks for a password, it will NOT accept your account password.
REM Use a Personal Access Token instead:
REM   GitHub -> Settings -> Developer settings -> Personal access tokens
REM   -> Tokens (classic) -> Generate new token -> tick "repo"

cd /d "%~dp0"

set /p REPO=Repository URL (example: https://github.com/USER/hearing-loss-cross-endpoint-scDRS.git):
if "%REPO%"=="" (
  echo No URL entered. Nothing was pushed.
  pause
  exit /b 1
)

echo.
echo Remote: %REPO%
echo.

git remote remove origin 2>nul
git remote add origin %REPO%
if errorlevel 1 (
  echo Failed to add remote.
  pause
  exit /b 1
)

git push -u origin main
if errorlevel 1 (
  echo.
  echo Push failed. Common causes:
  echo   1. Repository does not exist yet, or URL typo
  echo   2. Authentication rejected - use a Personal Access Token, not your password
  echo   3. Repository is private - Zenodo needs it to be PUBLIC
) else (
  echo.
  echo Push succeeded. Open the repository page on github.com to confirm 36 files are there.
)

pause
