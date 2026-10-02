@echo off
REM Open Git Bash in THIS folder.
REM This machine has no standard Git for Windows, so there is no
REM "Git Bash Here" right-click entry. This launcher opens the portable
REM Git Bash that ships with WorkBuddy, already cd'd into this folder.
cd /d "%~dp0"
"C:\Users\docta\.workbuddy\binaries\PortableGit\versions\1.2.0\bin\bash.exe" --login -i
