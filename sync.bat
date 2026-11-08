@echo off
rem Copies code, raw data and results between this computer and the lab server.
rem
rem   sync pull    get the latest code, raw data and results
rem   sync push    send this computer's results to the server
rem
rem Add --dry-run to either one to list what would be copied without copying.
rem The raw data is mirrored from the server, so local raw files that are not on
rem the server are deleted. Results are copied in both directions without deleting
rem anything, and a file is replaced only when the incoming copy is newer.

setlocal
set "SERVER=R:\Basic_Sciences\CDB\OzbudakLab\Arya CDB\MEK_oscillation_project"
set "RESULTS=analysis_scripts\ihc_results"
for %%I in ("%~dp0..") do set "LOCAL=%%~fI"

if not exist "%SERVER%\" (
    echo Cannot reach the server folder: %SERVER%
    exit /b 1
)
if /I "%LOCAL%"=="%SERVER%" (
    echo This copy is on the server. Run sync from a copy on this computer.
    exit /b 1
)

set "DRY="
if /I "%~2"=="--dry-run" set "DRY=/L"
set "OPTS=/MT:16 /NP /NDL /NJH %DRY%"
set "SKIP=/XD models runtime viewer_settings"

if /I "%~1"=="pull" goto pull
if /I "%~1"=="push" goto push
echo Usage: sync pull ^| sync push   (add --dry-run to only list changes)
exit /b 1

:pull
if not defined DRY (
    git -C "%~dp0." pull --ff-only
    if errorlevel 1 exit /b 1
)
echo Raw data:
robocopy "%SERVER%" "%LOCAL%" /MIR %OPTS% /XD analysis_scripts analysis_scripts_backup_* .venv /XF *.ims
if errorlevel 8 exit /b 1
echo Results:
robocopy "%SERVER%\%RESULTS%" "%LOCAL%\%RESULTS%" /E /XO %OPTS% %SKIP%
if errorlevel 8 exit /b 1
exit /b 0

:push
echo Results:
robocopy "%LOCAL%\%RESULTS%" "%SERVER%\%RESULTS%" /E /XO %OPTS% %SKIP%
if errorlevel 8 exit /b 1
echo Code changes not yet committed:
git -C "%~dp0." status --short
exit /b 0
