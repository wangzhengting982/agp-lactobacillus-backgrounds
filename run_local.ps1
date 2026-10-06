param([string]$PythonExe='python', [switch]$Setup)
$ErrorActionPreference='Stop'
$projectRoot=$PSScriptRoot
$venvPython=Join-Path $projectRoot '.venv/Scripts/python.exe'
if ($Setup -or -not (Test-Path -LiteralPath $venvPython)) {
    & $PythonExe -c "import sys; assert sys.version_info[:2]==(3,12), 'Python 3.12 required'"
    if ($LASTEXITCODE -ne 0) { throw 'Select a Python 3.12 interpreter with -PythonExe.' }
    & $PythonExe -m venv (Join-Path $projectRoot '.venv')
    if ($LASTEXITCODE -ne 0) { throw 'Virtual environment creation failed.' }
    & $venvPython -m pip install -r (Join-Path $projectRoot 'requirements-lock.txt')
    if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed.' }
}
& $venvPython -m pip check
if ($LASTEXITCODE -ne 0) { throw 'Dependency verification failed.' }
& $venvPython (Join-Path $projectRoot 'scripts/run_all.py')
if ($LASTEXITCODE -ne 0) { throw 'Reproduction failed. See reports/full_run_status.json and its log paths.' }
