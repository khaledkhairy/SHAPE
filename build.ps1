# Build everything on Windows: venv, C++ core extension, and the double-clickable exe.
# Usage (from the project root, in PowerShell):
#     .\build.ps1
$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
Set-Location $root

if (-not (Test-Path ".venv")) {
    Write-Host "Creating virtual environment..." -ForegroundColor Cyan
    python -m venv .venv
}
$py = ".\.venv\Scripts\python.exe"

Write-Host "Installing Python dependencies..." -ForegroundColor Cyan
& $py -m pip install --upgrade pip
& $py -m pip install -r requirements.txt

Write-Host "Building the C++ core (shp_core) with MSVC..." -ForegroundColor Cyan
& $py setup.py build_ext --inplace

Write-Host "Packaging the double-clickable application..." -ForegroundColor Cyan
# stop a previous instance so the dist folder is not locked
Get-Process SHAPE -ErrorAction SilentlyContinue | Stop-Process -Force
& .\.venv\Scripts\pyinstaller.exe shape.spec --noconfirm --distpath dist --workpath build_pyi

Write-Host ""
Write-Host "Done. Double-click:  dist\SHAPE\SHAPE.exe" -ForegroundColor Green
