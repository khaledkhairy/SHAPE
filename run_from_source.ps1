# Run the app directly from source (no packaging) for quick iteration.
#     .\run_from_source.ps1
$root = $PSScriptRoot
Set-Location $root
if (-not (Test-Path "shp_core*.pyd")) {
    Write-Host "Building the C++ core first..." -ForegroundColor Cyan
    & .\.venv\Scripts\python.exe setup.py build_ext --inplace
}
& .\.venv\Scripts\python.exe app\shape_app.py @args
