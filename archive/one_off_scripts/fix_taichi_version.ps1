# Fix Taichi Version for RTX 5070 Compatibility
# Taichi 1.7+ has memory allocation bugs with new GPUs

Write-Host "================================" -ForegroundColor Cyan
Write-Host "Taichi Version Fix" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""

Write-Host "Problem: Taichi 1.7+ has GPU memory allocation issues" -ForegroundColor Yellow
Write-Host "Solution: Downgrade to stable Taichi 1.6.0" -ForegroundColor Green
Write-Host ""

# Check current version
Write-Host "Current Taichi version:" -ForegroundColor White
python -c "import taichi as ti; print(f'  {ti.__version__}')" 2>$null

Write-Host ""
Write-Host "Uninstalling current Taichi..." -ForegroundColor Yellow
pip uninstall taichi -y

Write-Host ""
Write-Host "Installing Taichi 1.6.0 (stable)..." -ForegroundColor Yellow
pip install taichi==1.6.0

Write-Host ""
Write-Host "Verifying installation..." -ForegroundColor White
python -c "import taichi as ti; print(f'Installed: Taichi {ti.__version__}')"

Write-Host ""
Write-Host "================================" -ForegroundColor Cyan
Write-Host "Complete!" -ForegroundColor Green
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Now you can run:" -ForegroundColor White
Write-Host "  .\run_hybrid_test.ps1" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press any key to exit..." -ForegroundColor Gray
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")

