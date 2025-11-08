# CPU Performance Benchmark
# Compares different CPU thread configurations

Write-Host "================================" -ForegroundColor Cyan
Write-Host "CPU Performance Benchmark" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Testing CPU with different thread counts" -ForegroundColor White
Write-Host "This avoids GPU memory issues" -ForegroundColor Yellow
Write-Host ""

# Check if virtual environment exists
if (Test-Path ".\live_env\Scripts\Activate.ps1") {
    Write-Host "Activating virtual environment..." -ForegroundColor Yellow
    & .\live_env\Scripts\Activate.ps1
}

python tests\benchmark_hybrid.py --modes pure_cpu --steps 100

Write-Host ""
Write-Host "================================" -ForegroundColor Cyan
Write-Host "Complete!" -ForegroundColor Green
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press any key to exit..." -ForegroundColor Gray
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")

