# Test Hybrid GPU+CPU Stepper
# Runs tests and benchmark to verify hybrid mode works

Write-Host "================================" -ForegroundColor Cyan
Write-Host "Hybrid GPU+CPU Mode - Test Suite" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""

# Check if virtual environment exists
if (Test-Path ".\live_env\Scripts\Activate.ps1") {
    Write-Host "Activating virtual environment..." -ForegroundColor Yellow
    & .\live_env\Scripts\Activate.ps1
} else {
    Write-Host "Warning: Virtual environment not found. Using system Python." -ForegroundColor Yellow
}

Write-Host ""
Write-Host "Step 1: Running functionality tests..." -ForegroundColor Green
Write-Host ""

python tests\test_hybrid_stepper.py

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "❌ Tests failed! Fix errors before running benchmark." -ForegroundColor Red
    Write-Host ""
    Write-Host "Press any key to exit..." -ForegroundColor Gray
    $null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
    exit 1
}

Write-Host ""
Write-Host "✅ Tests passed!" -ForegroundColor Green
Write-Host ""
Write-Host "Step 2: Running performance benchmark..." -ForegroundColor Green
Write-Host "This compares: Pure GPU vs Pure CPU vs Hybrid" -ForegroundColor White
Write-Host "Estimated time: 3-5 minutes" -ForegroundColor Yellow
Write-Host ""

python tests\benchmark_hybrid.py --steps 200

Write-Host ""
Write-Host "================================" -ForegroundColor Cyan
Write-Host "Complete!" -ForegroundColor Green
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Results saved to: tests\hybrid_benchmark_results.json" -ForegroundColor White
Write-Host ""
Write-Host "Next steps:" -ForegroundColor Yellow
Write-Host "1. Check benchmark results above" -ForegroundColor White
Write-Host "2. If Hybrid is faster, use HybridSimulationStepper in your code" -ForegroundColor White
Write-Host "3. See docs\HYBRID_GPU_CPU_GUIDE.md for usage instructions" -ForegroundColor White
Write-Host ""
Write-Host "Press any key to exit..." -ForegroundColor Gray
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")

