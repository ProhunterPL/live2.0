# Benchmark GPU vs CPU Performance
# Run this to test which backend is faster for your system

Write-Host "================================" -ForegroundColor Cyan
Write-Host "Live 2.0 - GPU vs CPU Benchmark" -ForegroundColor Cyan
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
Write-Host "Starting benchmark..." -ForegroundColor Green
Write-Host "This will test:" -ForegroundColor White
Write-Host "  1. GPU (CUDA) if available" -ForegroundColor White
Write-Host "  2. CPU with different thread counts" -ForegroundColor White
Write-Host ""
Write-Host "Estimated time: 5-10 minutes" -ForegroundColor Yellow
Write-Host ""

# Run benchmark
python tests\benchmark_gpu_vs_cpu.py

Write-Host ""
Write-Host "================================" -ForegroundColor Cyan
Write-Host "Benchmark Complete!" -ForegroundColor Green
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Results saved to: tests\benchmark_results.json" -ForegroundColor White
Write-Host ""
Write-Host "Press any key to exit..." -ForegroundColor Gray
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")

