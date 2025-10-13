# Phase 2 Overnight Test Runner
# Uruchamia testową symulację Miller-Urey na noc
# ~10 godzin runtime dla 10M kroków

Write-Host "================================================" -ForegroundColor Cyan
Write-Host "PHASE 2 OVERNIGHT TEST - MILLER-UREY" -ForegroundColor Cyan
Write-Host "================================================" -ForegroundColor Cyan
Write-Host ""

$timestamp = Get-Date -Format "yyyy-MM-dd_HH-mm-ss"
$output_dir = "results/overnight_test_$timestamp"
$log_file = "$output_dir/simulation.log"

Write-Host "Konfiguracja:" -ForegroundColor Yellow
Write-Host "  Scenariusz: Miller-Urey (test scale)" -ForegroundColor White
Write-Host "  Kroki: 10,000,000" -ForegroundColor White
Write-Host "  Molekuły: 180 (650 atomów)" -ForegroundColor White
Write-Host "  Output: $output_dir" -ForegroundColor White
Write-Host "  Log: $log_file" -ForegroundColor White
Write-Host "  Przewidywany czas: ~10 godzin" -ForegroundColor White
Write-Host ""

# Utworz katalog output
New-Item -ItemType Directory -Force -Path $output_dir | Out-Null

Write-Host "Uruchamiam symulację..." -ForegroundColor Green
Write-Host ""
Write-Host "UWAGA: Symulacja będzie działać w tle!" -ForegroundColor Yellow
Write-Host "Sprawdź postęp: Get-Content $log_file -Tail 20 -Wait" -ForegroundColor Yellow
Write-Host ""

# Zapisz info o starcie
$start_info = @"
Overnight Test Started
======================
Start Time: $(Get-Date)
Configuration: configs/phase2_miller_urey_test.yaml
Steps: 10,000,000
Expected Duration: ~10 hours
Expected Completion: $(Get-Date).AddHours(10)

Progress Monitoring:
  Get-Content $log_file -Tail 20 -Wait

Stop if needed:
  Get-Process python | Where-Object {`$_.CommandLine -like "*run_phase2_full*"} | Stop-Process
"@

$start_info | Out-File "$output_dir/START_INFO.txt"

# Uruchom symulację w tle
Start-Process python -ArgumentList @(
    "scripts/run_phase2_full.py",
    "--config", "configs/phase2_miller_urey_test.yaml",
    "--output", $output_dir,
    "--steps", "10000000",
    "--seed", "42"
) -RedirectStandardOutput "$output_dir/stdout.log" -RedirectStandardError "$output_dir/stderr.log" -NoNewWindow

Write-Host "✓ Symulacja uruchomiona!" -ForegroundColor Green
Write-Host ""
Write-Host "Sprawdź logi:" -ForegroundColor Cyan
Write-Host "  Get-Content $output_dir\stdout.log -Tail 20 -Wait" -ForegroundColor White
Write-Host "  Get-Content $output_dir\stderr.log -Tail 20 -Wait" -ForegroundColor White
Write-Host ""
Write-Host "Zatrzymaj jeśli potrzeba:" -ForegroundColor Cyan
Write-Host "  Get-Process python | Stop-Process" -ForegroundColor White
Write-Host ""
Write-Host "Dobranoc! Sprawdzimy wyniki rano! 🌙" -ForegroundColor Green

