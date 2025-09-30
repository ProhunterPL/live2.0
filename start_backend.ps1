# Prosty skrypt do uruchamiania backend Live 2.0

Write-Host "Uruchamianie backend Live 2.0..." -ForegroundColor Green

# Sprawdź czy port 8001 jest zajęty
$portCheck = netstat -ano | Select-String ":8001 "
if ($portCheck) {
    Write-Host "⚠️  Port 8001 jest już zajęty!" -ForegroundColor Yellow
    Write-Host "Uruchamianie skryptu czyszczącego..." -ForegroundColor Cyan
    & ".\cleanup_processes.ps1"
    
    # Poczekaj chwilę na czyszczenie
    Start-Sleep -Seconds 2
    
    # Sprawdź ponownie
    $portCheck2 = netstat -ano | Select-String ":8001 "
    if ($portCheck2) {
        Write-Host "❌ Port 8001 nadal jest zajęty. Uruchom ręcznie: .\cleanup_processes.ps1" -ForegroundColor Red
        exit 1
    }
}

# Sprawdź czy środowisko conda istnieje
if (-not (Test-Path "D:\conda_envs\live\python.exe")) {
    Write-Host "BAD Środowisko conda 'live' nie istnieje!" -ForegroundColor Red
    Write-Host "Uruchom najpierw: .\activate_live_env.ps1" -ForegroundColor Yellow
    exit 1
}

$python = "D:\conda_envs\live\python.exe"

# Przejdź do katalogu backend
Set-Location backend

# Ustaw PYTHONPATH
$env:PYTHONPATH = (Get-Location).Path

Write-Host ""
Write-Host "Uruchamianie serwera backend..." -ForegroundColor Cyan
Write-Host "Lokalizacja: http://localhost:8001" -ForegroundColor White
Write-Host "API Docs: http://localhost:8001/docs" -ForegroundColor White
Write-Host ""
Write-Host "Aby zatrzymać serwer: Ctrl+C" -ForegroundColor Yellow
Write-Host "------------------------------" -ForegroundColor Gray

try {
    & $python api/server.py
}
catch {
    Write-Host "BAD Błąd podczas uruchamiania serwera:" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Yellow
}
finally {
    # Powrót do katalogu głównego
    Set-Location ..
}

Write-Host "Serwer zatrzymany" -ForegroundColor Green
