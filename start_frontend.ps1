# Skrypt do uruchamiania frontend  
# Użycie: .\start_frontend.ps1

Write-Host "Uruchamianie frontend Live 2.0..." -ForegroundColor Green

# Sprawdź czy node_modules istnieje
if (-not (Test-Path "frontend\node_modules")) {
    Write-Host "Instalowanie zależności frontend..." -ForegroundColor Yellow
    Set-Location frontend
    npm install
    Set-Location ..
}

# Sprawdź czy package.json istnieje
if (-not (Test-Path "frontend\package.json")) {
    Write-Host "BAD Nie znaleziono frontend/package.json!" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "Uruchamianie serwera frontend..." -ForegroundColor Cyan
Write-Host "Lokalizacja: http://localhost:3000" -ForegroundColor White
Write-Host ""
Write-Host "Aby zatrzymać serwer: Ctrl+C" -ForegroundColor Yellow
Write-Host "----------------------------------------" -ForegroundColor Gray

try {
    Set-Location frontend
    npm run dev
}
catch {
    Write-Host "BAD Błąd podczas uruchamiania frontend:" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Rozwiązanie:" -ForegroundColor Cyan
    Write-Host "1. Sprawdź czy port 3000 nie jest zajęty" -ForegroundColor White
    Write-Host "2. Uruchom: netstat -ano | findstr :3000" -ForegroundColor White
    Write-Host "3. Spróbuj: npm run build && npm run preview" -ForegroundColor White
} finally {
    Set-Location ..
}

Write-Host "Serwer frontend zatrzymany" -ForegroundColor Green
