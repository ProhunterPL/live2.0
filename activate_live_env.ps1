# Skrypt do aktywacji środowiska conda 'live'
# Użycie: .\activate_live_env.ps1

[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

Write-Host "Aktywacja środowiska conda 'live'..." -ForegroundColor Green

# Sprawdź czy conda jest dostępna
if (-not (Get-Command "D:\conda\Scripts\conda.exe" -ErrorAction SilentlyContinue)) {
    Write-Host "BAD Conda nie znaleziona w oczekiwanej lokalizacji" -ForegroundColor Red
    Write-Host "Sprawdź czy miniconda jest zainstalowana w D:\conda" -ForegroundColor Yellow
    exit 1
}

# Zainicjalizuj conda dla PowerShell
Write-Host "Inicjalizacja conda..." -ForegroundColor Yellow
D:\conda\Scripts\conda.exe init powershell

# Ustaw zmienne środowiskowe dla bieżącej sesji
Write-Host "Konfiguracja zmiennych środowiskowych..." -ForegroundColor Yellow

# Dodaj conda do PATH
$condaScriptsPath = "D:\conda\Scripts"
$condaEnvsScriptsPath = "D:\conda_envs\live\Scripts"

if ($env:PATH -notlike "*$condaScriptsPath*") {
    $env:PATH = "$condaScriptsPath;$env:PATH"
}

if ($env:PATH -notlike "*$condaEnvsScriptsPath*") {
    $env:PATH = "$condaEnvsScriptsPath;$env:PATH"
}

# Ustaw CONDA_PREFIX
$env:CONDA_PREFIX = "D:\conda_envs\live"

# Ustaw CONDA_PROMPT_MODIFIER
$env:CONDA_PROMPT_MODIFIER = "(live) "

# Ustaw PYTHONPATH na środowisko conda
$env:PYTHONPATH = "D:\conda_envs\live\Lib\site-packages"

Write-Host "OK Środowisko 'live' zostalo skonfigurowane!" -ForegroundColor Green
Write-Host ""
Write-Host "Dostepne komendy:" -ForegroundColor Cyan
Write-Host "  python -> D:\conda_envs\live\python.exe" -ForegroundColor White
Write-Host "  pip    -> D:\conda_envs\live\Scripts\pip.exe" -ForegroundColor White
Write-Host ""
Write-Host "Uruchomienie aplikacji:" -ForegroundColor Cyan
Write-Host "  Backend:  D:\conda_envs\live\python.exe -m backend.api.server" -ForegroundColor White
Write-Host "  Frontend: npm run dev (w katalogu frontend)" -ForegroundColor White
Write-Host ""
Write-Host 'Test instalacji:' -ForegroundColor Cyan
Write-Host ' - python -c "import taichi"' -ForegroundColor White
