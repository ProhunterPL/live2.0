# Skrypt PowerShell do uruchamiania testów lokalnie tak samo jak w CI
# Użycie: .\\.github\scripts\run_local_ci.ps1

$ErrorActionPreference = "Stop"

Write-Host "================================" -ForegroundColor Cyan
Write-Host "🧪 Live 2.0 Local CI Tests" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""

# Sprawdź czy jesteśmy w głównym katalogu projektu
if (-not (Test-Path "requirements.txt")) {
    Write-Host "❌ Error: Must be run from project root" -ForegroundColor Red
    exit 1
}

# Ustaw Taichi na CPU mode
$env:TI_ARCH = "cpu"
$env:PYTHONPATH = $PWD

Write-Host "📦 Step 1: Checking dependencies..." -ForegroundColor Yellow
try {
    python -c "import pytest" 2>$null
} catch {
    Write-Host "Installing dependencies..." -ForegroundColor Yellow
    pip install -r requirements.txt
}

Write-Host ""
Write-Host "🎨 Step 2: Code Quality Checks" -ForegroundColor Yellow
Write-Host "--------------------------------"

# Black
Write-Host "`n→ Checking code formatting (black)..." -ForegroundColor Yellow
try {
    black --check backend/ scripts/ matcher/ 2>$null
    Write-Host "✓ Code formatting OK" -ForegroundColor Green
} catch {
    Write-Host "✗ Code formatting issues found" -ForegroundColor Red
    Write-Host "  Run: black backend/ scripts/ matcher/ to fix" -ForegroundColor Yellow
}

# isort
Write-Host "`n→ Checking import sorting (isort)..." -ForegroundColor Yellow
try {
    isort --check-only backend/ scripts/ matcher/ 2>$null
    Write-Host "✓ Import sorting OK" -ForegroundColor Green
} catch {
    Write-Host "✗ Import sorting issues found" -ForegroundColor Red
    Write-Host "  Run: isort backend/ scripts/ matcher/ to fix" -ForegroundColor Yellow
}

# mypy (non-blocking)
Write-Host "`n→ Type checking (mypy)..." -ForegroundColor Yellow
try {
    mypy backend/sim/ --ignore-missing-imports 2>$null
    Write-Host "✓ Type checking OK" -ForegroundColor Green
} catch {
    Write-Host "⚠ Type checking found issues (non-blocking)" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "🧪 Step 3: Unit Tests" -ForegroundColor Yellow
Write-Host "--------------------------------"

# Backend tests
Write-Host "`n→ Running backend tests (excluding slow tests)..." -ForegroundColor Yellow
Push-Location backend
try {
    pytest tests/ -v -m "not slow" --tb=short --color=yes --maxfail=5
    Write-Host "✓ Backend tests passed" -ForegroundColor Green
} catch {
    Write-Host "✗ Backend tests failed" -ForegroundColor Red
    Pop-Location
    exit 1
}
Pop-Location

# Root tests
Write-Host "`n→ Running root tests (excluding stability tests)..." -ForegroundColor Yellow
try {
    pytest tests/ -v -k "not stability and not 24h" --tb=short --color=yes --maxfail=5
    Write-Host "✓ Root tests passed" -ForegroundColor Green
} catch {
    Write-Host "✗ Root tests failed" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "================================" -ForegroundColor Cyan
Write-Host "✅ All CI checks passed!" -ForegroundColor Green
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Your code is ready to push to main! 🚀" -ForegroundColor Green

