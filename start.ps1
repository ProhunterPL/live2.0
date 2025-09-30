# Quick Start Script for Live 2.0
# Run this script to start the application

Write-Host "🚀 Starting Live 2.0 Simulation..." -ForegroundColor Green

# Check if Docker is running
try {
    docker version | Out-Null
    Write-Host "✅ Docker is running" -ForegroundColor Green
} catch {
    Write-Host "❌ Docker is not running. Please start Docker Desktop first." -ForegroundColor Red
    exit 1
}

# Check if docker-compose.yml exists
if (-not (Test-Path "docker-compose.yml")) {
    Write-Host "❌ docker-compose.yml not found. Please run this script from the project root." -ForegroundColor Red
    exit 1
}

Write-Host "📦 Building and starting containers..." -ForegroundColor Yellow

# Start the application
docker compose up -d --build

if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ Application started successfully!" -ForegroundColor Green
    Write-Host ""
    Write-Host "🌐 Access the application:" -ForegroundColor Cyan
    Write-Host "   Frontend: http://localhost:3000" -ForegroundColor White
    Write-Host "   Backend:  http://localhost:8001" -ForegroundColor White
    Write-Host "   API Docs: http://localhost:8001/docs" -ForegroundColor White
    Write-Host ""
    Write-Host "📊 Check status: docker compose ps" -ForegroundColor Yellow
    Write-Host "📋 View logs:   docker compose logs -f" -ForegroundColor Yellow
    Write-Host "🛑 Stop app:    docker compose down" -ForegroundColor Yellow
} else {
    Write-Host "❌ Failed to start application. Check logs with: docker compose logs" -ForegroundColor Red
}
