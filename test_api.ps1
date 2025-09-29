# Test script to check API functionality

Write-Host "🧪 Testing Live 2.0 API..." -ForegroundColor Yellow

$baseUrl = "http://localhost:8000"

# Test 1: Health check
Write-Host "`n1. Testing health endpoint..." -ForegroundColor Cyan
try {
    $response = Invoke-WebRequest -Uri "$baseUrl/" -Method GET
    $health = $response.Content | ConvertFrom-Json
    Write-Host "✅ Health check OK: $($health.message) v$($health.version)" -ForegroundColor Green
} catch {
    Write-Host "❌ Health check failed: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}

# Test 2: Create simulation
Write-Host "`n2. Testing simulation creation..." -ForegroundColor Cyan
$config = @{
    config = @{
        grid_height = 128
        grid_width = 128
        mode = "open_chemistry"
        max_particles = 5000
        max_time = 1000
        dt = 0.01
        energy_decay = 0.95
        energy_threshold = 0.1
        particle_radius = 0.5
        binding_threshold = 0.8
        unbinding_threshold = 0.2
        novelty_window = 100
        min_cluster_size = 2
        vis_frequency = 5
        log_frequency = 100
        seed = 123456
    }
    mode = "open_chemistry"
}

try {
    $body = $config | ConvertTo-Json -Depth 10
    $response = Invoke-WebRequest -Uri "$baseUrl/simulation/create" -Method POST -Body $body -ContentType "application/json"
    $result = $response.Content | ConvertFrom-Json
    
    if ($result.success) {
        Write-Host "✅ Simulation created: $($result.simulation_id)" -ForegroundColor Green
        $simId = $result.simulation_id
    } else {
        Write-Host "❌ Failed to create simulation: $($result.message)" -ForegroundColor Red
        exit 1
    }
} catch {
    Write-Host "❌ Create simulation failed: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}

# Test 3: Get simulation status
Write-Host "`n3. Testing simulation status..." -ForegroundColor Cyan
try {
    $response = Invoke-WebRequest -Uri "$baseUrl/simulation/$simId/status" -Method GET
    $status = $response.Content | ConvertFrom-Json
    Write-Host "✅ Status: Running=$($status.is_running), Paused=$($status.is_paused), Particles=$($status.particle_count)" -ForegroundColor Green
} catch {
    Write-Host "❌ Get status failed: $($_.Exception.Message)" -ForegroundColor Red
}

# Test 4: Start simulation
Write-Host "`n4. Testing simulation start..." -ForegroundColor Cyan
try {
    $response = Invoke-WebRequest -Uri "$baseUrl/simulation/$simId/start" -Method POST
    Write-Host "✅ Simulation started successfully" -ForegroundColor Green
} catch {
    Write-Host "❌ Start simulation failed: $($_.Exception.Message)" -ForegroundColor Red
}

# Test 5: Check status after start
Write-Host "`n5. Testing status after start..." -ForegroundColor Cyan
Start-Sleep -Seconds 2
try {
    $response = Invoke-WebRequest -Uri "$baseUrl/simulation/$simId/status" -Method GET
    $status = $response.Content | ConvertFrom-Json
    Write-Host "✅ Status: Running=$($status.is_running), Paused=$($status.is_paused), Steps=$($status.step_count)" -ForegroundColor Green
} catch {
    Write-Host "❌ Get status after start failed: $($_.Exception.Message)" -ForegroundColor Red
}

Write-Host "`n🎉 API tests completed!" -ForegroundColor Green
Write-Host "Simulation ID: $simId" -ForegroundColor White
