Start-Sleep -Seconds 3
$status = Invoke-RestMethod -Uri 'http://localhost:8000/simulation/sim_1759007240781/status'
Write-Host "Backend Status - Time: $($status.current_time), Running: $($status.is_running), Steps: $($status.step_count)"
Write-Host "Frontend should be available at: http://localhost:3000"
Write-Host "Backend API available at: http://localhost:8000"
