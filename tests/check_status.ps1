Start-Sleep -Seconds 2
$status = Invoke-RestMethod -Uri 'http://localhost:8000/simulation/sim_1759004454738/status'
Write-Host "Time: $($status.current_time), Running: $($status.is_running), Steps: $($status.step_count)"
