$headers = @{'Content-Type' = 'application/json'}
$body = @{
    config = @{
        grid_height = 256
        grid_width = 256
        mode = 'open_chemistry'
        max_particles = 10000
        max_time = 1000
        dt = 0.01
        energy_decay = 0.95
        energy_threshold = 0.1
        particle_radius = 0.5
        binding_threshold = 0.8
        unbinding_threshold = 0.2
        novelty_window = 100
        min_cluster_size = 2
        vis_frequency = 10
        log_frequency = 100
    }
    mode = 'open_chemistry'
} | ConvertTo-Json -Depth 3

$response = Invoke-RestMethod -Uri 'http://localhost:8000/simulation/create' -Method POST -Headers $headers -Body $body
$id = $response.simulation_id
Write-Host "Created simulation: $id"

Invoke-RestMethod -Uri "http://localhost:8000/simulation/$id/start" -Method POST

1..3 | ForEach-Object {
    Start-Sleep -Milliseconds 500
    $status = Invoke-RestMethod -Uri "http://localhost:8000/simulation/$id/status"
    Write-Host "Time: $($status.current_time), Running: $($status.is_running), Steps: $($status.step_count)"
}
