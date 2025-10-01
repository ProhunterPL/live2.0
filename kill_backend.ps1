# Quick kill script for Live 2.0 backend processes
# Use this when you need to quickly stop the backend

Write-Host "🛑 Stopping Live 2.0 backend..." -ForegroundColor Red

# Kill processes using port 8000
$connections = netstat -ano | Select-String ":8000 "
if ($connections) {
    $pids = $connections | ForEach-Object { 
        ($_ -split '\s+')[-1] | Where-Object { $_ -match '^\d+$' }
    } | Sort-Object -Unique
    
    foreach ($processId in $pids) {
        try {
            $process = Get-Process -Id $processId -ErrorAction SilentlyContinue
            if ($process -and $process.ProcessName -eq "python") {
                Write-Host "Killing Python process: PID $processId" -ForegroundColor Yellow
                Stop-Process -Id $processId -Force
            }
        }
        catch {
            # Ignore errors
        }
    }
}

# Kill any remaining Python processes that might be the backend
Get-Process -Name "python" -ErrorAction SilentlyContinue | ForEach-Object {
    try {
        Write-Host "Killing Python process: PID $($_.Id)" -ForegroundColor Yellow
        Stop-Process -Id $_.Id -Force
    }
    catch {
        # Ignore errors
    }
}

Write-Host "✅ Backend stopped!" -ForegroundColor Green
