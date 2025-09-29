# Cleanup script for Live 2.0 - kills hanging Python/uvicorn processes
# Run this if you get "port 8000 is already in use" errors

Write-Host "🧹 Cleaning up Live 2.0 processes..." -ForegroundColor Yellow

# Function to kill processes by port
function Kill-ProcessByPort {
    param([int]$Port)
    
    try {
        $connections = netstat -ano | Select-String ":$Port "
        if ($connections) {
            Write-Host "Found processes using port $Port:" -ForegroundColor Cyan
            $connections | ForEach-Object { Write-Host "  $_" -ForegroundColor Gray }
            
            $pids = $connections | ForEach-Object { 
                ($_ -split '\s+')[-1] | Where-Object { $_ -match '^\d+$' }
            } | Sort-Object -Unique
            
            foreach ($pid in $pids) {
                try {
                    $process = Get-Process -Id $pid -ErrorAction SilentlyContinue
                    if ($process) {
                        Write-Host "Killing process: $($process.ProcessName) (PID: $pid)" -ForegroundColor Red
                        Stop-Process -Id $pid -Force
                    }
                }
                catch {
                    Write-Host "Could not kill process $pid: $($_.Exception.Message)" -ForegroundColor Yellow
                }
            }
        } else {
            Write-Host "No processes found using port $Port" -ForegroundColor Green
        }
    }
    catch {
        Write-Host "Error checking port $Port: $($_.Exception.Message)" -ForegroundColor Red
    }
}

# Function to kill Python processes by name
function Kill-PythonProcesses {
    try {
        $pythonProcesses = Get-Process -Name "python" -ErrorAction SilentlyContinue
        if ($pythonProcesses) {
            Write-Host "Found Python processes:" -ForegroundColor Cyan
            $pythonProcesses | ForEach-Object { 
                Write-Host "  $($_.ProcessName) (PID: $($_.Id))" -ForegroundColor Gray 
            }
            
            $pythonProcesses | ForEach-Object {
                try {
                    Write-Host "Killing Python process: PID $($_.Id)" -ForegroundColor Red
                    Stop-Process -Id $_.Id -Force
                }
                catch {
                    Write-Host "Could not kill Python process $($_.Id): $($_.Exception.Message)" -ForegroundColor Yellow
                }
            }
        } else {
            Write-Host "No Python processes found" -ForegroundColor Green
        }
    }
    catch {
        Write-Host "Error checking Python processes: $($_.Exception.Message)" -ForegroundColor Red
    }
}

# Function to kill uvicorn processes
function Kill-UvicornProcesses {
    try {
        $uvicornProcesses = Get-Process | Where-Object { $_.ProcessName -eq "python" -and $_.CommandLine -like "*uvicorn*" } -ErrorAction SilentlyContinue
        if ($uvicornProcesses) {
            Write-Host "Found uvicorn processes:" -ForegroundColor Cyan
            $uvicornProcesses | ForEach-Object { 
                Write-Host "  $($_.ProcessName) (PID: $($_.Id))" -ForegroundColor Gray 
            }
            
            $uvicornProcesses | ForEach-Object {
                try {
                    Write-Host "Killing uvicorn process: PID $($_.Id)" -ForegroundColor Red
                    Stop-Process -Id $_.Id -Force
                }
                catch {
                    Write-Host "Could not kill uvicorn process $($_.Id): $($_.Exception.Message)" -ForegroundColor Yellow
                }
            }
        } else {
            Write-Host "No uvicorn processes found" -ForegroundColor Green
        }
    }
    catch {
        Write-Host "Error checking uvicorn processes: $($_.Exception.Message)" -ForegroundColor Red
    }
}

# Main cleanup
Write-Host "`n🔍 Checking for processes using port 8000..." -ForegroundColor Yellow
Kill-ProcessByPort -Port 8000

Write-Host "`n🔍 Checking for Python processes..." -ForegroundColor Yellow
Kill-PythonProcesses

Write-Host "`n🔍 Checking for uvicorn processes..." -ForegroundColor Yellow
Kill-UvicornProcesses

# Wait a moment for processes to terminate
Start-Sleep -Seconds 2

# Verify cleanup
Write-Host "`n✅ Verifying cleanup..." -ForegroundColor Yellow
$remainingConnections = netstat -ano | Select-String ":8000 "
if ($remainingConnections) {
    Write-Host "⚠️  Warning: Some processes may still be using port 8000:" -ForegroundColor Yellow
    $remainingConnections | ForEach-Object { Write-Host "  $_" -ForegroundColor Gray }
} else {
    Write-Host "✅ Port 8000 is now free!" -ForegroundColor Green
}

Write-Host "`n🎉 Cleanup complete! You can now start the backend." -ForegroundColor Green
Write-Host "Run: .\start_backend.ps1" -ForegroundColor Cyan
