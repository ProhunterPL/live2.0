# Watch and Auto-Analyze
# ======================
# Monitors simulation progress and automatically runs analysis when complete
#
# Usage:
#   .\scripts\watch_and_analyze.ps1 -ResultDir "results/overnight_test_2025-10-13_18-17-09"

param(
    [Parameter(Mandatory=$true)]
    [string]$ResultDir,
    
    [int]$CheckInterval = 60,  # Check every 60 seconds
    
    [switch]$FullAnalysis = $false  # Include PubChem matching
)

$ErrorActionPreference = "Continue"

Write-Host ""
Write-Host "=" * 70 -ForegroundColor Cyan
Write-Host "SIMULATION MONITOR & AUTO-ANALYZER" -ForegroundColor Cyan
Write-Host "=" * 70 -ForegroundColor Cyan
Write-Host ""
Write-Host "Monitoring: $ResultDir" -ForegroundColor Yellow
Write-Host "Check interval: $CheckInterval seconds" -ForegroundColor Yellow
Write-Host "Full analysis: $FullAnalysis" -ForegroundColor Yellow
Write-Host ""

# Paths
$logFile = Join-Path $ResultDir "stderr.log"
$doneFile = Join-Path $ResultDir ".analysis_complete"

if (-not (Test-Path $ResultDir)) {
    Write-Host "ERROR: Result directory not found: $ResultDir" -ForegroundColor Red
    exit 1
}

# Function to check simulation status
function Get-SimulationStatus {
    if (-not (Test-Path $logFile)) {
        return @{ Status = "Unknown"; Step = 0; Total = 0; Progress = 0 }
    }
    
    # Get last 50 lines
    $recentLines = Get-Content $logFile -Tail 50
    
    # Look for step completion messages
    $stepLines = $recentLines | Where-Object { $_ -match "Step (\d+) completed" }
    if ($stepLines) {
        $lastStep = $stepLines[-1]
        if ($lastStep -match "Step (\d+) completed") {
            $currentStep = [int]$Matches[1]
        }
    } else {
        $currentStep = 0
    }
    
    # Look for total steps (from start of log)
    $totalSteps = 100000  # Default
    $headerLines = Get-Content $logFile -Head 100
    foreach ($line in $headerLines) {
        if ($line -match "max_steps[:\s]+(\d+)") {
            $totalSteps = [int]$Matches[1]
            break
        }
    }
    
    # Check if completed
    $isComplete = $false
    $completionLines = $recentLines | Where-Object { $_ -match "(Simulation complete|SIMULATION COMPLETE|Final snapshot)" }
    if ($completionLines -or $currentStep -ge $totalSteps) {
        $isComplete = $true
    }
    
    $progress = if ($totalSteps -gt 0) { ($currentStep / $totalSteps) * 100 } else { 0 }
    
    return @{
        Status = if ($isComplete) { "Complete" } else { "Running" }
        Step = $currentStep
        Total = $totalSteps
        Progress = $progress
        IsComplete = $isComplete
    }
}

# Function to estimate time remaining
function Get-TimeRemaining {
    param($status, $startTime)
    
    if ($status.Progress -le 0) {
        return "Unknown"
    }
    
    $elapsed = (Get-Date) - $startTime
    $totalTime = $elapsed.TotalSeconds / ($status.Progress / 100.0)
    $remaining = $totalTime - $elapsed.TotalSeconds
    
    if ($remaining -lt 0) {
        return "Finishing..."
    }
    
    $hours = [Math]::Floor($remaining / 3600)
    $minutes = [Math]::Floor(($remaining % 3600) / 60)
    
    return "{0:D2}:{1:D2}" -f $hours, $minutes
}

# Check if analysis already done
if (Test-Path $doneFile) {
    Write-Host ""
    Write-Host "Analysis already completed for this run." -ForegroundColor Green
    Write-Host "Delete $doneFile to re-run analysis." -ForegroundColor Yellow
    Write-Host ""
    exit 0
}

# Main monitoring loop
$startTime = Get-Date
$lastStep = 0
$stuckCount = 0

Write-Host "Monitoring started at: $(Get-Date -Format 'HH:mm:ss')" -ForegroundColor Cyan
Write-Host "Press Ctrl+C to stop monitoring" -ForegroundColor Gray
Write-Host ""

while ($true) {
    $status = Get-SimulationStatus
    
    # Check if stuck
    if ($status.Step -eq $lastStep -and $status.Status -eq "Running") {
        $stuckCount++
    } else {
        $stuckCount = 0
    }
    
    $lastStep = $status.Step
    
    # Display status
    $timeRemaining = Get-TimeRemaining -status $status -startTime $startTime
    $progressBar = "[" + ("=" * [Math]::Floor($status.Progress / 2)) + (" " * (50 - [Math]::Floor($status.Progress / 2))) + "]"
    
    Write-Host "`r$progressBar $($status.Progress.ToString('F1'))% - Step $($status.Step)/$($status.Total) - ETA: $timeRemaining" -NoNewline -ForegroundColor Cyan
    
    # Check if stuck for too long
    if ($stuckCount -gt 10) {
        Write-Host ""
        Write-Host ""
        Write-Host "WARNING: Simulation appears stuck (no progress for $($stuckCount * $CheckInterval) seconds)" -ForegroundColor Red
        Write-Host "Last step: $($status.Step)" -ForegroundColor Yellow
        Write-Host ""
        Write-Host "Options:" -ForegroundColor Yellow
        Write-Host "  1. Check if process is still running" -ForegroundColor Gray
        Write-Host "  2. Review error log: $logFile" -ForegroundColor Gray
        Write-Host "  3. Press Ctrl+C to stop monitoring" -ForegroundColor Gray
        Write-Host ""
    }
    
    # Check if complete
    if ($status.IsComplete) {
        Write-Host ""
        Write-Host ""
        Write-Host "=" * 70 -ForegroundColor Green
        Write-Host "SIMULATION COMPLETE!" -ForegroundColor Green
        Write-Host "=" * 70 -ForegroundColor Green
        Write-Host ""
        Write-Host "Final step: $($status.Step)" -ForegroundColor Cyan
        Write-Host "Total time: $(((Get-Date) - $startTime).ToString('hh\:mm\:ss'))" -ForegroundColor Cyan
        Write-Host ""
        
        # Run analysis
        Write-Host "Starting automatic analysis..." -ForegroundColor Yellow
        Write-Host ""
        
        # Quick analysis
        Write-Host "[1/4] Quick analysis..." -ForegroundColor Cyan
        if ($FullAnalysis) {
            python scripts/quick_analyze.py $ResultDir --full
        } else {
            python scripts/quick_analyze.py $ResultDir
        }
        
        # Network analysis
        Write-Host ""
        Write-Host "[2/4] Building reaction network..." -ForegroundColor Cyan
        python scripts/reaction_network_analyzer.py $ResultDir --export both
        
        # Check if network was created
        $networkFile = "analysis/reaction_network/reaction_network.json"
        if (Test-Path $networkFile) {
            # Autocatalytic cycles
            Write-Host ""
            Write-Host "[3/4] Detecting autocatalytic cycles..." -ForegroundColor Cyan
            python scripts/autocatalytic_detector.py $networkFile
            
            # Visualization
            Write-Host ""
            Write-Host "[4/4] Creating visualizations..." -ForegroundColor Cyan
            $cyclesFile = "analysis/autocatalytic_cycles/autocatalytic_cycles.json"
            if (Test-Path $cyclesFile) {
                python scripts/network_visualizer.py $networkFile --cycles $cyclesFile --interactive
            } else {
                python scripts/network_visualizer.py $networkFile --interactive
            }
        } else {
            Write-Host "WARNING: Network analysis failed, skipping cycles and visualization" -ForegroundColor Yellow
        }
        
        # Mark as done
        "Analysis completed at $(Get-Date)" | Out-File $doneFile
        
        Write-Host ""
        Write-Host "=" * 70 -ForegroundColor Green
        Write-Host "ANALYSIS COMPLETE!" -ForegroundColor Green
        Write-Host "=" * 70 -ForegroundColor Green
        Write-Host ""
        Write-Host "Results:" -ForegroundColor Cyan
        Write-Host "  - Quick summary: analysis/summary.txt" -ForegroundColor Gray
        Write-Host "  - Network data: analysis/reaction_network/" -ForegroundColor Gray
        Write-Host "  - Cycles: analysis/autocatalytic_cycles/" -ForegroundColor Gray
        Write-Host "  - Visualizations: analysis/visualizations/" -ForegroundColor Gray
        Write-Host ""
        
        break
    }
    
    # Wait before next check
    Start-Sleep -Seconds $CheckInterval
}

Write-Host "Monitoring stopped." -ForegroundColor Cyan
Write-Host ""

