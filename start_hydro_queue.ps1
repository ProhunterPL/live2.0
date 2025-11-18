# Quick Start - Hydrothermal Queue Runner
# ========================================
# Uruchamia symulacje hydrothermal od run_10 do run_1

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "HYDROTHERMAL QUEUE - QUICK START" -ForegroundColor Cyan
Write-Host "========================================`n" -ForegroundColor Cyan

# Check if virtual environment should be activated
if (Test-Path ".\activate_live_env.ps1") {
    Write-Host "Activating virtual environment..." -ForegroundColor Yellow
    & .\activate_live_env.ps1
}

# Check system info
Write-Host "`n[SYSTEM] System Info:" -ForegroundColor Yellow
$cpuCores = (Get-WmiObject Win32_ComputerSystem).NumberOfLogicalProcessors
Write-Host "   CPU Cores: $cpuCores" -ForegroundColor White

try {
    $gpuInfo = nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader
    Write-Host "   GPU: $gpuInfo" -ForegroundColor White
    Write-Host "`n[OK] GPU detected, but CPU is FASTER for this workload!" -ForegroundColor Green
    Write-Host "   (Tests showed CPU outperforms GPU for chemistry-heavy simulations)" -ForegroundColor Cyan
} catch {
    Write-Host "`n[OK] No GPU - will use CPU (optimal for this task!)" -ForegroundColor Green
}

# Menu
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Wybierz opcje:" -ForegroundColor White
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "1) [TEST] Test CPU (10K steps, ~5 minut)" -ForegroundColor White
Write-Host "2) [RUN]  Pelna kolejka CPU (run 10->1, ~10h) >>ZALECANE<<" -ForegroundColor Green
Write-Host "3) [HYBRID] Pelna kolejka HYBRID (run 10->1, ~7.5h) EKSPERYMENTALNE" -ForegroundColor Yellow
Write-Host "4) [STATUS] Sprawdz status" -ForegroundColor White
Write-Host "5) [MONITOR] Monitoruj zasoby (CPU/GPU)" -ForegroundColor White
Write-Host "6) [EXIT] Wyjdz" -ForegroundColor White
Write-Host "========================================`n" -ForegroundColor Cyan

$choice = Read-Host "Wybor (1-6)"

switch ($choice) {
    "1" {
        Write-Host "`n[TEST] Uruchamiam test CPU..." -ForegroundColor Green
        Write-Host "Uzywam CPU z $cpuCores rdzeniami`n" -ForegroundColor Cyan
        python scripts/run_phase2_full.py `
            --config aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml `
            --output results/test_hydro_local `
            --steps 10000 `
            --seed 42 `
            --force-cpu
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "`n[OK] Test zakonczony sukcesem!" -ForegroundColor Green
            Write-Host "`nMozesz teraz uruchomic pelna kolejke CPU:" -ForegroundColor Yellow
            Write-Host "python run_phase2b_hydro_queue.py --start 10 --end 1`n" -ForegroundColor Cyan
        } else {
            Write-Host "`n[ERROR] Test sie nie powiodl. Sprawdz logi." -ForegroundColor Red
        }
    }
    
    "2" {
        Write-Host "`n[RUN] Uruchamiam pelna kolejke CPU (10 runs)..." -ForegroundColor Green
        Write-Host "Mode: CPU z $cpuCores rdzeniami" -ForegroundColor Cyan
        Write-Host "Szacowany czas: ~10 godzin`n" -ForegroundColor Yellow
        
        $confirm = Read-Host "Kontynuowac? (y/n)"
        if ($confirm -eq "y" -or $confirm -eq "Y") {
            Write-Host "`nStarting queue... Mozesz przerwac przez Ctrl+C (bezpiecznie - restart bedzie kontynowal)`n" -ForegroundColor Yellow
            python run_phase2b_hydro_queue.py --start 10 --end 1
        } else {
            Write-Host "`nAnulowano." -ForegroundColor Yellow
        }
    }
    
    "3" {
        Write-Host "`n[HYBRID] Uruchamiam pelna kolejke HYBRID (10 runs)..." -ForegroundColor Yellow
        Write-Host "Mode: GPU (physics) + CPU (chemistry)" -ForegroundColor Cyan
        Write-Host "Szacowany czas: ~7.5 godzin" -ForegroundColor Yellow
        Write-Host "[WARNING] EKSPERYMENTALNE - moze wymagac dostosowania`n" -ForegroundColor Yellow
        
        $confirm = Read-Host "Kontynuowac? (y/n)"
        if ($confirm -eq "y" -or $confirm -eq "Y") {
            Write-Host "`nStarting hybrid queue...`n" -ForegroundColor Yellow
            python run_phase2b_hydro_queue.py --start 10 --end 1 --hybrid
        } else {
            Write-Host "`nAnulowano." -ForegroundColor Yellow
        }
    }
    
    "4" {
        Write-Host "`n[STATUS] Status ukonczonych runs:" -ForegroundColor Green
        
        $completed = Get-ChildItem -Path "results/phase2b_local/hydrothermal/*/results.json" -ErrorAction SilentlyContinue
        if ($completed) {
            Write-Host "`n[OK] Ukonczone runs:" -ForegroundColor Green
            foreach ($file in $completed) {
                $runDir = Split-Path -Parent $file
                $runName = Split-Path -Leaf $runDir
                Write-Host "  - $runName" -ForegroundColor White
            }
            Write-Host "`nTotal: $($completed.Count) / 10 runs`n" -ForegroundColor Cyan
        } else {
            Write-Host "`n[WARNING] Brak ukonczonych runs." -ForegroundColor Yellow
        }
        
        # Check if any logs exist
        $logs = Get-ChildItem -Path "logs/hydro_queue_*.log" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
        if ($logs) {
            Write-Host "Ostatni log: $($logs.Name)" -ForegroundColor Yellow
            Write-Host "`nOstatnie 10 linii logu:" -ForegroundColor Yellow
            Get-Content $logs.FullName -Tail 10
        }
    }
    
    "5" {
        Write-Host "`n[MONITOR] Monitorowanie zasobow..." -ForegroundColor Green
        Write-Host "`nCPU Usage (odswieza sie co 5 sekund, Ctrl+C zeby wyjsc):`n" -ForegroundColor Cyan
        
        Write-Host "Uwaga: Dla GPU monitoring uruchom w osobnym oknie: nvidia-smi -l 5`n" -ForegroundColor Yellow
        
        # Monitor CPU
        while ($true) {
            $cpu = Get-WmiObject Win32_Processor | Measure-Object -Property LoadPercentage -Average | Select-Object -ExpandProperty Average
            $mem = Get-WmiObject Win32_OperatingSystem
            $memUsed = [math]::Round(($mem.TotalVisibleMemorySize - $mem.FreePhysicalMemory) / 1MB, 2)
            $memTotal = [math]::Round($mem.TotalVisibleMemorySize / 1MB, 2)
            $memPercent = [math]::Round(($memUsed / $memTotal) * 100, 1)
            
            Write-Host "`r[$(Get-Date -Format 'HH:mm:ss')] CPU: $cpu% | RAM: $memUsed/$memTotal GB ($memPercent%)" -NoNewline
            Start-Sleep -Seconds 5
        }
    }
    
    "6" {
        Write-Host "`nDo widzenia!`n" -ForegroundColor Cyan
        exit
    }
    
    default {
        Write-Host "`n[ERROR] Nieprawidlowy wybor (1-6).`n" -ForegroundColor Red
    }
}

Write-Host "`n[TIP] CPU jest szybsze niz GPU dla tego typu symulacji!" -ForegroundColor Cyan

Write-Host "`n" 

