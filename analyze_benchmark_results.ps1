# Analiza wyników benchmarku Hybrid Mode
# Sprawdza wyniki i daje rekomendacje

Write-Host "================================" -ForegroundColor Cyan
Write-Host "Analiza wyników benchmarku" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""

$resultsFile = "tests\hybrid_benchmark_results.json"

if (-not (Test-Path $resultsFile)) {
    Write-Host "❌ Plik z wynikami nie został znaleziony: $resultsFile" -ForegroundColor Red
    Write-Host ""
    Write-Host "Uruchom najpierw benchmark:" -ForegroundColor Yellow
    Write-Host "  .\run_hybrid_test.ps1" -ForegroundColor White
    Write-Host ""
    exit 1
}

Write-Host "📊 Wczytuję wyniki z: $resultsFile" -ForegroundColor Green
Write-Host ""

try {
    $results = Get-Content $resultsFile | ConvertFrom-Json
    
    # Filtruj tylko udane testy
    $successful = @($results | Where-Object { $_.success -eq $true })
    
    if ($successful.Count -eq 0) {
        Write-Host "❌ Brak udanych testów w wynikach" -ForegroundColor Red
        exit 1
    }
    
    # Sortuj wedlug wydajnosci symulacji
    $sortedBySim = $successful | Sort-Object -Property steps_per_sec -Descending
    $sortedByViz = $successful | Sort-Object -Property viz_avg_ms
    
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "WYNIKI WYDAJNOSCI" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host ""
    
    Write-Host "📈 Wydajnosc symulacji (steps/sec):" -ForegroundColor Yellow
    Write-Host ("-" * 70) -ForegroundColor Gray
    $slowest = $sortedBySim[-1].steps_per_sec
    foreach ($result in $sortedBySim) {
        $speedup = [math]::Round($result.steps_per_sec / $slowest, 1)
        $stepMs = [math]::Round($result.avg_step_ms, 1)
        $msg = "  {0,-40} {1,8:F1} steps/sec ({2,5:F1}ms/step) [{3}x]" -f $result.name, $result.steps_per_sec, $stepMs, $speedup
        Write-Host $msg -ForegroundColor White
    }
    
    Write-Host ""
    Write-Host "🎨 Wydajnosc wizualizacji (ms/frame, mniej = lepiej):" -ForegroundColor Yellow
    Write-Host ("-" * 70) -ForegroundColor Gray
    $slowestViz = $sortedByViz[-1].viz_avg_ms
    foreach ($result in $sortedByViz) {
        $speedup = [math]::Round($slowestViz / $result.viz_avg_ms, 1)
        $vizMs = [math]::Round($result.viz_avg_ms, 1)
        $msg = "  {0,-40} {1,6:F1}ms avg [{2}x szybsze]" -f $result.name, $vizMs, $speedup
        Write-Host $msg -ForegroundColor White
    }
    
    Write-Host ""
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "REKOMENDACJE" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host ""
    
    $bestSim = $sortedBySim[0]
    $bestViz = $sortedByViz[0]
    
    Write-Host "🚀 Najlepsze dla symulacji: $($bestSim.name)" -ForegroundColor Green
    $msg = "   Predkosc: {0:F1} steps/sec ({1:F1}ms/step)" -f $bestSim.steps_per_sec, $bestSim.avg_step_ms
    Write-Host $msg -ForegroundColor White
    Write-Host ""
    
    Write-Host "🎨 Najlepsze dla wizualizacji: $($bestViz.name)" -ForegroundColor Green
    $msg = "   Predkosc: {0:F1}ms na klatke" -f $bestViz.viz_avg_ms
    Write-Host $msg -ForegroundColor White
    Write-Host ""
    
    # Oblicz zwyciezce ogolnego
    if ($bestSim.name -eq $bestViz.name) {
        Write-Host "✅ ZWYCIĘZCA: $($bestSim.name)" -ForegroundColor Green
        Write-Host ""
        Write-Host "   Ten tryb jest najlepszy zarowno dla symulacji jak i wizualizacji!" -ForegroundColor White
    } else {
        # Oblicz zlozony wynik (70% symulacja, 30% wizualizacja)
        $bestSimScore = $bestSim.steps_per_sec
        $bestVizScore = $sortedByViz[0].viz_avg_ms
        
        $combinedScores = @()
        foreach ($result in $successful) {
            $simScore = $result.steps_per_sec / $bestSimScore
            $vizScore = $bestVizScore / $result.viz_avg_ms
            $combined = 0.7 * $simScore + 0.3 * $vizScore
            $combinedScores += [PSCustomObject]@{
                Name = $result.name
                Score = $combined
                Result = $result
            }
        }
        
        $bestOverall = ($combinedScores | Sort-Object -Property Score -Descending)[0]
        Write-Host "✅ ZWYCIĘZCA OGÓLNY: $($bestOverall.Name)" -ForegroundColor Green
        $msg = "   Zlozony wynik: {0:F2}" -f $bestOverall.Score
        Write-Host $msg -ForegroundColor White
        Write-Host ""
        Write-Host "   Uwaga: Rozne tryby sa najlepsze dla roznych zadan." -ForegroundColor Yellow
        Write-Host "   - Symulacja: $($bestSim.name)" -ForegroundColor White
        Write-Host "   - Wizualizacja: $($bestViz.name)" -ForegroundColor White
    }
    
    Write-Host ""
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "SZCZEGOLOWA ANALIZA" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host ""
    
    # Sprawdz czy Hybrid jest dostepny
    $hybridResult = $successful | Where-Object { $_.name -like "*Hybrid*" }
    if ($hybridResult) {
        Write-Host "🔬 Analiza trybu Hybrid:" -ForegroundColor Yellow
        Write-Host ""
        
        $gpuResult = $successful | Where-Object { $_.name -like "*Pure GPU*" -or $_.name -like "*CUDA*" }
        $cpuResult = $successful | Where-Object { $_.name -like "*Pure CPU*" }
        
        if ($gpuResult) {
            $hybridVsGpu = [math]::Round(($hybridResult.steps_per_sec / $gpuResult.steps_per_sec) * 100, 1)
            if ($hybridVsGpu -ge 95) {
                $msg = "   ✅ Hybrid jest prawie tak szybki jak Pure GPU ({0}%)" -f $hybridVsGpu
                Write-Host $msg -ForegroundColor Green
                Write-Host "      + Dodatkowo wykonuje analize chemiczna na CPU!" -ForegroundColor White
            } elseif ($hybridVsGpu -ge 80) {
                $hybridDiff = 100 - $hybridVsGpu
                $msg = "   ⚠️  Hybrid jest {0}% wolniejszy niz Pure GPU" -f $hybridDiff
                Write-Host $msg -ForegroundColor Yellow
                Write-Host "      Ale nadal wykonuje analize chemiczna rownolegle!" -ForegroundColor White
            } else {
                $msg = "   ❌ Hybrid jest znacznie wolniejszy ({0}% wydajnosci GPU)" -f $hybridVsGpu
                Write-Host $msg -ForegroundColor Red
                Write-Host "      Rozwaz uzycie Pure GPU jesli nie potrzebujesz analizy chemicznej" -ForegroundColor White
            }
        }
        
        if ($cpuResult) {
            $hybridVsCpu = [math]::Round(($hybridResult.steps_per_sec / $cpuResult.steps_per_sec), 1)
            Write-Host ""
            Write-Host "   📊 Hybrid vs Pure CPU: {0}x szybszy" -f $hybridVsCpu -ForegroundColor Cyan
        }
        
        Write-Host ""
        Write-Host "   💡 Zalecenia dla Hybrid:" -ForegroundColor Yellow
        Write-Host "      - Uzywaj gdy potrzebujesz pelnej analizy chemicznej" -ForegroundColor White
        Write-Host "      - Idealny dla symulacji z wiazaniami i klastrami" -ForegroundColor White
        Write-Host "      - Dobry kompromis miedzy wydajnoscia a funkcjonalnoscia" -ForegroundColor White
    } else {
        Write-Host "⚠️  Tryb Hybrid nie zostal przetestowany lub nie powiodl sie" -ForegroundColor Yellow
    }
    
    Write-Host ""
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "NASTEPNE KROKI" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host ""
    
    if ($bestSim.name -like "*Hybrid*" -or $bestViz.name -like "*Hybrid*") {
        Write-Host "1. ✅ Hybrid Mode jest rekomendowany!" -ForegroundColor Green
        Write-Host "   Uzyj HybridSimulationStepper w swoim kodzie" -ForegroundColor White
        Write-Host ""
        Write-Host "2. 📖 Zobacz dokumentacje:" -ForegroundColor Yellow
        Write-Host "   docs\HYBRID_GPU_CPU_GUIDE.md" -ForegroundColor White
        Write-Host ""
        Write-Host "3. 💻 Przykład użycia:" -ForegroundColor Yellow
        Write-Host "   from backend.sim.core.hybrid_stepper import HybridSimulationStepper" -ForegroundColor Gray
        Write-Host "   stepper = HybridSimulationStepper(config)" -ForegroundColor Gray
    } else {
        Write-Host "1. ⚠️  Hybrid Mode nie jest najszybszy dla tego przypadku" -ForegroundColor Yellow
        Write-Host "   Rozwaz uzycie: $($bestSim.name)" -ForegroundColor White
        Write-Host ""
        Write-Host "2. 💡 Hybrid moze byc nadal uzyteczny jesli potrzebujesz:" -ForegroundColor Yellow
        Write-Host "   - Pelnej analizy chemicznej - wiazania, klastry" -ForegroundColor White
        Write-Host "   - Wykrywania nowych substancji" -ForegroundColor White
        Write-Host "   - Zaawansowanej diagnostyki" -ForegroundColor White
    }
    
    Write-Host ""
    Write-Host "📊 Szczegolowe wyniki zapisane w: $resultsFile" -ForegroundColor Cyan
    Write-Host ""
    
} catch {
    Write-Host "❌ Blad podczas analizy wynikow: $_" -ForegroundColor Red
    Write-Host ""
    exit 1
}

