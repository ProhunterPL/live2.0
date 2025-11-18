# 🌊 Hydrothermal Queue - Lokalne Uruchomienie

## 🎯 Cel

Uruchomić symulacje hydrothermal lokalnie (run_10 → run_1) podczas gdy AWS kończy miller_urey.

## ⚡ Szybki Start

### 1. Test Pojedynczej Symulacji (5-10 minut)

Najpierw przetestuj jedną krótką symulację żeby sprawdzić że wszystko działa:

```powershell
python scripts/run_phase2_full.py `
  --config aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml `
  --output results/test_hydro_local `
  --steps 10000 `
  --seed 42
```

**Oczekiwany czas**: 5-10 minut  
**Sukces jeśli**: Pojawi się `results/test_hydro_local/results.json`

### 2. Uruchom Pełną Kolejkę (15 godzin)

Jeśli test działa, uruchom pełną kolejkę:

```powershell
python run_phase2b_hydro_queue.py --start 10 --end 1
```

**To uruchomi**: run_10, run_9, run_8, ..., run_2, run_1  
**Czas trwania**: ~15 godzin (90 min × 10 runs)  
**Output**: `results/phase2b_local/hydrothermal/run_XX/`

## 📊 Status Podczas Uruchamiania

### Sprawdź Co Się Dzieje

```powershell
# Sprawdź logi w czasie rzeczywistym
Get-Content logs/hydro_queue_*.log -Wait -Tail 20

# Lista ukończonych runs
Get-ChildItem results/phase2b_local/hydrothermal/*/results.json

# Sprawdź ile zostało
(Get-ChildItem results/phase2b_local/hydrothermal/ -Directory).Count
```

### Monitorowanie GPU

```powershell
# W osobnym oknie PowerShell
nvidia-smi -l 5  # Odświeża co 5 sekund
```

**Oczekiwane**: GPU Utilization: 80-100%

## 🎛️ Opcje Uruchomienia

### Wszystkie Runs (10 → 1)

```powershell
python run_phase2b_hydro_queue.py --start 10 --end 1
```

### Tylko Część (np. 10 → 5)

```powershell
python run_phase2b_hydro_queue.py --start 10 --end 5
```

### Dokończ Przerwane (np. od 7 → 1)

```powershell
python run_phase2b_hydro_queue.py --start 7 --end 1
```

**Skrypt automatycznie pominie ukończone runs!**

## ⏱️ Timeline

| Etap | Runs | Czas | Status |
|------|------|------|--------|
| **Test** | 1 (10K steps) | 5-10 min | Przed startem |
| **Run 10** | 1 | ~90 min | T+0h |
| **Run 9** | 1 | ~90 min | T+1.5h |
| **Run 8** | 1 | ~90 min | T+3h |
| **Run 7** | 1 | ~90 min | T+4.5h |
| **Run 6** | 1 | ~90 min | T+6h |
| **Run 5** | 1 | ~90 min | T+7.5h |
| **Run 4** | 1 | ~90 min | T+9h |
| **Run 3** | 1 | ~90 min | T+10.5h |
| **Run 2** | 1 | ~90 min | T+12h |
| **Run 1** | 1 | ~90 min | T+13.5h |
| **TOTAL** | 10 | **~15h** | **Done!** |

## 📁 Struktura Wyników

```
results/phase2b_local/hydrothermal/
├── run_10/
│   ├── results.json          # Główne wyniki
│   ├── molecules.json        # Wykryte molekuły
│   ├── simulation.log        # Log symulacji
│   ├── summary.txt           # Podsumowanie
│   └── snapshots/            # Snapshots co 50K steps
│       ├── step_00050000.json
│       ├── step_00100000.json
│       └── ... (10 snapshots total)
├── run_09/
├── run_08/
└── ... (10 runs total)
```

## 🔍 Sprawdź Wyniki Po Zakończeniu

### Podsumowanie Wszystkich Runs

```powershell
# Lista wszystkich ukończonych
Get-ChildItem results/phase2b_local/hydrothermal/*/results.json

# Zlicz molekuły we wszystkich runs
Get-ChildItem results/phase2b_local/hydrothermal/*/molecules.json | 
  ForEach-Object { (Get-Content $_ | ConvertFrom-Json).unique_molecules.Count }
```

### Analiza Wyników

```powershell
# Gdy wszystkie runs są gotowe
python scripts/analyze_phase2b_complete.py --scenario hydrothermal
```

## 🚨 Rozwiązywanie Problemów

### GPU Out of Memory

Jeśli GPU ma za mało pamięci:

1. Zmniejsz `n_particles` w config:
   ```yaml
   n_particles: 800  # Zamiast 1000
   ```

2. Lub użyj CPU (wolniejsze ale stabilniejsze):
   ```powershell
   python scripts/run_phase2_full.py `
     --config aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml `
     --output results/phase2b_local/hydrothermal/run_10 `
     --steps 500000 `
     --seed 109 `
     --force-cpu
   ```

### Symulacja "Zawiesiła Się"

Jeśli symulacja nie robi postępu (sprawdź w logu):

1. Ctrl+C żeby przerwać
2. Usuń częściowo ukończony run:
   ```powershell
   Remove-Item -Recurse results/phase2b_local/hydrothermal/run_XX
   ```
3. Uruchom ponownie - skrypt zacznie od tego runa

### Komputer Się Wyłącza / Restart

Skrypt jest **restart-safe**:

```powershell
# Po restarcie, po prostu uruchom ponownie tę samą komendę
python run_phase2b_hydro_queue.py --start 10 --end 1

# Skrypt automatycznie:
# - Znajdzie ukończone runs
# - Pominie je
# - Kontynuuje od pierwszego nieukończonego
```

## 💡 Tips & Tricks

### Uruchom w Tle (Overnight)

```powershell
# Metoda 1: Start-Job (PowerShell background)
Start-Job -ScriptBlock {
    cd C:\Users\user\Desktop\live2.0
    python run_phase2b_hydro_queue.py --start 10 --end 1
}

# Sprawdź status
Get-Job

# Metoda 2: nohup (jeśli masz Git Bash)
nohup python run_phase2b_hydro_queue.py --start 10 --end 1 > hydro_queue.log 2>&1 &
```

### Zapobiegaj Uśpieniu Komputera

```powershell
# PowerShell - zapobiega uśpieniu na 24h
powercfg /change standby-timeout-ac 1440
powercfg /change monitor-timeout-ac 0
```

**Przywróć normalnie po zakończeniu**:
```powershell
powercfg /change standby-timeout-ac 30
powercfg /change monitor-timeout-ac 10
```

### Równoległe Uruchamianie (Zaawansowane)

Jeśli masz 2 GPU lub dużo RAM, możesz uruchomić 2 symulacje równolegle:

```powershell
# Terminal 1
python run_phase2b_hydro_queue.py --start 10 --end 6

# Terminal 2
python run_phase2b_hydro_queue.py --start 5 --end 1
```

**Uwaga**: To wymaga ~16GB RAM + 2 GPU lub silny CPU

## 📊 Porównanie z AWS

| Aspekt | AWS (Miller-Urey) | Local (Hydrothermal) |
|--------|-------------------|---------------------|
| **Scenariusz** | miller_urey | hydrothermal |
| **Runs** | 11 done, 3 in progress | 10 planned |
| **Koszt** | ~$50-100 | $0 (lokalnie) |
| **Czas** | ~24h | ~15h |
| **Równoległość** | 4 runs jednocześnie | 1 run (lub 2 jeśli masz 2 GPU) |
| **Status** | Prawie skończone | Zaczynamy |

## 🎯 Następne Kroki Po Zakończeniu

1. **Sprawdź ile mamy unikalnych molekuł**:
   ```powershell
   python scripts/analyze_phase2b_complete.py --scenario hydrothermal
   ```

2. **Połącz z wynikami AWS**:
   - AWS: ~14 miller_urey runs
   - Local: 10 hydrothermal runs
   - **Total**: ~24 runs (cel: 30)

3. **Zdecyduj o formamide**:
   - Jeśli mamy już >100 unikalnych molekuł → możemy zacząć pisać paper
   - Jeśli <100 → uruchom formamide (10 runs, ~15h)

## ✅ Checklist

Przed rozpoczęciem:
- [ ] GPU test działa (`nvidia-smi`)
- [ ] Test 10K steps zakończony sukcesem
- [ ] Sprawdziłeś że masz ~50GB wolnego miejsca
- [ ] Ustawiłeś zapobieganie uśpieniu
- [ ] Gotowy na ~15h obliczeń

Podczas uruchamiania:
- [ ] Sprawdź GPU utilization (80-100%)
- [ ] Monitoruj logi co 1-2h
- [ ] Sprawdź postęp: ile runs ukończone

Po zakończeniu:
- [ ] Wszystkie 10 runs mają `results.json`
- [ ] Analiza wyników
- [ ] Połączenie z AWS results

---

## 🚀 START HERE

```powershell
# 1. Test (5 minut)
python scripts/run_phase2_full.py `
  --config aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml `
  --output results/test_hydro_local `
  --steps 10000 `
  --seed 42

# 2. Jeśli test OK, uruchom pełną kolejkę
python run_phase2b_hydro_queue.py --start 10 --end 1
```

**Szacowany czas**: 15 godzin  
**Szacowany koszt**: $0  
**Rezultat**: 10 hydrothermal runs gotowych do analizy

---

**Powodzenia! 🌊🔥**

