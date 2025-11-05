# 🚀 Phase 2B - Lokalne Uruchomienie na RTX 5070

## ✅ Status: GOTOWE DO URUCHOMIENIA

**Data**: 24 października 2025  
**GPU**: RTX 5070 (zweryfikowane)  
**Konfiguracje**: 3 scenariusze (miller_urey, hydrothermal, formamide)  
**Symulacje**: 30 total (10 per scenario, 500K steps each)

---

## ⚡ Szybki Start (5 minut)

### Krok 1: Test GPU (1 minuta)

```powershell
python test_gpu_quick.py
```

Oczekiwany output: `✅ ALL TESTS PASSED!`

### Krok 2: Test Symulacji (5 minut)

```powershell
python scripts/run_phase2_full.py `
  --config aws_test/configs/phase2_miller_urey_extended.yaml `
  --output results/test_local_miller_urey `
  --steps 10000 `
  --seed 42
```

### Krok 3: Jeśli Test Działa - Uruchom Pełny Batch

```powershell
# 10 symulacji Miller-Urey (10-20 godzin)
python run_phase2b_local.py --scenario miller_urey --runs 10
```

---

## 📊 Szczegóły Uruchomienia

### Pliki Utworzone

✅ `run_phase2b_local.py` - główny skrypt batch  
✅ `test_gpu_quick.py` - test GPU  
✅ `test_config.py` - test konfiguracji  
✅ `RUCHAM_LOCALNY_PHASE2B.md` - pełna dokumentacja  
✅ `QUICKSTART_LOCAL_PHASE2B.md` - szybki przewodnik  

### Konfiguracje

✅ `aws_test/configs/phase2_miller_urey_extended.yaml` - Miller-Urey  
✅ `aws_test/configs/phase2_hydrothermal_extended.yaml` - Hydrothermal  
✅ `aws_test/configs/phase2_formamide_extended.yaml` - Formamide  

Wszystkie z 500K kroków, zoptymalizowane dla Phase 2B.

---

## 🎯 Komendy Uruchomienia

### Testowe Uruchomienia

```powershell
# 1. Test GPU
python test_gpu_quick.py

# 2. Test konfiguracji
python test_config.py

# 3. Test symulacji (10K kroków, ~5 min)
python scripts/run_phase2_full.py --config aws_test/configs/phase2_miller_urey_extended.yaml --output results/test_local_miller_urey --steps 10000 --seed 42
```

### Pełne Uruchomienia

```powershell
# 1 symulacja - Miller-Urey (1-2h)
python scripts/run_phase2_full.py --config aws_test/configs/phase2_miller_urey_extended.yaml --output results/phase2b_local/miller_urey/run_01 --steps 500000 --seed 100

# 10 symulacji - Miller-Urey (10-20h)
python run_phase2b_local.py --scenario miller_urey --runs 10

# Wszystkie 30 symulacji (2-3 dni)
python run_phase2b_local.py --all --runs 10
```

---

## 📈 Oczekiwane Wyniki

### Po 30 Symulacjach:

| Metric | Obecnie | Po Phase 2B | Cel |
|--------|---------|-------------|-----|
| **Unikalne molekuły** | 6 | 80-140 | 100+ |
| **Cykle autokatalityczne** | 0 | 5-20 | 10+ |
| **Miller-Urey** | 0 | 30-50 | 30+ |
| **Hydrothermal** | 0 | 30-50 | 30+ |
| **Formamide** | 0 | 20-40 | 30+ |

### Dodatkowe:
- ✅ Reaction network
- ✅ Time-series danych  
- ✅ Molecular diversity metrics
- ✅ Autocatalytic cycles detection

---

## ⏱️ Timeline

| Etap | Symulacje | Czas | Opis |
|------|-----------|------|------|
| **Test** | 1 | 5 min | Test GPU + symulacja |
| **1 Pełna** | 1 | 1-2h | Pełna symulacja 500K steps |
| **Batch 10** | 10 | 10-20h | 10 symulacji jednego scenariusza |
| **Full 30** | 30 | 2-3 dni | Wszystkie 3 scenariusze |

---

## 💾 Wymagania Systemowe

### Hardwarowe:
- ✅ **GPU**: RTX 5070 (zweryfikowane działa)
- 💾 **RAM**: 32+ GB (zalecane)
- 💾 **Dysk**: 10-50 GB wolnego miejsca

### Softwarowe:
- ✅ Python 3.9+
- ✅ Taichi 1.7+ (CUDA)
- ✅ CUDA toolkit
- ✅ Wszystkie zależności z `requirements.txt`

---

## 🔍 Monitoring

### Sprawdź Postęp:

```powershell
# Lista ukończonych
Get-ChildItem -Recurse results/phase2b_local -Filter "results.json" | Select-Object DirectoryName

# Sprawdź logi
Get-Content results/phase2b_local/*/run_*/simulation.log -Tail 20
```

### Sprawdź Wyniki:

```powershell
# Podsumowanie
Get-ChildItem -Recurse results/phase2b_local -Filter "summary.txt"

# Molekuły
Get-ChildItem -Recurse results/phase2b_local -Filter "molecules.json"
```

---

## 📁 Struktura Wyników

```
results/phase2b_local/
├── miller_urey/
│   ├── run_01/
│   │   ├── results.json          # Główne wyniki
│   │   ├── molecules.json        # Molekuły
│   │   ├── summary.txt           # Podsumowanie
│   │   ├── simulation.log        # Log
│   │   └── snapshots/            # Snapshots
│   ├── run_02/
│   └── ... (10 runs total)
├── hydrothermal/
│   └── ... (10 runs)
└── formamide/
    └── ... (10 runs)
```

---

## 🎯 Następne Kroki Po Ukończeniu

### 1. Analiza Wyników

```powershell
python scripts/analyze_phase2_results.py results/phase2b_local
```

### 2. Generowanie Wykresów

```powershell
python scripts/generate_all_figures.py
```

### 3. Decyzja

- ✅ **100+ molekuł** → Przejdź do Phase 3 (pisanie papera)
- ⚠️ **<100 molekuł** → Rozważ więcej uruchomień
- 📊 Analiza cykli autokatalitycznych

---

## ⚠️ Rozwiązywanie Problemów

### GPU nie działa:

```powershell
python -c "import taichi as ti; ti.init(arch=ti.cuda); print('OK')"
```

Jeśli fail:
- Sprawdź sterowniki NVIDIA
- Sprawdź CUDA toolkit
- Sprawdź `nvidia-smi`

### Brak pamięci GPU:

Zmniejsz `n_particles` w konfiguracji (z 2000 do 1500)

### Spowolnienia:

Wyłącz `diagnostics.enabled` w konfiguracji

---

## 💡 Tips

### Równoległe Uruchamianie:

```powershell
# PowerShell background jobs
Start-Job { python run_phase2b_local.py --scenario miller_urey --runs 10 }
Start-Job { python run_phase2b_local.py --scenario hydrothermal --runs 10 }

# Sprawdź status
Get-Job
```

### Przerwanie i Kontynuacja:

Skrypty automatycznie przeskakują ukończone uruchomienia:

```powershell
# Restart (skoczy ukończone)
python run_phase2b_local.py --scenario miller_urey --runs 10
```

---

## 📚 Dokumentacja

- `RUCHAM_LOCALNY_PHASE2B.md` - pełna dokumentacja
- `QUICKSTART_LOCAL_PHASE2B.md` - szybki przewodnik
- `RUN_LOCAL_RTX5070.md` - oryginalna dokumentacja RTX 5070
- `PHASE2_NEXT_STEPS_PL.md` - następne kroki z analizy

---

## ✅ Status Checklist

- [x] GPU test - PASSED
- [x] Konfiguracje - READY
- [x] Skrypty - READY
- [x] Dokumentacja - READY
- [ ] Test symulacji (10K steps)
- [ ] 1 pełna symulacja (500K steps)
- [ ] Batch 10 symulacji
- [ ] Pełne 30 symulacji
- [ ] Analiza wyników

---

## 🚀 START HERE

```powershell
# Najpierw test (5 minut)
python scripts/run_phase2_full.py --config aws_test/configs/phase2_miller_urey_extended.yaml --output results/test_local_miller_urey --steps 10000 --seed 42
```

**Jeśli działa - uruchom pełny batch:**

```powershell
python run_phase2b_local.py --all --runs 10
```

**Szacowany czas**: 2-3 dni  
**Szacowany koszt**: $0 (lokalnie)  
**Oczekiwane wyniki**: 100+ molekuł, cykle autokatalityczne

---

**Powodzenia! 🎉**

