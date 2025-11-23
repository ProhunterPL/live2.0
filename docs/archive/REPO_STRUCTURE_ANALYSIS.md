---
date: 2025-11-23
label: analysis
---

# Analiza Struktury Repozytorium Live 2.0

**Data analizy:** 2025-11-23  
**Źródło:** `structure.txt` (22429 linii)  
**Cel:** Identyfikacja duplikatów, jednorazowych skryptów, starych dokumentów, chaosu od agenta, katalogów do archiwizacji

---

## 📋 1. KATALOGI, KTÓRE WYGLĄDAJĄ NA DUPLIKATY

### 1.1 Potencjalne duplikaty wyników
- **`aws_results/`** (root) - zawiera `miller_urey_extended/` i `run_1/` do `run_18/`
- **`results/`** (root) - zawiera pliki `.txt` i `.md`
- **`aws_test/results/`** - zawiera wyniki z AWS
- **`aws_test/results_16_completed/`** - archiwum ukończonych
- **`aws_test/results_28_completed/`** - archiwum ukończonych
- **`aws_test/results_all_completed/`** - archiwum ukończonych
- **`phase2b_aws_results/`** (w structure.txt) - kolejne archiwum wyników

**Obserwacja:** Wiele miejsc na wyniki - `aws_results/` w root wygląda na duplikat `aws_test/results/` lub `results/phase2b_additional/`

### 1.2 Duplikaty diagnostyki
- **`diagnostics/`** (root) - katalog diagnostyczny
- **`backend/diagnostics/`** - diagnostyka w backendzie

**Obserwacja:** Dwa katalogi diagnostyczne - wymaga weryfikacji czy mają różne przeznaczenie

### 1.3 Duplikaty konfiguracji
- **`configs/`** (root) - 24 pliki YAML (testy, optymalizacje, scenariusze)
- **`aws_test/configs/`** - 14 plików YAML (produkcyjne, SUPER_FAST)

**Obserwacja:** `configs/` w root wygląda na testowe/eksperymentalne, `aws_test/configs/` na produkcyjne - może być OK, ale wymaga weryfikacji

### 1.4 Katalog `live2.0/` (root)
- **`live2.0/`** - zawiera `backend/` i inne podkatalogi

**Obserwacja:** Wygląda na **błąd agenta** - prawdopodobnie próba utworzenia zagnieżdżonej struktury projektu. To jest **prawdopodobny duplikat całego projektu**.

---

## 🔧 2. SKRYPTY JEDNORAZOWE / DEBUGUJĄCE (w root)

### 2.1 Skrypty diagnostyczne/debug
- `check_real_clusters.py` - jednorazowy check
- `diagnose_chemistry.py` - diagnostyka
- `diagnose_round1.sh` - diagnostyka konkretnej rundy
- `fix_catalog_timeline.py` - fix jednorazowy
- `force_cluster_detection.py` - debug/test
- `QUICK_RUN_PHASE2.py` - quick test

### 2.2 Skrypty AWS (jednorazowe/emergency)
- `AWS_EMERGENCY_FIX.txt` - notatka emergency
- `aws_minimal_setup.txt` - notatka setup
- `AWS_RECOMMENDED_ACTION.txt` - notatka
- `AWS_ROUND2_COMMANDS.txt` - notatka
- `aws_start_missing_9.sh` - fix jednorazowy
- `CHECK_AWS_RESULTS.sh` - jednorazowy check
- `copy_fix_to_aws.ps1` - jednorazowy fix
- `copy_to_aws.ps1` - jednorazowy copy
- `setup_aws_instance.sh` - setup (może być OK)
- `test_aws_instance.sh` - test (może być OK)

### 2.3 Skrypty start/stop (może być OK, ale wiele wersji)
- `start.ps1`
- `start_backend.ps1`
- `start_backend_simple.ps1` - wersja "simple"
- `start_frontend.ps1`
- `start_hydro_queue.ps1`
- `kill_backend.ps1`

**Obserwacja:** Wiele wersji start scripts - `start_backend_simple.ps1` wygląda na jednorazową wersję

### 2.4 Skrypty benchmark/test
- `analyze_benchmark_results.ps1` - analiza benchmarków
- `run_benchmark.ps1` - benchmark
- `run_cpu_benchmark.ps1` - benchmark CPU
- `run_hybrid_test.ps1` - test hybrid

### 2.5 Skrypty Phase 2B (może być OK)
- `run_phase2b_hydro_queue.py` - może być OK
- `run_phase2b_local.py` - może być OK

### 2.6 Skrypty fix/cleanup
- `cleanup_processes.ps1` - cleanup
- `fix_taichi_version.ps1` - fix jednorazowy
- `fix_taichi_version.sh` - fix jednorazowy (duplikat .ps1?)
- `create_new_sim.ps1` - utility (może być OK)
- `monitor_aws_runs.sh` - monitoring (może być OK)

**Obserwacja:** Większość tych skryptów wygląda na jednorazowe/debug i powinna być w `scripts/` lub `archive/one_off_scripts/`

---

## 📄 3. STARE DOKUMENTY .MD POZA `docs/`

### 3.1 W root
- **`README.md`** - OK, standardowy plik projektu
- **`docs_structure.txt`** - struktura dokumentacji (może być do archiwizacji lub przeniesienia do `docs/`)

### 3.2 W innych katalogach (do weryfikacji)
- `paper/` - zawiera dokumenty .md (OK, to jest katalog paper)
- `scripts/` - może zawierać .md (do weryfikacji)
- `aws_test/scripts/` - zawiera 2 pliki .md (do weryfikacji)

**Obserwacja:** `docs_structure.txt` w root wygląda na tymczasowy plik pomocniczy - może być przeniesiony do `docs/` lub archiwum

---

## 🤖 4. MIEJSCA, GDZIE AGENT PRAWDOPODOBNIE NATWORZYŁ CHAOS

### 4.1 Katalog `live2.0/` w root
**Problem:** Zagnieżdżony katalog projektu w samym projekcie  
**Prawdopodobna przyczyna:** Błąd agenta przy tworzeniu struktury  
**Obserwacja:** Wygląda na duplikat całego projektu - wymaga pilnej weryfikacji

### 4.2 Wiele plików .txt w root (notatki AWS)
- `AWS_EMERGENCY_FIX.txt`
- `aws_minimal_setup.txt`
- `AWS_RECOMMENDED_ACTION.txt`
- `AWS_ROUND2_COMMANDS.txt`
- `HYDRO_SETUP_COMPLETE.txt`

**Problem:** Notatki tekstowe zamiast dokumentacji w `docs/`  
**Obserwacja:** Te pliki powinny być w `docs/aws_test/` lub `archive/old_docs/`

### 4.3 Duplikaty skryptów start
- `start_backend.ps1` vs `start_backend_simple.ps1`
- `fix_taichi_version.ps1` vs `fix_taichi_version.sh`

**Problem:** Wiele wersji tego samego skryptu  
**Obserwacja:** Wymaga wyboru kanonicznej wersji

### 4.4 Skrypty w root zamiast w `scripts/`
**Problem:** ~30+ skryptów w root, które powinny być w `scripts/`  
**Obserwacja:** Agent prawdopodobnie tworzył skrypty ad-hoc w root zamiast w odpowiednim katalogu

### 4.5 Wiele katalogów wyników
**Problem:** `aws_results/`, `results/`, `aws_test/results/`, `aws_test/results_*_completed/`, `phase2b_aws_results/`  
**Obserwacja:** Agent prawdopodobnie tworzył nowe katalogi zamiast używać istniejących

---

## 📦 5. KATALOGI, KTÓRE NADAJĄ SIĘ DO PRZENIESIENIA DO `archive/`

### 5.1 Do `archive/one_off_scripts/`
**Wszystkie skrypty z sekcji 2** (diagnostyczne, debug, jednorazowe fixes):
- `check_real_clusters.py`
- `diagnose_chemistry.py`
- `diagnose_round1.sh`
- `fix_catalog_timeline.py`
- `force_cluster_detection.py`
- `QUICK_RUN_PHASE2.py`
- `aws_start_missing_9.sh`
- `CHECK_AWS_RESULTS.sh`
- `copy_fix_to_aws.ps1`
- `copy_to_aws.ps1`
- `start_backend_simple.ps1` (jeśli `start_backend.ps1` jest kanoniczny)
- `analyze_benchmark_results.ps1`
- `run_benchmark.ps1`
- `run_cpu_benchmark.ps1`
- `run_hybrid_test.ps1`
- `cleanup_processes.ps1`
- `fix_taichi_version.ps1` i `fix_taichi_version.sh` (jeśli już nieużywane)
- `kill_backend.ps1` (jeśli nie jest częścią standardowego workflow)

### 5.2 Do `archive/old_docs/`
- `docs_structure.txt` (jeśli zastąpiony przez `docs/INDEX.md` lub `docs/NAVIGATION_GUIDE.md`)
- `AWS_EMERGENCY_FIX.txt` (jeśli informacje przeniesione do `docs/`)
- `aws_minimal_setup.txt` (jeśli zastąpiony przez `docs/aws_test/`)
- `AWS_RECOMMENDED_ACTION.txt` (jeśli zastąpiony)
- `AWS_ROUND2_COMMANDS.txt` (jeśli zastąpiony)
- `HYDRO_SETUP_COMPLETE.txt` (jeśli zastąpiony)

### 5.3 Do `archive/experiments/`
- `configs/` (root) - jeśli to eksperymentalne konfiguracje, a `aws_test/configs/` są produkcyjne
- `diagnostics/` (root) - jeśli to eksperymentalna diagnostyka, a `backend/diagnostics/` jest kanoniczna

### 5.4 Do `archive/deprecated/`
- `live2.0/` (root) - jeśli to błąd i duplikat projektu
- `start_backend_simple.ps1` - jeśli `start_backend.ps1` jest nowszą wersją
- `fix_taichi_version.sh` - jeśli `fix_taichi_version.ps1` jest nowszą wersją (lub odwrotnie)

### 5.5 Do `archive/tmp_results/` (lub poza repo jeśli duże)
- `aws_results/` (root) - jeśli to duplikat `aws_test/results/` lub `results/phase2b_additional/`
- `aws_test/results_16_completed/` - jeśli już nieużywane
- `aws_test/results_28_completed/` - jeśli już nieużywane
- `aws_test/results_all_completed/` - jeśli już nieużywane
- `phase2b_aws_results/` - jeśli już nieużywane

---

## ⚠️ UWAGI I OSTRZEŻENIA

### 6.1 Read-Only Zones (NIE RUSZAĆ)
Zgodnie z `docs/NAVIGATION_GUIDE.md` i `docs/ARCHIVE_POLICY.md`, **NIE WOLNO** archiwizować:
- `backend/sim/core/**`
- `backend/sim/chemistry/**`
- `scripts/run_phase2_full.py`
- `aws_test/configs/phase2_*.yaml`
- `aws_test/configs/*SUPER_FAST*.yaml`
- `docs/phase2b/**`
- `docs/technical/**`
- `results/phase2b_additional/**` (ukończone runy)

### 6.2 Wymagana weryfikacja przed archiwizacją
Przed przeniesieniem do `archive/` należy zweryfikować:
1. Czy `live2.0/` to rzeczywiście duplikat czy ma jakieś unikalne pliki
2. Czy `aws_results/` to duplikat czy zawiera unikalne dane
3. Czy `configs/` (root) to eksperymenty czy aktywnie używane konfiguracje
4. Czy `diagnostics/` (root) i `backend/diagnostics/` mają różne przeznaczenie
5. Czy skrypty w root są rzeczywiście jednorazowe czy część workflow

### 6.3 Potencjalne problemy
- **Duplikacja wyników:** Wiele miejsc na wyniki może prowadzić do konfuzji
- **Skrypty w root:** Utrudnia nawigację, powinny być w `scripts/`
- **Notatki .txt:** Powinny być w `docs/` jako właściwa dokumentacja
- **Katalog `live2.0/`:** Jeśli to duplikat, zajmuje miejsce i wprowadza chaos

---

## 📊 PODSUMOWANIE STATYSTYK

- **Skrypty jednorazowe w root:** ~30+
- **Pliki .txt (notatki) w root:** 5
- **Potencjalne duplikaty katalogów:** 4-5
- **Katalogi wyników:** 6+ różnych lokalizacji
- **Katalogi do archiwizacji:** ~10-15

---

**Uwaga:** Ta analiza jest **tylko obserwacją**. Przed wykonaniem jakichkolwiek zmian należy:
1. Zweryfikować każdy element
2. Upewnić się, że nie archiwizujemy read-only zones
3. Stworzyć plan archiwizacji zgodnie z `docs/ARCHIVE_POLICY.md`
4. Uzyskać akceptację przed wykonaniem `git mv`

---

*Analiza wykonana na podstawie `structure.txt` z 2025-11-23*

