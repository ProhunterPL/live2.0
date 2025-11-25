---
date: 2025-11-23
label: analysis
---

# Głęboka Analiza Struktury Repozytorium Live 2.0

**Data analizy:** 2025-11-23  
**Źródło:** `structure.txt` (22429 linii)  
**Cel:** Identyfikacja duplikatów, chaosu strukturalnego, kandydatów do archiwizacji

---

## 🔍 1. DUPLIKATY PLIKÓW I KATALOGÓW

### 1.1 Duplikaty całych katalogów (KRYTYCZNE)

#### `live2.0/backend/` vs `backend/`
**Lokalizacja:** `live2.0/backend/` (root)  
**Zawartość:** Pełna struktura backendu (api/, sim/, tests/, snapshots/)  
**Obserwacja:** 
- Katalog `live2.0/` zawiera tylko `backend/`
- Struktura jest identyczna z głównym `backend/`
- Prawdopodobny błąd agenta przy tworzeniu zagnieżdżonej struktury
- **Status:** Prawdopodobny duplikat całego backendu

**Szczegóły duplikatu:**
- `live2.0/backend/sim/core/` - duplikat `backend/sim/core/`
- `live2.0/backend/tests/` - duplikat `backend/tests/`
- `live2.0/backend/api/` - duplikat `backend/api/`
- `live2.0/backend/snapshots/` - zawiera debug snapshots (może być unikalne)

### 1.2 Duplikaty skryptów testowych

#### `check_status*.ps1` - 8 wersji w dwóch lokalizacjach
**Lokalizacje:**
- `backend/tests/check_status.ps1` do `check_status8.ps1` (8 plików)
- `tests/check_status.ps1` do `check_status8.ps1` (8 plików)

**Obserwacja:**
- Identyczne nazwy w dwóch katalogach
- 16 plików łącznie (8+8)
- Prawdopodobnie iteracyjne wersje testowe
- **Status:** Duplikaty lub wersje rozwojowe

### 1.3 Duplikaty skryptów Phase 2B

#### `run_phase2b_*.py` w różnych lokalizacjach
**Lokalizacje:**
- Root: `run_phase2b_hydro_queue.py`, `run_phase2b_local.py`
- `aws_test/`: `run_phase2b_aws.sh`, `run_phase2b_master.py`, `run_phase2b_additional.py`

**Obserwacja:**
- Różne implementacje Phase 2B w różnych miejscach
- Może być zamierzone (różne środowiska), ale wymaga weryfikacji
- **Status:** Potencjalne duplikaty lub różne wersje

### 1.4 Duplikaty matcher

#### `matcher.py` vs `matcher_v2.py`
**Lokalizacja:** `matcher/`  
**Obserwacja:**
- Dwie wersje matchera w tym samym katalogu
- `matcher_v2.py` prawdopodobnie nowsza wersja
- **Status:** Stara wersja może być do archiwizacji

---

## 🗂️ 2. CHAOS STRUKTURALNY

### 2.1 Katalog `live2.0/` w root (KRYTYCZNY CHAOS)

**Problem:** Zagnieżdżony katalog projektu w samym projekcie  
**Struktura:**
```
live2.0/
└── backend/
    ├── api/
    ├── sim/
    ├── tests/
    └── snapshots/
```

**Obserwacja:**
- Wygląda na błąd agenta przy tworzeniu struktury
- Zawiera duplikat `backend/`
- Może zawierać unikalne snapshots debug (wymaga weryfikacji)
- **Wpływ:** Wprowadza konfuzję, zajmuje miejsce, może powodować błędy importów

### 2.2 Wiele lokalizacji wyników (CHAOS ORGANIZACYJNY)

**Obserwowane lokalizacje:**
1. `aws_results/` (root) - zawiera `miller_urey_extended/` z run_1 do run_18
2. `results/` (root) - zawiera pliki `.txt` i `.md`
3. `aws_test/results/` - wyniki z AWS
4. `aws_test/results_16_completed/` - archiwum
5. `aws_test/results_28_completed/` - archiwum
6. `aws_test/results_all_completed/` - archiwum
7. `results/phase2b_additional/` - główne wyniki Phase 2B
8. `phase2b_aws_results/` - kolejne archiwum

**Obserwacja:**
- 8 różnych lokalizacji na wyniki
- Brak jasnej hierarchii
- Trudno określić, które są aktywne, które archiwalne
- **Wpływ:** Konfuzja, ryzyko duplikacji danych, trudność w zarządzaniu

### 2.3 Skrypty w root zamiast w `scripts/` (CHAOS ORGANIZACYJNY)

**Obserwacja:** ~30+ skryptów w root, które powinny być w `scripts/`:
- Diagnostyczne: `check_real_clusters.py`, `diagnose_chemistry.py`
- AWS: `aws_start_missing_9.sh`, `CHECK_AWS_RESULTS.sh`
- Benchmark: `run_benchmark.ps1`, `run_cpu_benchmark.ps1`
- Start/Stop: `start_backend.ps1`, `start_backend_simple.ps1`, `kill_backend.ps1`
- Fix: `fix_taichi_version.ps1`, `fix_taichi_version.sh`

**Wpływ:**
- Utrudnia nawigację
- Root jest zaśmiecony
- Trudno znaleźć właściwy skrypt

### 2.4 Pliki .txt w root (CHAOS DOKUMENTACYJNY)

**Obserwowane pliki:**
- `AWS_EMERGENCY_FIX.txt`
- `aws_minimal_setup.txt`
- `AWS_RECOMMENDED_ACTION.txt`
- `AWS_ROUND2_COMMANDS.txt`
- `HYDRO_SETUP_COMPLETE.txt`
- `docs_structure.txt`

**Obserwacja:**
- Notatki tekstowe zamiast dokumentacji w `docs/`
- Brak struktury dokumentacyjnej
- **Wpływ:** Trudno znaleźć dokumentację, notatki mogą być przestarzałe

### 2.5 Duplikaty konfiguracji

**Obserwowane lokalizacje:**
- `configs/` (root) - 24 pliki YAML (testy, optymalizacje)
- `aws_test/configs/` - 14 plików YAML (produkcyjne, SUPER_FAST)

**Obserwacja:**
- Dwie lokalizacje konfiguracji
- Może być zamierzone (testowe vs produkcyjne), ale wymaga weryfikacji
- **Wpływ:** Możliwa konfuzja, które konfiguracje są aktualne

### 2.6 Duplikaty diagnostyki

**Obserwowane lokalizacje:**
- `diagnostics/` (root)
- `backend/diagnostics/`

**Obserwacja:**
- Dwa katalogi diagnostyczne
- Wymaga weryfikacji czy mają różne przeznaczenie
- **Wpływ:** Możliwa konfuzja

---

## 📦 3. KANDYDACI DO ARCHIWIZACJI

### 3.1 KANDYDACI DO `archive/one_off_scripts/`

#### Skrypty diagnostyczne/debug (PEWNE)
- `check_real_clusters.py` - jednorazowy check
- `diagnose_chemistry.py` - diagnostyka
- `diagnose_round1.sh` - diagnostyka konkretnej rundy
- `fix_catalog_timeline.py` - fix jednorazowy
- `force_cluster_detection.py` - debug/test
- `QUICK_RUN_PHASE2.py` - quick test

#### Skrypty AWS emergency (PEWNE)
- `aws_start_missing_9.sh` - fix jednorazowy
- `CHECK_AWS_RESULTS.sh` - jednorazowy check
- `copy_fix_to_aws.ps1` - jednorazowy fix
- `copy_to_aws.ps1` - jednorazowy copy

#### Skrypty benchmark/test (PEWNE)
- `analyze_benchmark_results.ps1`
- `run_benchmark.ps1`
- `run_cpu_benchmark.ps1`
- `run_hybrid_test.ps1`

#### Skrypty fix (PEWNE)
- `cleanup_processes.ps1`
- `fix_taichi_version.ps1` i `fix_taichi_version.sh` (jeśli już zastosowane)

#### Skrypty testowe - wersje iteracyjne (DO WERYFIKACJI)
- `backend/tests/check_status2.ps1` do `check_status8.ps1` (7 plików)
- `tests/check_status2.ps1` do `check_status8.ps1` (7 plików)
- **Uwaga:** `check_status.ps1` może być kanoniczny, wersje 2-8 to iteracje

### 3.2 KANDYDACI DO `archive/old_docs/`

#### Pliki .txt w root (DO WERYFIKACJI)
- `AWS_EMERGENCY_FIX.txt` - jeśli informacje przeniesione do `docs/`
- `aws_minimal_setup.txt` - jeśli zastąpiony przez `docs/aws_test/`
- `AWS_RECOMMENDED_ACTION.txt` - jeśli zastąpiony
- `AWS_ROUND2_COMMANDS.txt` - jeśli zastąpiony
- `HYDRO_SETUP_COMPLETE.txt` - jeśli zastąpiony
- `docs_structure.txt` - jeśli zastąpiony przez `docs/INDEX.md` lub `docs/NAVIGATION_GUIDE.md`

### 3.3 KANDYDACI DO `archive/experiments/`

#### Eksperymentalne wersje (DO WERYFIKACJI)
- `matcher/matcher.py` - jeśli `matcher_v2.py` jest kanoniczny
- `configs/` (root) - jeśli to eksperymentalne konfiguracje testowe
- `diagnostics/` (root) - jeśli to eksperymentalna diagnostyka

#### Alternatywne implementacje Phase 2B (DO WERYFIKACJI)
- `run_phase2b_hydro_queue.py` (root) - jeśli zastąpiony przez wersję w `scripts/`
- `run_phase2b_local.py` (root) - jeśli zastąpiony przez wersję w `scripts/`

### 3.4 KANDYDACI DO `archive/deprecated/`

#### Duplikaty całych katalogów (WYSOKI PRIORYTET)
- `live2.0/` (root) - **KRYTYCZNY** - prawdopodobny duplikat `backend/`
  - **Uwaga:** Wymaga weryfikacji czy `live2.0/backend/snapshots/` zawiera unikalne dane

#### Zastąpione wersje skryptów (DO WERYFIKACJI)
- `start_backend_simple.ps1` - jeśli `start_backend.ps1` jest kanoniczny
- `fix_taichi_version.sh` - jeśli `fix_taichi_version.ps1` jest kanoniczny (lub odwrotnie)

### 3.5 KANDYDACI DO `archive/tmp_results/`

#### Archiwalne wyniki (DO WERYFIKACJI)
- `aws_results/` (root) - jeśli duplikat `aws_test/results/` lub `results/phase2b_additional/`
- `aws_test/results_16_completed/` - jeśli już nieużywane
- `aws_test/results_28_completed/` - jeśli już nieużywane
- `aws_test/results_all_completed/` - jeśli już nieużywane
- `phase2b_aws_results/` - jeśli już nieużywane

**Uwaga:** Jeśli te katalogi są duże (>2MB), rozważyć przeniesienie poza repo z README wskazującym lokalizację.

---

## 🎯 4. PRIORYTETY WERYFIKACJI

### WYSOKI PRIORYTET (KRYTYCZNE)
1. **`live2.0/backend/`** - duplikat całego backendu
   - Weryfikacja: Czy zawiera unikalne pliki (szczególnie snapshots)?
   - Działanie: Jeśli duplikat → `archive/deprecated/`

### ŚREDNI PRIORYTET (WAŻNE)
2. **`aws_results/` (root)** - duplikat wyników?
   - Weryfikacja: Czy zawiera unikalne dane czy duplikat?
   - Działanie: Jeśli duplikat → `archive/tmp_results/`

3. **`configs/` (root) vs `aws_test/configs/`**
   - Weryfikacja: Czy `configs/` to eksperymenty czy aktywnie używane?
   - Działanie: Jeśli eksperymenty → `archive/experiments/`

4. **`check_status*.ps1`** - 16 plików (8+8)
   - Weryfikacja: Która wersja jest kanoniczna?
   - Działanie: Wersje 2-8 → `archive/one_off_scripts/`

### NISKI PRIORYTET (ORGANIZACYJNE)
5. **Skrypty w root** - ~30+ plików
   - Weryfikacja: Które są aktywnie używane?
   - Działanie: Jednorazowe → `archive/one_off_scripts/`

6. **Pliki .txt w root** - 6 plików
   - Weryfikacja: Czy informacje przeniesione do `docs/`?
   - Działanie: Jeśli zastąpione → `archive/old_docs/`

---

## 📊 STATYSTYKI CHAOSU

### Duplikaty:
- **Katalogi:** 1 (live2.0/backend/)
- **Skrypty:** ~16 (check_status*.ps1 w dwóch lokalizacjach)
- **Wersje:** 2 (matcher.py vs matcher_v2.py)

### Chaos strukturalny:
- **Lokalizacje wyników:** 8 różnych
- **Skrypty w root:** ~30+
- **Pliki .txt w root:** 6
- **Lokalizacje konfiguracji:** 2
- **Lokalizacje diagnostyki:** 2

### Kandydaci do archiwizacji:
- **one_off_scripts:** ~20+ plików (pewnych) + ~14 (do weryfikacji)
- **old_docs:** 6 plików (do weryfikacji)
- **experiments:** 3-5 elementów (do weryfikacji)
- **deprecated:** 1 katalog (wysoki priorytet) + 2 pliki (do weryfikacji)
- **tmp_results:** 5 katalogów (do weryfikacji)

---

## ⚠️ UWAGI I OSTRZEŻENIA

### Read-Only Zones (NIE DOTYKAMY)
Zgodnie z zasadami, **NIE PROponujemy** archiwizacji:
- `backend/sim/core/**`
- `backend/sim/chemistry/**`
- `scripts/run_phase2_full.py`
- `aws_test/configs/**`
- `docs/phase2b/**`
- `docs/technical/**`
- `results/phase2b_additional/**` (ukończone runy)

### Wymagana weryfikacja przed archiwizacją
1. **`live2.0/backend/snapshots/`** - może zawierać unikalne debug snapshots
2. **`aws_results/`** - może zawierać unikalne dane
3. **`configs/` (root)** - może być aktywnie używane lokalnie
4. **Skrypty start/utility** - mogą być częścią workflow
5. **Pliki .txt** - mogą być używane jako szybkie referencje

---

## 🎯 REKOMENDACJE

### Natychmiastowe działania (po weryfikacji):
1. Zweryfikować `live2.0/backend/` - jeśli duplikat → archiwum
2. Zweryfikować `aws_results/` - jeśli duplikat → archiwum
3. Przenieść skrypty jednorazowe z root do `archive/one_off_scripts/`

### Długoterminowe:
1. Ujednolicić lokalizacje wyników
2. Przenieść dokumentację .txt do `docs/`
3. Uporządkować skrypty (root → `scripts/` lub archiwum)
4. Wybrać kanoniczne wersje duplikatów

---

*Analiza wykonana na podstawie `structure.txt` z 2025-11-23*

