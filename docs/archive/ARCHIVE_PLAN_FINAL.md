---
date: 2025-11-23
label: plan
---

# PLAN ARCHIWIZACJI - Live 2.0 (FINAL)

**Data:** 2025-11-23  
**Podstawa:** `docs/STRUCTURE_DEEP_ANALYSIS.md`  
**Zasady:** Zgodnie z `docs/ARCHIVE_POLICY.md` i `docs/NAVIGATION_GUIDE.md`

---

## ⚠️ READ-ONLY ZONES (NIE DOTYKAMY)

Zgodnie z `ARCHIVE_POLICY.md` sekcja 3, **NIE PROponujemy** archiwizacji:
- `backend/sim/core/**`
- `backend/sim/chemistry/**`
- `backend/sim/config.py`
- `backend/sim/molecule_extractor.py`
- `scripts/run_phase2_full.py`
- `scripts/phase2.py`
- `aws_test/scripts/**`
- `aws_test/configs/phase2_*.yaml`
- `aws_test/configs/SUPER_FAST.yaml`
- `docs/phase2b/**`
- `docs/technical/**`
- `docs/VALIDATION_ROADMAP.md`
- `docs/SCIENTIFIC_VALIDITY_ANALYSIS.md`
- `docs/NAVIGATION_GUIDE_LIVE2.md`
- `results/**`
- `phase2b_results/**`
- `run/snapshots/**`
- `run_*/checkpoints/**`

---

## A. ARCHIVE/ONE_OFF_SCRIPTS → Skrypty jednorazowe, debug, testy chwilowe

### A.1 Skrypty diagnostyczne/debug (PEWNE - zgodnie z ARCHIVE_POLICY 4.1)

| Ścieżka | Nowa lokalizacja | Uzasadnienie |
|---------|------------------|--------------|
| `check_real_clusters.py` | `archive/one_off_scripts/check_real_clusters.py` | Jednorazowy skrypt do weryfikacji klastrów, użyty raz podczas debugowania |
| `diagnose_chemistry.py` | `archive/one_off_scripts/diagnose_chemistry.py` | Skrypt diagnostyczny do analizy chemii, jednorazowe użycie |
| `diagnose_round1.sh` | `archive/one_off_scripts/diagnose_round1.sh` | Diagnostyka konkretnej rundy (round1), jednorazowe |
| `fix_catalog_timeline.py` | `archive/one_off_scripts/fix_catalog_timeline.py` | Fix jednorazowy dla catalog timeline, prawdopodobnie już zastosowany |
| `force_cluster_detection.py` | `archive/one_off_scripts/force_cluster_detection.py` | Debug/test cluster detection, jednorazowe |
| `QUICK_RUN_PHASE2.py` | `archive/one_off_scripts/QUICK_RUN_PHASE2.py` | Quick test script, jednorazowe użycie |

### A.2 Skrypty AWS emergency/fix (PEWNE - zgodnie z ARCHIVE_POLICY 4.1)

| Ścieżka | Nowa lokalizacja | Uzasadnienie |
|---------|------------------|--------------|
| `aws_start_missing_9.sh` | `archive/one_off_scripts/aws_start_missing_9.sh` | Fix jednorazowy dla konkretnego problemu (missing run 9) |
| `CHECK_AWS_RESULTS.sh` | `archive/one_off_scripts/CHECK_AWS_RESULTS.sh` | Jednorazowy check wyników AWS, prawdopodobnie zastąpiony przez `aws_test/scripts/check_*` |
| `copy_fix_to_aws.ps1` | `archive/one_off_scripts/copy_fix_to_aws.ps1` | Jednorazowy fix do kopiowania na AWS |
| `copy_to_aws.ps1` | `archive/one_off_scripts/copy_to_aws.ps1` | Jednorazowy skrypt kopiowania, prawdopodobnie zastąpiony przez lepsze rozwiązanie |

### A.3 Skrypty benchmark/test (PEWNE - zgodnie z ARCHIVE_POLICY 4.1)

| Ścieżka | Nowa lokalizacja | Uzasadnienie |
|---------|------------------|--------------|
| `analyze_benchmark_results.ps1` | `archive/one_off_scripts/analyze_benchmark_results.ps1` | Analiza wyników benchmarków, jednorazowe użycie |
| `run_benchmark.ps1` | `archive/one_off_scripts/run_benchmark.ps1` | Skrypt benchmark, jednorazowe testy |
| `run_cpu_benchmark.ps1` | `archive/one_off_scripts/run_cpu_benchmark.ps1` | Benchmark CPU, jednorazowe testy |
| `run_hybrid_test.ps1` | `archive/one_off_scripts/run_hybrid_test.ps1` | Test hybrid mode, jednorazowe |

### A.4 Skrypty fix/cleanup (PEWNE - zgodnie z ARCHIVE_POLICY 4.1)

| Ścieżka | Nowa lokalizacja | Uzasadnienie |
|---------|------------------|--------------|
| `cleanup_processes.ps1` | `archive/one_off_scripts/cleanup_processes.ps1` | Cleanup procesów, jednorazowe użycie |
| `fix_taichi_version.ps1` | `archive/one_off_scripts/fix_taichi_version.ps1` | Fix wersji Taichi, jednorazowe (jeśli już zastosowany) |
| `fix_taichi_version.sh` | `archive/one_off_scripts/fix_taichi_version.sh` | Duplikat fix_taichi_version.ps1, wersja bash |

### A.5 Skrypty testowe - wersje iteracyjne (DO POTWIERDZENIA - zgodnie z ARCHIVE_POLICY 6.4)

| Ścieżka | Nowa lokalizacja | Uzasadnienie | Status |
|---------|------------------|--------------|--------|
| `backend/tests/check_status2.ps1` | `archive/one_off_scripts/backend_tests_check_status2.ps1` | Wersja iteracyjna testu, check_status.ps1 może być kanoniczny | ⚠️ DO POTWIERDZENIA |
| `backend/tests/check_status3.ps1` | `archive/one_off_scripts/backend_tests_check_status3.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `backend/tests/check_status4.ps1` | `archive/one_off_scripts/backend_tests_check_status4.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `backend/tests/check_status5.ps1` | `archive/one_off_scripts/backend_tests_check_status5.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `backend/tests/check_status6.ps1` | `archive/one_off_scripts/backend_tests_check_status6.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `backend/tests/check_status7.ps1` | `archive/one_off_scripts/backend_tests_check_status7.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `backend/tests/check_status8.ps1` | `archive/one_off_scripts/backend_tests_check_status8.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `tests/check_status2.ps1` | `archive/one_off_scripts/tests_check_status2.ps1` | Wersja iteracyjna testu (duplikat z backend/tests/) | ⚠️ DO POTWIERDZENIA |
| `tests/check_status3.ps1` | `archive/one_off_scripts/tests_check_status3.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `tests/check_status4.ps1` | `archive/one_off_scripts/tests_check_status4.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `tests/check_status5.ps1` | `archive/one_off_scripts/tests_check_status5.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `tests/check_status6.ps1` | `archive/one_off_scripts/tests_check_status6.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `tests/check_status7.ps1` | `archive/one_off_scripts/tests_check_status7.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |
| `tests/check_status8.ps1` | `archive/one_off_scripts/tests_check_status8.ps1` | Wersja iteracyjna testu | ⚠️ DO POTWIERDZENIA |

**Uwaga:** `check_status.ps1` w obu lokalizacjach może być kanoniczny - wymaga weryfikacji przed archiwizacją wersji 2-8.

### A.6 Skrypty start/stop/utility (DO POTWIERDZENIA)

| Ścieżka | Nowa lokalizacja | Uzasadnienie | Status |
|---------|------------------|--------------|--------|
| `start_backend_simple.ps1` | `archive/one_off_scripts/start_backend_simple.ps1` | Wersja "simple" start_backend, prawdopodobnie zastąpiona przez `start_backend.ps1` | ⚠️ DO POTWIERDZENIA |
| `kill_backend.ps1` | `archive/one_off_scripts/kill_backend.ps1` | Kill backend script - wymaga weryfikacji czy jest częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `create_new_sim.ps1` | `archive/one_off_scripts/create_new_sim.ps1` | Utility do tworzenia nowej symulacji - może być aktywnie używane | ⚠️ DO POTWIERDZENIA |
| `monitor_aws_runs.sh` | `archive/one_off_scripts/monitor_aws_runs.sh` | Monitoring AWS runs - może być aktywnie używane | ⚠️ DO POTWIERDZENIA |
| `setup_aws_instance.sh` | `archive/one_off_scripts/setup_aws_instance.sh` | Setup AWS instance - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `test_aws_instance.sh` | `archive/one_off_scripts/test_aws_instance.sh` | Test AWS instance - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `start.ps1` | `archive/one_off_scripts/start.ps1` | Główny skrypt start - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `start_backend.ps1` | `archive/one_off_scripts/start_backend.ps1` | Start backend - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `start_frontend.ps1` | `archive/one_off_scripts/start_frontend.ps1` | Start frontend - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `start_hydro_queue.ps1` | `archive/one_off_scripts/start_hydro_queue.ps1` | Start hydro queue - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `run_aws_production.sh` | `archive/one_off_scripts/run_aws_production.sh` | Run AWS production - może być aktywnie używane | ⚠️ DO POTWIERDZENIA |

### A.7 Skrypty Phase 2B (DO POTWIERDZENIA - mogą być aktywnie używane)

| Ścieżka | Nowa lokalizacja | Uzasadnienie | Status |
|---------|------------------|--------------|--------|
| `run_phase2b_hydro_queue.py` | `archive/one_off_scripts/run_phase2b_hydro_queue.py` | Runner Phase 2B dla hydrothermal - może być aktywnie używany lub zastąpiony przez wersję w `scripts/` | ⚠️ DO POTWIERDZENIA |
| `run_phase2b_local.py` | `archive/one_off_scripts/run_phase2b_local.py` | Runner Phase 2B lokalny - może być aktywnie używany lub zastąpiony przez wersję w `scripts/` | ⚠️ DO POTWIERDZENIA |

**Podsumowanie sekcji A:**
- ✅ **PEWNE do archiwum:** 18 plików (A.1-A.4)
- ⚠️ **DO POTWIERDZENIA:** 23 pliki (A.5-A.7)

---

## B. ARCHIVE/OLD_DOCS → Dokumenty zastąpione nowszymi wersjami

### B.1 Pliki .txt w root (DO POTWIERDZENIA - zgodnie z ARCHIVE_POLICY 4.2)

| Ścieżka | Nowa lokalizacja | Uzasadnienie | Status |
|---------|------------------|--------------|--------|
| `AWS_EMERGENCY_FIX.txt` | `archive/old_docs/AWS_EMERGENCY_FIX.txt` | Notatka emergency fix - wymaga weryfikacji czy informacje przeniesione do `docs/aws_test/` | ⚠️ DO POTWIERDZENIA |
| `aws_minimal_setup.txt` | `archive/old_docs/aws_minimal_setup.txt` | Notatka minimal setup - wymaga weryfikacji czy zastąpiona przez `docs/aws_test/` | ⚠️ DO POTWIERDZENIA |
| `AWS_RECOMMENDED_ACTION.txt` | `archive/old_docs/AWS_RECOMMENDED_ACTION.txt` | Notatka recommended action - wymaga weryfikacji czy zastąpiona | ⚠️ DO POTWIERDZENIA |
| `AWS_ROUND2_COMMANDS.txt` | `archive/old_docs/AWS_ROUND2_COMMANDS.txt` | Notatka commands round 2 - wymaga weryfikacji czy zastąpiona | ⚠️ DO POTWIERDZENIA |
| `HYDRO_SETUP_COMPLETE.txt` | `archive/old_docs/HYDRO_SETUP_COMPLETE.txt` | Notatka hydro setup - wymaga weryfikacji czy zastąpiona przez `docs/` | ⚠️ DO POTWIERDZENIA |
| `docs_structure.txt` | `archive/old_docs/docs_structure.txt` | Struktura dokumentacji - wymaga weryfikacji czy zastąpiony przez `docs/INDEX.md` lub `docs/NAVIGATION_GUIDE.md` | ⚠️ DO POTWIERDZENIA |

**Uwaga:** Wszystkie pliki .txt wymagają weryfikacji czy:
1. Informacje zostały przeniesione do właściwej dokumentacji w `docs/`
2. Pliki są już nieaktualne
3. Nie są używane jako szybkie referencje

**Podsumowanie sekcji B:**
- ⚠️ **DO POTWIERDZENIA:** 6 plików

---

## C. ARCHIVE/EXPERIMENTS → Prototypy, próbne wersje algorytmów, alternatywne pipeline'y

### C.1 Eksperymentalne wersje (DO POTWIERDZENIA - zgodnie z ARCHIVE_POLICY 4.3)

| Ścieżka | Nowa lokalizacja | Uzasadnienie | Status |
|---------|------------------|--------------|--------|
| `matcher/matcher.py` | `archive/experiments/matcher/matcher.py` | Stara wersja matchera - jeśli `matcher_v2.py` jest kanoniczny, to `matcher.py` to eksperyment | ⚠️ DO POTWIERDZENIA |
| `configs/` (root) | `archive/experiments/configs/` | 24 pliki YAML - wymaga weryfikacji czy to eksperymentalne konfiguracje testowe (testy, optymalizacje) czy aktywnie używane lokalnie. Jeśli `aws_test/configs/` są produkcyjne, to `configs/` może być eksperymentalne | ⚠️ DO POTWIERDZENIA |
| `diagnostics/` (root) | `archive/experiments/diagnostics/` | Katalog diagnostyczny - wymaga weryfikacji czy to eksperymentalna diagnostyka czy ma różne przeznaczenie niż `backend/diagnostics/` | ⚠️ DO POTWIERDZENIA |

**Uwaga:** Wszystkie elementy wymagają weryfikacji:
1. Czy są aktywnie używane
2. Czy mają unikalne funkcje
3. Czy są eksperymentami czy częścią produkcyjnego workflow

**Podsumowanie sekcji C:**
- ⚠️ **DO POTWIERDZENIA:** 3 elementy (1 plik + 2 katalogi)

---

## D. ARCHIVE/TMP_RESULTS → Wyniki nieużywane do publikacji / analizy

### D.1 Katalogi wyników AWS (DO POTWIERDZENIA - zgodnie z ARCHIVE_POLICY 4.4)

| Ścieżka | Nowa lokalizacja | Uzasadnienie | Status |
|---------|------------------|--------------|--------|
| `aws_results/` (root) | `archive/tmp_results/aws_results/` | Zawiera `miller_urey_extended/` z run_1 do run_18 - wymaga weryfikacji czy to duplikat `aws_test/results/` lub `results/phase2b_additional/`. Jeśli duplikat → archiwum. Jeśli unikalne dane → NIE RUSZAĆ | ⚠️ DO POTWIERDZENIA |
| `aws_test/results_16_completed/` | `archive/tmp_results/aws_test_results_16_completed/` | Archiwum 16 ukończonych runów - wymaga weryfikacji czy już nieużywane do analizy | ⚠️ DO POTWIERDZENIA |
| `aws_test/results_28_completed/` | `archive/tmp_results/aws_test_results_28_completed/` | Archiwum 28 ukończonych runów - wymaga weryfikacji czy już nieużywane do analizy | ⚠️ DO POTWIERDZENIA |
| `aws_test/results_all_completed/` | `archive/tmp_results/aws_test_results_all_completed/` | Archiwum wszystkich ukończonych - wymaga weryfikacji czy już nieużywane do analizy | ⚠️ DO POTWIERDZENIA |
| `phase2b_aws_results/` | `archive/tmp_results/phase2b_aws_results/` | Kolejne archiwum wyników Phase 2B - wymaga weryfikacji czy już nieużywane | ⚠️ DO POTWIERDZENIA |

**Uwaga:** 
- Jeśli te katalogi są duże (>2MB), rozważyć przeniesienie poza repo (np. lokalnie / cloud storage) i wstawić w repo jedynie link / README o lokalizacji (zgodnie z ARCHIVE_POLICY 4.4).
- Wszystkie katalogi wyników wymagają weryfikacji:
  1. Czy zawierają unikalne dane czy są duplikatami
  2. Czy są używane do analizy/publikacji
  3. Czy można je bezpiecznie zarchiwizować

**Podsumowanie sekcji D:**
- ⚠️ **DO POTWIERDZENIA:** 5 katalogów

---

## E. ARCHIVE/DEPRECATED → Kod/konfiguracje zastąpione finalną wersją

### E.1 Katalog `live2.0/` (DO POTWIERDZENIA - WYSOKI PRIORYTET - zgodnie z ARCHIVE_POLICY 6.4)

| Ścieżka | Nowa lokalizacja | Uzasadnienie | Status |
|---------|------------------|--------------|--------|
| `live2.0/` (root) | `archive/deprecated/live2.0/` | Zagnieżdżony katalog projektu - wygląda na błąd agenta. Zawiera `backend/` z pełną strukturą (api/, sim/, tests/, snapshots/). Wymaga **pilnej weryfikacji** czy to duplikat całego `backend/` czy ma unikalne pliki (szczególnie `live2.0/backend/snapshots/` może zawierać unikalne debug snapshots). Jeśli duplikat → archiwum. Jeśli ma unikalne pliki → wymaga decyzji | ⚠️ DO POTWIERDZENIA (WYSOKI PRIORYTET) |

**Uwaga:** 
- **KRYTYCZNE** - wymaga weryfikacji przed archiwizacją
- Sprawdzić czy `live2.0/backend/snapshots/` zawiera unikalne dane
- Jeśli duplikat → archiwum zgodnie z ARCHIVE_POLICY 6.4 (zduplikowany kod → wybierasz kanoniczny i archiwizujesz resztę)

### E.2 Zastąpione wersje skryptów (DO POTWIERDZENIA - zgodnie z ARCHIVE_POLICY 6.4)

| Ścieżka | Nowa lokalizacja | Uzasadnienie | Status |
|---------|------------------|--------------|--------|
| `start_backend_simple.ps1` | `archive/deprecated/start_backend_simple.ps1` | Wersja "simple" - jeśli `start_backend.ps1` jest kanoniczny, to `start_backend_simple.ps1` jest zastąpiony | ⚠️ DO POTWIERDZENIA |
| `fix_taichi_version.sh` | `archive/deprecated/fix_taichi_version.sh` | Jeśli `fix_taichi_version.ps1` jest kanoniczny (lub odwrotnie) - jedna wersja jest zastąpiona | ⚠️ DO POTWIERDZENIA |

**Podsumowanie sekcji E:**
- ⚠️ **DO POTWIERDZENIA:** 1 katalog (wysoki priorytet) + 2 pliki

---

## 📊 PODSUMOWANIE PLANU

### Statystyki:
- ✅ **PEWNE do archiwum:** 18 plików (sekcja A.1-A.4)
- ⚠️ **DO POTWIERDZENIA:** 39 elementów (pliki + katalogi)

### Rozkład kategorii:
- **A. one_off_scripts:** 18 pewnych + 23 do potwierdzenia = 41 plików
- **B. old_docs:** 6 plików do potwierdzenia
- **C. experiments:** 3 elementy do potwierdzenia (1 plik + 2 katalogi)
- **D. tmp_results:** 5 katalogów do potwierdzenia
- **E. deprecated:** 1 katalog (wysoki priorytet) + 2 pliki do potwierdzenia

### Priorytety weryfikacji:
1. **WYSOKI:** `live2.0/` (root) - może być duplikat całego backendu
2. **ŚREDNI:** `aws_results/` (root) - może być duplikat wyników
3. **ŚREDNI:** `configs/` (root) vs `aws_test/configs/` - weryfikacja przeznaczenia
4. **ŚREDNI:** `check_status*.ps1` - 14 plików iteracyjnych (wersje 2-8)
5. **NISKI:** Skrypty start/utility - weryfikacja czy część workflow
6. **NISKI:** Pliki .txt - weryfikacja czy zastąpione przez `docs/`

---

## 🔧 PROCES WYKONANIA (zgodnie z ARCHIVE_POLICY sekcja 5)

### Krok 1: ✅ Identyfikacja plików (WYKONANE)
Plan archiwizacji został wygenerowany na podstawie analizy struktury.

### Krok 2: ✅ Plan archiwizacji (WYKONANE)
Plan zawiera:
1. ✅ Dokładne ścieżki plików
2. ✅ Proponowane nowe lokalizacje
3. ✅ Uzasadnienie kategorii (one-off, old-docs, experiments, results, deprecated)

### Krok 3: ⏳ Czekanie na akceptację
**Agent nie rusza plików bez akceptacji!**

### Krok 4: ⏳ Po akceptacji - wykonanie operacji
Operacje będą wykonane przez `git mv`:
```bash
git mv check_real_clusters.py archive/one_off_scripts/
git mv diagnose_chemistry.py archive/one_off_scripts/
# ... itd dla wszystkich zatwierdzonych plików
```

### Krok 5: ⏳ Utworzenie ARCHIVE_LOG.md
Po wykonaniu operacji zostanie utworzony `archive/ARCHIVE_LOG.md` z wpisami w formacie:
```
[2025-11-23] Archived: check_real_clusters.py
Reason: One-off debugging script, used once during cluster debugging.
Moved to: archive/one_off_scripts/
```

---

## ⚠️ WAŻNE UWAGI

1. **Przed wykonaniem:** Wszystkie elementy oznaczone jako "DO POTWIERDZENIA" wymagają weryfikacji przez użytkownika
2. **Read-Only Zones:** Żaden element z read-only zones nie jest proponowany do archiwizacji
3. **Git mv:** Wszystkie operacje będą wykonane przez `git mv` aby zachować historię (zgodnie z ARCHIVE_POLICY sekcja 5)
4. **ARCHIVE_LOG.md:** Po wykonaniu planu zostanie utworzony log archiwizacji (zgodnie z ARCHIVE_POLICY sekcja 5, krok 4)
5. **Zasada kanoniczna:** Dla duplikatów wybieramy kanoniczną wersję i archiwizujemy resztę (zgodnie z ARCHIVE_POLICY 6.4)

---

## ❓ CZY ZATWIERDZASZ PLAN?

**Plan zawiera:**
- 18 plików **PEWNYCH** do archiwizacji (sekcja A.1-A.4)
- 39 elementów **DO POTWIERDZENIA** przed archiwizacją

**Proponowane działanie:**
1. **Zatwierdź 18 pewnych plików** → wykonanie `git mv` dla sekcji A.1-A.4
2. **Zweryfikuj elementy "DO POTWIERDZENIA"** → decyzja o każdym z osobna
3. **Po weryfikacji** → wykonanie `git mv` dla zatwierdzonych elementów
4. **Utworzenie ARCHIVE_LOG.md** → dokumentacja wszystkich operacji

**Czy zatwierdzasz plan?**

