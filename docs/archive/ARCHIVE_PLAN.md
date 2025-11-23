---
date: 2025-11-23
label: plan
---

# PLAN ARCHIWIZACJI - Live 2.0

**Data:** 2025-11-23  
**Podstawa:** `docs/REPO_STRUCTURE_ANALYSIS.md`  
**Zasady:** Zgodnie z `docs/ARCHIVE_POLICY.md` i `docs/NAVIGATION_GUIDE.md`

---

## ⚠️ READ-ONLY ZONES (NIE DOTYKAMY)

Zgodnie z zasadami, **NIE PROponujemy** archiwizacji:
- `backend/sim/core/**`
- `backend/sim/chemistry/**`
- `scripts/run_phase2_full.py`
- `aws_test/configs/**` (wszystkie konfiguracje AWS)
- `docs/phase2b/**`
- `docs/technical/**`
- `results/**` (wszystkie wyniki w głównym katalogu results/)
- `phase2b_*_results/**`

---

## A. ARCHIVE/ONE_OFF_SCRIPTS → Skrypty jednorazowe

### A.1 Skrypty diagnostyczne/debug (PEWNE)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `check_real_clusters.py` | Jednorazowy skrypt do weryfikacji klastrów, prawdopodobnie użyty raz podczas debugowania | ✅ DO ARCHIWUM |
| `diagnose_chemistry.py` | Skrypt diagnostyczny do analizy chemii, jednorazowe użycie | ✅ DO ARCHIWUM |
| `diagnose_round1.sh` | Diagnostyka konkretnej rundy (round1), jednorazowe | ✅ DO ARCHIWUM |
| `fix_catalog_timeline.py` | Fix jednorazowy dla catalog timeline, prawdopodobnie już zastosowany | ✅ DO ARCHIWUM |
| `force_cluster_detection.py` | Debug/test cluster detection, jednorazowe | ✅ DO ARCHIWUM |
| `QUICK_RUN_PHASE2.py` | Quick test script, jednorazowe użycie | ✅ DO ARCHIWUM |

### A.2 Skrypty AWS emergency/fix (PEWNE)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `aws_start_missing_9.sh` | Fix jednorazowy dla konkretnego problemu (missing run 9) | ✅ DO ARCHIWUM |
| `CHECK_AWS_RESULTS.sh` | Jednorazowy check wyników AWS, prawdopodobnie zastąpiony przez `aws_test/scripts/check_*` | ✅ DO ARCHIWUM |
| `copy_fix_to_aws.ps1` | Jednorazowy fix do kopiowania na AWS | ✅ DO ARCHIWUM |
| `copy_to_aws.ps1` | Jednorazowy skrypt kopiowania, prawdopodobnie zastąpiony przez lepsze rozwiązanie | ✅ DO ARCHIWUM |

### A.3 Skrypty benchmark/test (PEWNE)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `analyze_benchmark_results.ps1` | Analiza wyników benchmarków, jednorazowe użycie | ✅ DO ARCHIWUM |
| `run_benchmark.ps1` | Skrypt benchmark, jednorazowe testy | ✅ DO ARCHIWUM |
| `run_cpu_benchmark.ps1` | Benchmark CPU, jednorazowe testy | ✅ DO ARCHIWUM |
| `run_hybrid_test.ps1` | Test hybrid mode, jednorazowe | ✅ DO ARCHIWUM |

### A.4 Skrypty fix/cleanup (PEWNE)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `cleanup_processes.ps1` | Cleanup procesów, jednorazowe użycie | ✅ DO ARCHIWUM |
| `fix_taichi_version.ps1` | Fix wersji Taichi, jednorazowe (jeśli już zastosowany) | ✅ DO ARCHIWUM |
| `fix_taichi_version.sh` | Duplikat fix_taichi_version.ps1, wersja bash | ✅ DO ARCHIWUM |

### A.5 Skrypty start/stop (DO POTWIERDZENIA)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `start_backend_simple.ps1` | Wersja "simple" start_backend, prawdopodobnie zastąpiona przez `start_backend.ps1` | ⚠️ DO POTWIERDZENIA |
| `kill_backend.ps1` | Kill backend script - wymaga weryfikacji czy jest częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |

### A.6 Skrypty utility (DO POTWIERDZENIA)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `create_new_sim.ps1` | Utility do tworzenia nowej symulacji - może być aktywnie używane | ⚠️ DO POTWIERDZENIA |
| `monitor_aws_runs.sh` | Monitoring AWS runs - może być aktywnie używane | ⚠️ DO POTWIERDZENIA |
| `setup_aws_instance.sh` | Setup AWS instance - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `test_aws_instance.sh` | Test AWS instance - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |

### A.7 Skrypty Phase 2B (DO POTWIERDZENIA)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `run_phase2b_hydro_queue.py` | Runner Phase 2B dla hydrothermal - może być aktywnie używany | ⚠️ DO POTWIERDZENIA |
| `run_phase2b_local.py` | Runner Phase 2B lokalny - może być aktywnie używany | ⚠️ DO POTWIERDZENIA |

### A.8 Skrypty start (DO POTWIERDZENIA - mogą być częścią workflow)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `start.ps1` | Główny skrypt start - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `start_backend.ps1` | Start backend - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `start_frontend.ps1` | Start frontend - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |
| `start_hydro_queue.ps1` | Start hydro queue - może być częścią standardowego workflow | ⚠️ DO POTWIERDZENIA |

**Podsumowanie sekcji A:**
- ✅ **PEWNE do archiwum:** 18 plików
- ⚠️ **DO POTWIERDZENIA:** 10 plików

---

## B. ARCHIVE/OLD_DOCS → Dokumenty zduplikowane, starsze wersje

### B.1 Pliki .txt w root (notatki AWS) - DO POTWIERDZENIA
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `AWS_EMERGENCY_FIX.txt` | Notatka emergency fix - wymaga weryfikacji czy informacje przeniesione do `docs/aws_test/` | ⚠️ DO POTWIERDZENIA |
| `aws_minimal_setup.txt` | Notatka minimal setup - wymaga weryfikacji czy zastąpiona przez `docs/aws_test/` | ⚠️ DO POTWIERDZENIA |
| `AWS_RECOMMENDED_ACTION.txt` | Notatka recommended action - wymaga weryfikacji czy zastąpiona | ⚠️ DO POTWIERDZENIA |
| `AWS_ROUND2_COMMANDS.txt` | Notatka commands round 2 - wymaga weryfikacji czy zastąpiona | ⚠️ DO POTWIERDZENIA |
| `HYDRO_SETUP_COMPLETE.txt` | Notatka hydro setup - wymaga weryfikacji czy zastąpiona przez `docs/` | ⚠️ DO POTWIERDZENIA |

### B.2 Pliki struktury (DO POTWIERDZENIA)
| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `docs_structure.txt` | Struktura dokumentacji - wymaga weryfikacji czy zastąpiony przez `docs/INDEX.md` lub `docs/NAVIGATION_GUIDE.md` | ⚠️ DO POTWIERDZENIA |

**Podsumowanie sekcji B:**
- ⚠️ **DO POTWIERDZENIA:** 6 plików

**Uwaga:** Wszystkie pliki .txt wymagają weryfikacji czy:
1. Informacje zostały przeniesione do właściwej dokumentacji w `docs/`
2. Pliki są już nieaktualne
3. Nie są używane jako szybkie referencje

---

## C. ARCHIVE/EXPERIMENTS → Prototypy, alternatywne pipeline'y

### C.1 Katalogi eksperymentalne (DO POTWIERDZENIA)
| Katalog | Uzasadnienie | Status |
|---------|--------------|--------|
| `configs/` (root) | 24 pliki YAML - wymaga weryfikacji czy to eksperymentalne konfiguracje (testy, optymalizacje) czy aktywnie używane. Jeśli `aws_test/configs/` są produkcyjne, to `configs/` może być eksperymentalne | ⚠️ DO POTWIERDZENIA |
| `diagnostics/` (root) | Katalog diagnostyczny - wymaga weryfikacji czy to eksperymentalna diagnostyka czy ma różne przeznaczenie niż `backend/diagnostics/` | ⚠️ DO POTWIERDZENIA |

**Podsumowanie sekcji C:**
- ⚠️ **DO POTWIERDZENIA:** 2 katalogi

**Uwaga:** Oba katalogi wymagają weryfikacji:
1. Czy są aktywnie używane
2. Czy mają unikalne funkcje
3. Czy są eksperymentami czy częścią produkcyjnego workflow

---

## D. ARCHIVE/TMP_RESULTS → Dane tymczasowe / nieużywane

### D.1 Katalogi wyników AWS (DO POTWIERDZENIA)
| Katalog | Uzasadnienie | Status |
|---------|--------------|--------|
| `aws_results/` (root) | Zawiera `miller_urey_extended/` i `run_1/` do `run_18/` - wymaga weryfikacji czy to duplikat `aws_test/results/` lub `results/phase2b_additional/`. Jeśli duplikat → archiwum. Jeśli unikalne dane → NIE RUSZAĆ | ⚠️ DO POTWIERDZENIA |
| `aws_test/results_16_completed/` | Archiwum 16 ukończonych runów - wymaga weryfikacji czy już nieużywane do analizy | ⚠️ DO POTWIERDZENIA |
| `aws_test/results_28_completed/` | Archiwum 28 ukończonych runów - wymaga weryfikacji czy już nieużywane do analizy | ⚠️ DO POTWIERDZENIA |
| `aws_test/results_all_completed/` | Archiwum wszystkich ukończonych - wymaga weryfikacji czy już nieużywane do analizy | ⚠️ DO POTWIERDZENIA |

**Uwaga:** Jeśli te katalogi są duże (>2MB), rozważyć przeniesienie poza repo z README wskazującym lokalizację.

**Podsumowanie sekcji D:**
- ⚠️ **DO POTWIERDZENIA:** 4 katalogi

**Uwaga:** Wszystkie katalogi wyników wymagają weryfikacji:
1. Czy zawierają unikalne dane czy są duplikatami
2. Czy są używane do analizy/publikacji
3. Czy można je bezpiecznie zarchiwizować

---

## E. ARCHIVE/DEPRECATED → Kod/konfiguracje zastąpione finalną wersją

### E.1 Katalog `live2.0/` (DO POTWIERDZENIA - WYSOKI PRIORYTET)
| Katalog | Uzasadnienie | Status |
|---------|--------------|--------|
| `live2.0/` (root) | Zagnieżdżony katalog projektu - wygląda na błąd agenta. Zawiera `backend/` i inne podkatalogi. Wymaga **pilnej weryfikacji** czy to duplikat całego projektu czy ma unikalne pliki. Jeśli duplikat → archiwum. Jeśli ma unikalne pliki → wymaga decyzji | ⚠️ DO POTWIERDZENIA (WYSOKI PRIORYTET) |

**Podsumowanie sekcji E:**
- ⚠️ **DO POTWIERDZENIA:** 1 katalog (WYSOKI PRIORYTET)

---

## 📊 PODSUMOWANIE PLANU

### Statystyki:
- ✅ **PEWNE do archiwum:** 18 plików (sekcja A.1-A.4)
- ⚠️ **DO POTWIERDZENIA:** 23 elementy (pliki + katalogi)

### Rozkład kategorii:
- **A. one_off_scripts:** 18 pewnych + 10 do potwierdzenia = 28 plików
- **B. old_docs:** 6 plików do potwierdzenia
- **C. experiments:** 2 katalogi do potwierdzenia
- **D. tmp_results:** 4 katalogi do potwierdzenia
- **E. deprecated:** 1 katalog do potwierdzenia (wysoki priorytet)

### Priorytety weryfikacji:
1. **WYSOKI:** `live2.0/` (root) - może być duplikat całego projektu
2. **ŚREDNI:** `aws_results/` (root) - może być duplikat wyników
3. **ŚREDNI:** `configs/` (root) vs `aws_test/configs/` - weryfikacja przeznaczenia
4. **NISKI:** Skrypty start/utility - weryfikacja czy część workflow

---

## ⚠️ WAŻNE UWAGI

1. **Przed wykonaniem:** Wszystkie elementy oznaczone jako "DO POTWIERDZENIA" wymagają weryfikacji przez użytkownika
2. **Read-Only Zones:** Żaden element z read-only zones nie jest proponowany do archiwizacji
3. **Git mv:** Wszystkie operacje będą wykonane przez `git mv` aby zachować historię
4. **ARCHIVE_LOG.md:** Po wykonaniu planu zostanie utworzony `archive/ARCHIVE_LOG.md` z wpisami

---

## ❓ CZY ZATWIERDZASZ PLAN?

**Plan zawiera:**
- 18 plików **PEWNYCH** do archiwizacji (sekcja A.1-A.4)
- 23 elementy **DO POTWIERDZENIA** przed archiwizacją

**Proponowane działanie:**
1. Zatwierdź 18 pewnych plików → wykonanie `git mv` dla sekcji A.1-A.4
2. Zweryfikuj elementy "DO POTWIERDZENIA" → decyzja o każdym z osobna
3. Po weryfikacji → wykonanie `git mv` dla zatwierdzonych elementów

**Czy zatwierdzasz plan?**

