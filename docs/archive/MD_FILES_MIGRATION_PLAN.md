---
date: 2025-11-23
label: plan
---

# PLAN PRZENIESIENIA PLIKÓW .MD POZA `docs/`

**Data:** 2025-11-23  
**Zasady:** Zgodnie z `docs/ARCHIVE_POLICY.md` i `docs/NAVIGATION_GUIDE.md`  
**Cel:** Wszystkie dokumenty .md powinny być w `docs/` lub `archive/old_docs/` jeśli zastąpione

---

## ⚠️ READ-ONLY ZONES (NIE DOTYKAMY)

Zgodnie z `ARCHIVE_POLICY.md`, **NIE PROponujemy** przenoszenia:
- `docs/phase2b/**`
- `docs/technical/**`
- `results/phase2b_additional/**` (ukończone runy)

---

## 📋 ANALIZA PLIKÓW .MD POZA `docs/`

### Kategorie plików:

1. **README.md** - standardowe pliki README (zostają na miejscu)
2. **Dokumentacja techniczna** - powinna być w `docs/technical/`
3. **Dokumentacja AWS** - powinna być w `docs/aws_test/`
4. **Analizy wyników** - powinny być w `docs/analysis/` lub `docs/phase2b/`
5. **Raporty wyników** - mogą być w `docs/` lub `archive/old_docs/` jeśli zastąpione
6. **Dokumenty paper** - zostają w `paper/` (to jest katalog publikacji)
7. **Dokumentacja CI/CD** - może być w `docs/` lub zostaje w `.github/`

---

## A. PLIKI DO ZOSTAWIENIA NA MIEJSCU (OK)

### A.1 Standardowe README (zostają)

| Plik | Lokalizacja | Uzasadnienie |
|------|-------------|--------------|
| `README.md` | Root | Standardowy README projektu - zostaje na miejscu |
| `archive/README.md` | archive/ | README archiwum - zostaje |
| `archive/deprecated/README.md` | archive/deprecated/ | README deprecated - zostaje |
| `archive/experiments/README.md` | archive/experiments/ | README experiments - zostaje |
| `archive/old_docs/README.md` | archive/old_docs/ | README old_docs - zostaje |
| `archive/one_off_scripts/README.md` | archive/one_off_scripts/ | README one_off_scripts - zostaje |
| `archive/tmp_results/README.md` | archive/tmp_results/ | README tmp_results - zostaje |
| `paper/README.md` | paper/ | README katalogu paper - zostaje |
| `scripts/README.md` | scripts/ | README katalogu scripts - zostaje (może być OK) |
| `analysis/phase2b_miller_urey/README.md` | analysis/ | README katalogu analizy - zostaje (może być OK) |
| `.github/workflows/README.md` | .github/workflows/ | README GitHub workflows - zostaje (standardowa lokalizacja) |
| `backend/.pytest_cache/README.md` | backend/.pytest_cache/ | README cache - zostaje (cache, nie dotykamy) |

**Podsumowanie A.1:** 12 plików README - **ZOSTAJĄ**

### A.2 Dokumenty paper (zostają w `paper/`)

| Plik | Lokalizacja | Uzasadnienie |
|------|-------------|--------------|
| `paper/CONCLUSIONS_STRUCTURE.md` | paper/ | Struktura conclusions dla publikacji - zostaje |
| `paper/DISCUSSION_STRUCTURE.md` | paper/ | Struktura discussion dla publikacji - zostaje |
| `paper/EXTENDED_SESSION_COMPLETE.md` | paper/ | Podsumowanie sesji dla publikacji - zostaje |
| `paper/FINAL_SESSION_SUMMARY.md` | paper/ | Finalne podsumowanie dla publikacji - zostaje |
| `paper/INTRODUCTION_REVIEW.md` | paper/ | Review introduction dla publikacji - zostaje |
| `paper/INTRODUCTION_SESSION_SUMMARY.md` | paper/ | Podsumowanie introduction dla publikacji - zostaje |
| `paper/METHODS_REVIEW.md` | paper/ | Review methods dla publikacji - zostaje |
| `paper/PIPELINE_QUICK_REFERENCE.md` | paper/ | Quick reference dla publikacji - zostaje |
| `paper/QUANTUM_AI_EXPANSION_ANALYSIS.md` | paper/ | Analiza expansion dla publikacji - zostaje |
| `paper/RESULTS_STRUCTURE.md` | paper/ | Struktura results dla publikacji - zostaje |
| `paper/SESSION_SUMMARY.md` | paper/ | Podsumowanie sesji dla publikacji - zostaje |
| `paper/TIER1_IMPLEMENTATION_GUIDE.md` | paper/ | Guide implementacji dla publikacji - zostaje |
| `paper/TODAY_COMPLETE_SUMMARY.md` | paper/ | Podsumowanie dzisiejsze dla publikacji - zostaje |
| `paper/WORK_PLAN.md` | paper/ | Plan pracy dla publikacji - zostaje |

**Uzasadnienie:** Katalog `paper/` jest dedykowanym katalogiem dla dokumentów publikacji - zgodnie z ARCHIVE_POLICY, dokumenty naukowe związane z publikacją mogą być w dedykowanych katalogach.

**Podsumowanie A.2:** 14 plików paper - **ZOSTAJĄ**

---

## B. PLIKI DO PRZENIESIENIA DO `docs/` (DOKUMENTACJA)

### B.1 Dokumentacja AWS (do `docs/aws_test/`)

| Plik | Nowa lokalizacja | Uzasadnienie |
|------|------------------|--------------|
| `aws_test/scripts/CPU_THREADS_GUIDE.md` | `docs/aws_test/CPU_THREADS_GUIDE.md` | Dokumentacja AWS - powinna być w `docs/aws_test/` zgodnie z NAVIGATION_GUIDE |
| `aws_test/scripts/MONITORING_SCRIPTS_COMPARISON.md` | `docs/aws_test/MONITORING_SCRIPTS_COMPARISON.md` | Dokumentacja AWS - powinna być w `docs/aws_test/` |

**Uzasadnienie:** Zgodnie z `docs/NAVIGATION_GUIDE.md` sekcja 3.5, dokumentacja AWS powinna być w `docs/aws_test/`, nie w `aws_test/scripts/`.

**Podsumowanie B.1:** 2 pliki - **DO PRZENIESIENIA do `docs/aws_test/`**

### B.2 Dokumentacja techniczna backend (do `docs/technical/`)

| Plik | Nowa lokalizacja | Uzasadnienie |
|------|------------------|--------------|
| `backend/sim/io/schema.md` | `docs/technical/backend_sim_io_schema.md` | Dokumentacja schematu I/O - dokumentacja techniczna |
| `backend/tests/tests.md` | `docs/technical/backend_tests.md` | Dokumentacja testów backend - dokumentacja techniczna |
| `backend/tests/TEST_SUMMARY.md` | `docs/technical/backend_tests_summary.md` | Podsumowanie testów - dokumentacja techniczna |

**Uzasadnienie:** Zgodnie z `docs/NAVIGATION_GUIDE.md` sekcja 3.5, dokumentacja techniczna powinna być w `docs/technical/`.

**Podsumowanie B.2:** 3 pliki - **DO PRZENIESIENIA do `docs/technical/`**

### B.3 Dokumentacja CI/CD (do `docs/` lub zostaje)

| Plik | Nowa lokalizacja | Uzasadnienie | Status |
|------|------------------|--------------|--------|
| `.github/CI_CHEATSHEET.md` | `docs/CI_CHEATSHEET.md` | Cheatsheet CI/CD - może być w `docs/` | ⚠️ DO POTWIERDZENIA |

**Uzasadnienie:** Dokumentacja CI/CD może być w `docs/` lub zostaje w `.github/` (standardowa lokalizacja dla GitHub).

**Podsumowanie B.3:** 1 plik - **DO POTWIERDZENIA**

### B.4 Analizy wyników (do `docs/analysis/` lub `docs/phase2b/`)

| Plik | Nowa lokalizacja | Uzasadnienie |
|------|------------------|--------------|
| `analysis/phase2b_miller_urey/ANALYSIS_FINDINGS.md` | `docs/analysis/phase2b_miller_urey_ANALYSIS_FINDINGS.md` | Analiza wyników Phase 2B - powinna być w `docs/analysis/` |
| `analysis/phase2b_miller_urey/PAPER_SUMMARY.md` | `docs/analysis/phase2b_miller_urey_PAPER_SUMMARY.md` | Podsumowanie dla publikacji - powinna być w `docs/analysis/` |

**Uzasadnienie:** Analizy wyników powinny być w `docs/analysis/` zgodnie z NAVIGATION_GUIDE.

**Podsumowanie B.4:** 2 pliki - **DO PRZENIESIENIA do `docs/analysis/`**

### B.5 Dokumentacja skryptów (do `docs/` lub zostaje)

| Plik | Nowa lokalizacja | Uzasadnienie | Status |
|------|------------------|--------------|--------|
| `scripts/ANALYSIS_QUICK_REF.md` | `docs/scripts_ANALYSIS_QUICK_REF.md` | Quick reference analizy - może być w `docs/` | ⚠️ DO POTWIERDZENIA |

**Uzasadnienie:** Dokumentacja skryptów może być w `docs/` lub zostaje w `scripts/` (może być OK jako lokalna dokumentacja).

**Podsumowanie B.5:** 1 plik - **DO POTWIERDZENIA**

---

## C. PLIKI DO PRZENIESIENIA DO `archive/old_docs/` (ZASTĄPIONE)

### C.1 Raporty wyników w `results/` (DO POTWIERDZENIA - zgodnie z ARCHIVE_POLICY 4.2)

| Plik | Nowa lokalizacja | Uzasadnienie | Status |
|------|------------------|--------------|--------|
| `results/FINAL_ETA_REPORT.md` | `archive/old_docs/results_FINAL_ETA_REPORT.md` | Raport ETA - wymaga weryfikacji czy zastąpiony przez nowszy raport | ⚠️ DO POTWIERDZENIA |
| `results/FINAL_OPTIMIZATION_REPORT.md` | `archive/old_docs/results_FINAL_OPTIMIZATION_REPORT.md` | Raport optymalizacji - wymaga weryfikacji czy zastąpiony | ⚠️ DO POTWIERDZENIA |
| `results/OPTIMIZATION_SUMMARY.md` | `archive/old_docs/results_OPTIMIZATION_SUMMARY.md` | Podsumowanie optymalizacji - wymaga weryfikacji czy zastąpione | ⚠️ DO POTWIERDZENIA |
| `results/phase2b_additional/ROZWIAZANIE_KATALOG.md` | `archive/old_docs/results_phase2b_ROZWIAZANIE_KATALOG.md` | Rozwiązanie katalogu - wymaga weryfikacji czy zastąpione | ⚠️ DO POTWIERDZENIA |
| `results/phase2b_additional/run_1_ANALYSIS.md` | `archive/old_docs/results_phase2b_run_1_ANALYSIS.md` | Analiza run_1 - wymaga weryfikacji czy zastąpiona przez nowszą analizę | ⚠️ DO POTWIERDZENIA |
| `results/phase2b_additional/STATUS_SUMMARY.md` | `archive/old_docs/results_phase2b_STATUS_SUMMARY.md` | Podsumowanie statusu - wymaga weryfikacji czy zastąpione | ⚠️ DO POTWIERDZENIA |
| `results/phase2b_aws_results/FINAL_STATUS.md` | `archive/old_docs/results_phase2b_aws_FINAL_STATUS.md` | Finalny status - wymaga weryfikacji czy zastąpiony | ⚠️ DO POTWIERDZENIA |
| `results/phase2b_aws_results/phase2b_analysis_report.md` | `archive/old_docs/results_phase2b_aws_analysis_report.md` | Raport analizy - wymaga weryfikacji czy zastąpiony | ⚠️ DO POTWIERDZENIA |
| `results/phase2b_aws_results/phase2b_summary_report.md` | `archive/old_docs/results_phase2b_aws_summary_report.md` | Raport podsumowania - wymaga weryfikacji czy zastąpiony | ⚠️ DO POTWIERDZENIA |
| `results/phase2b_aws_results/STATUS.md` | `archive/old_docs/results_phase2b_aws_STATUS.md` | Status - wymaga weryfikacji czy zastąpiony | ⚠️ DO POTWIERDZENIA |
| `results/phase2b_aws_results/formamide_debug/formamide_debug_report.md` | `archive/old_docs/results_phase2b_aws_formamide_debug_report.md` | Raport debug formamide - wymaga weryfikacji czy zastąpiony | ⚠️ DO POTWIERDZENIA |
| `results/spatial_hash_test/PERFORMANCE_REPORT.md` | `archive/old_docs/results_spatial_hash_PERFORMANCE_REPORT.md` | Raport wydajności spatial hash - wymaga weryfikacji czy zastąpiony | ⚠️ DO POTWIERDZENIA |
| `results/test_formamide_10k/TEST_REPORT.md` | `archive/old_docs/results_test_formamide_10k_TEST_REPORT.md` | Raport testu formamide - wymaga weryfikacji czy zastąpiony | ⚠️ DO POTWIERDZENIA |

**Uzasadnienie:** Raporty wyników w `results/` mogą być:
1. Zastąpione przez nowsze raporty w `docs/` → `archive/old_docs/`
2. Aktywnie używane → przenieść do `docs/analysis/` lub `docs/phase2b/`
3. Częścią wyników Phase 2B → **NIE DOTYKAMY** (read-only zone)

**Uwaga:** `results/phase2b_additional/**` jest w read-only zone - wymaga szczególnej ostrożności.

**Podsumowanie C.1:** 13 plików - **DO POTWIERDZENIA** (wymaga weryfikacji czy zastąpione)

---

## D. PLIKI W READ-ONLY ZONES (NIE DOTYKAMY)

### D.1 Wyniki Phase 2B (NIE DOTYKAMY - zgodnie z ARCHIVE_POLICY sekcja 3)

| Plik | Lokalizacja | Status |
|------|-------------|--------|
| `results/phase2b_additional/ROZWIAZANIE_KATALOG.md` | results/phase2b_additional/ | ⚠️ **READ-ONLY ZONE** - wymaga weryfikacji przed przeniesieniem |
| `results/phase2b_additional/run_1_ANALYSIS.md` | results/phase2b_additional/ | ⚠️ **READ-ONLY ZONE** - wymaga weryfikacji przed przeniesieniem |
| `results/phase2b_additional/STATUS_SUMMARY.md` | results/phase2b_additional/ | ⚠️ **READ-ONLY ZONE** - wymaga weryfikacji przed przeniesieniem |

**Uzasadnienie:** Zgodnie z `ARCHIVE_POLICY.md` sekcja 3, `results/**` jest read-only. Przed przeniesieniem wymaga weryfikacji czy to dokumenty zastąpione czy aktywnie używane.

---

## 📊 PODSUMOWANIE PLANU

### Pliki do zostawienia (OK):
- **12 plików README** - standardowe README, zostają na miejscu
- **14 plików paper** - dokumenty publikacji, zostają w `paper/`

### Pliki do przeniesienia do `docs/`:
- **2 pliki AWS** → `docs/aws_test/`
- **3 pliki techniczne** → `docs/technical/`
- **2 pliki analiz** → `docs/analysis/`
- **2 pliki DO POTWIERDZENIA** → `docs/` (CI/CD, scripts)

### Pliki do przeniesienia do `archive/old_docs/`:
- **13 plików raportów** → `archive/old_docs/` (DO POTWIERDZENIA - wymaga weryfikacji czy zastąpione)

### Statystyki:
- ✅ **ZOSTAJĄ:** 26 plików (12 README + 14 paper)
- ✅ **DO PRZENIESIENIA do docs/:** 7 plików (5 pewnych + 2 do potwierdzenia)
- ⚠️ **DO PRZENIESIENIA do archive/old_docs/:** 13 plików (wszystkie DO POTWIERDZENIA)
- ⚠️ **READ-ONLY ZONE:** 3 pliki (wymagają weryfikacji)

---

## 🔧 SZCZEGÓŁOWY PLAN PRZENIESIENIA

### B.1 Dokumentacja AWS → `docs/aws_test/` (PEWNE)

| Operacja | Ścieżka źródłowa | Nowa lokalizacja |
|----------|-------------------|------------------|
| `git mv` | `aws_test/scripts/CPU_THREADS_GUIDE.md` | `docs/aws_test/CPU_THREADS_GUIDE.md` |
| `git mv` | `aws_test/scripts/MONITORING_SCRIPTS_COMPARISON.md` | `docs/aws_test/MONITORING_SCRIPTS_COMPARISON.md` |

**Uzasadnienie:** Zgodnie z NAVIGATION_GUIDE, dokumentacja AWS powinna być w `docs/aws_test/`.

### B.2 Dokumentacja techniczna → `docs/technical/` (PEWNE)

| Operacja | Ścieżka źródłowa | Nowa lokalizacja |
|----------|-------------------|------------------|
| `git mv` | `backend/sim/io/schema.md` | `docs/technical/backend_sim_io_schema.md` |
| `git mv` | `backend/tests/tests.md` | `docs/technical/backend_tests.md` |
| `git mv` | `backend/tests/TEST_SUMMARY.md` | `docs/technical/backend_tests_summary.md` |

**Uzasadnienie:** Zgodnie z NAVIGATION_GUIDE, dokumentacja techniczna powinna być w `docs/technical/`.

### B.4 Analizy wyników → `docs/analysis/` (PEWNE)

| Operacja | Ścieżka źródłowa | Nowa lokalizacja |
|----------|-------------------|------------------|
| `git mv` | `analysis/phase2b_miller_urey/ANALYSIS_FINDINGS.md` | `docs/analysis/phase2b_miller_urey_ANALYSIS_FINDINGS.md` |
| `git mv` | `analysis/phase2b_miller_urey/PAPER_SUMMARY.md` | `docs/analysis/phase2b_miller_urey_PAPER_SUMMARY.md` |

**Uzasadnienie:** Analizy wyników powinny być w `docs/analysis/` zgodnie z NAVIGATION_GUIDE.

### B.3, B.5 Dokumentacja CI/CD i scripts (DO POTWIERDZENIA)

| Operacja | Ścieżka źródłowa | Nowa lokalizacja | Status |
|----------|-------------------|------------------|--------|
| `git mv` | `.github/CI_CHEATSHEET.md` | `docs/CI_CHEATSHEET.md` | ⚠️ DO POTWIERDZENIA |
| `git mv` | `scripts/ANALYSIS_QUICK_REF.md` | `docs/scripts_ANALYSIS_QUICK_REF.md` | ⚠️ DO POTWIERDZENIA |

**Uzasadnienie:** Wymaga weryfikacji czy te pliki powinny być w `docs/` czy zostają na miejscu.

### C.1 Raporty wyników → `archive/old_docs/` (DO POTWIERDZENIA)

**Uwaga:** Wszystkie pliki w `results/` wymagają weryfikacji przed przeniesieniem:
1. Czy są zastąpione przez nowsze raporty w `docs/`?
2. Czy są aktywnie używane?
3. Czy są częścią read-only zone (`results/phase2b_additional/**`)?

**Proponowane przeniesienie (TYLKO jeśli zastąpione):**

| Operacja | Ścieżka źródłowa | Nowa lokalizacja | Status |
|----------|-------------------|------------------|--------|
| `git mv` | `results/FINAL_ETA_REPORT.md` | `archive/old_docs/results_FINAL_ETA_REPORT.md` | ⚠️ DO POTWIERDZENIA |
| `git mv` | `results/FINAL_OPTIMIZATION_REPORT.md` | `archive/old_docs/results_FINAL_OPTIMIZATION_REPORT.md` | ⚠️ DO POTWIERDZENIA |
| `git mv` | `results/OPTIMIZATION_SUMMARY.md` | `archive/old_docs/results_OPTIMIZATION_SUMMARY.md` | ⚠️ DO POTWIERDZENIA |
| `git mv` | `results/spatial_hash_test/PERFORMANCE_REPORT.md` | `archive/old_docs/results_spatial_hash_PERFORMANCE_REPORT.md` | ⚠️ DO POTWIERDZENIA |
| `git mv` | `results/test_formamide_10k/TEST_REPORT.md` | `archive/old_docs/results_test_formamide_10k_TEST_REPORT.md` | ⚠️ DO POTWIERDZENIA |

**UWAGA:** Pliki w `results/phase2b_additional/` i `results/phase2b_aws_results/` są w read-only zone - **NIE PROponujemy** przeniesienia bez wyraźnej weryfikacji.

---

## ⚠️ WAŻNE UWAGI

1. **Read-Only Zones:** Pliki w `results/phase2b_additional/**` są w read-only zone - wymagają szczególnej ostrożności
2. **Weryfikacja przed przeniesieniem:** Wszystkie pliki w `results/` wymagają weryfikacji:
   - Czy są zastąpione przez nowsze raporty?
   - Czy są aktywnie używane?
   - Czy można je bezpiecznie zarchiwizować?
3. **Git mv:** Wszystkie operacje będą wykonane przez `git mv` aby zachować historię
4. **ARCHIVE_LOG.md:** Po wykonaniu planu zostanie zaktualizowany log archiwizacji

---

## 📊 STATYSTYKI FINALNE

### Zostają na miejscu:
- **12 plików README** (standardowe)
- **14 plików paper** (dokumenty publikacji)

### Do przeniesienia:
- **7 plików do `docs/`** (5 pewnych + 2 do potwierdzenia)
- **13 plików do `archive/old_docs/`** (wszystkie DO POTWIERDZENIA)

### Priorytety:
1. **PEWNE:** 7 plików do `docs/` (B.1, B.2, B.4)
2. **DO POTWIERDZENIA:** 2 pliki do `docs/` (B.3, B.5)
3. **DO POTWIERDZENIA:** 13 plików do `archive/old_docs/` (C.1 - wymaga weryfikacji czy zastąpione)

---

## ❓ CZY ZATWIERDZASZ PLAN?

**Plan zawiera:**
- ✅ **7 plików PEWNYCH** do przeniesienia do `docs/` (B.1, B.2, B.4)
- ⚠️ **2 pliki DO POTWIERDZENIA** do przeniesienia do `docs/` (B.3, B.5)
- ⚠️ **13 plików DO POTWIERDZENIA** do przeniesienia do `archive/old_docs/` (C.1)

**Proponowane działanie:**
1. **Zatwierdź 7 pewnych plików** → wykonanie `git mv` dla B.1, B.2, B.4
2. **Zweryfikuj pliki DO POTWIERDZENIA** → decyzja o każdym z osobna
3. **Po weryfikacji** → wykonanie `git mv` dla zatwierdzonych elementów

**Czy zatwierdzasz plan?**

