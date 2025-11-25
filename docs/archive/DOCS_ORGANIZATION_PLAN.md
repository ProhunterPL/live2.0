---
date: 2025-11-23
label: plan
---

# PLAN SEGREGACJI PLIKÓW .MD W `docs/`

**Data:** 2025-11-23  
**Cel:** Uporządkowanie plików .md bezpośrednio w `docs/` poprzez segregację do odpowiednich podkatalogów lub archiwum

---

## 📊 ANALIZA PLIKÓW W `docs/`

Znaleziono **~100+ plików .md** bezpośrednio w `docs/`. Pliki zostały skategoryzowane na podstawie nazw i przeznaczenia.

---

## 📁 PROPOZOWANE NOWE PODKATALOGI

### 1. `docs/plans/` - Plany i roadmapy
### 2. `docs/optimization/` - Optymalizacje i performance
### 3. `docs/fixes/` - Fixy i rozwiązania problemów
### 4. `docs/archive/` - Plany archiwizacji (wykonane)

---

## A. PLIKI DO PRZENIESIENIA DO ISTNIEJĄCYCH PODKATALOGÓW

### A.1 → `docs/sessions/` (już istnieje)

| Plik | Uzasadnienie |
|------|--------------|
| `SESSION_SUMMARY_OCT13.md` | Podsumowanie sesji - powinno być w sessions/ |
| `SESSION_FINAL_OCT13_EVENING.md` | Podsumowanie sesji - powinno być w sessions/ |
| `SESSION_SUMMARY_2025-10-13.md` | Podsumowanie sesji z datą - powinno być w sessions/ |
| `SESSION_SUMMARY_2025-10-16_AWS_PIPELINE.md` | Podsumowanie sesji z datą - powinno być w sessions/ |
| `SESSION_SUMMARY_NOV8_2025.md` | Podsumowanie sesji z datą - powinno być w sessions/ |
| `FINAL_SESSION_SUMMARY_OCT13.md` | Finalne podsumowanie sesji - powinno być w sessions/ |
| `README_SESSION_OCT13.md` | README sesji - powinno być w sessions/ |
| `FINAL_PROGRESS_OCT13.md` | Raport postępu z sesji - powinno być w sessions/ |
| `TODAYS_FINAL_SUMMARY.md` | Podsumowanie dzisiejsze - powinno być w sessions/ |
| `WEEK2_DAY1_SUMMARY.md` | Podsumowanie tygodnia - powinno być w sessions/ |
| `WEEK4_COMPLETION.md` | Podsumowanie tygodnia - powinno być w sessions/ |

**Podsumowanie A.1:** 11 plików → `docs/sessions/`

### A.2 → `docs/aws_test/` (już istnieje)

| Plik | Uzasadnienie |
|------|--------------|
| `AWS_RESULTS_PIPELINE.md` | Dokumentacja AWS - powinna być w aws_test/ |
| `AWS_RUN_STATUS.md` | Status AWS - powinien być w aws_test/ |
| `CHECK_AWS_STATUS.md` | Check AWS - powinien być w aws_test/ |
| `DEBUG_AWS_PHASE2B.md` | Debug AWS Phase 2B - powinien być w aws_test/ |
| `FINAL_AWS_INSTRUCTIONS.md` | Instrukcje AWS - powinny być w aws_test/ |
| `FIX_AWS_SCRIPTS.md` | Fix AWS scripts - powinien być w aws_test/ |
| `INSTRUKCJA_AWS_UPD.md` | Instrukcja AWS - powinna być w aws_test/ |

**Podsumowanie A.2:** 7 plików → `docs/aws_test/`

### A.3 → `docs/troubleshooting/` (już istnieje)

| Plik | Uzasadnienie |
|------|--------------|
| `CLUSTER_DETECTION_ISSUE.md` | Problem z wykrywaniem klastrów - powinien być w troubleshooting/ |
| `CRASH_REPORT.md` | Raport crash - powinien być w troubleshooting/ |
| `FAQ_KLASTRY.md` | FAQ o klastrach - powinien być w troubleshooting/ |
| `GPU_MEMORY_ISSUE.md` | Problem z pamięcią GPU - powinien być w troubleshooting/ |
| `FIX_PATH_PROBLEM.md` | Fix problemu ścieżki - powinien być w troubleshooting/ |
| `FIX_PHYSICS_DATABASE.md` | Fix bazy fizyki - powinien być w troubleshooting/ |
| `WEBSOCKET_FIXES.md` | Fixy WebSocket - powinny być w troubleshooting/ |
| `MATCHER_BUTTON_FIX.md` | Fix przycisku matcher - powinien być w troubleshooting/ |
| `RUNTIME_TIMER_FIX.md` | Fix timera runtime - powinien być w troubleshooting/ |
| `MEMORY_PERFORMANCE_FIX.md` | Fix wydajności pamięci - powinien być w troubleshooting/ |

**Podsumowanie A.3:** 10 plików → `docs/troubleshooting/`

### A.4 → `docs/technical/` (już istnieje)

| Plik | Uzasadnienie |
|------|--------------|
| `PHYSICS_DATABASE.md` | Dokumentacja bazy fizyki - powinna być w technical/ |
| `PHYSICS_DB_INTEGRATION.md` | Integracja bazy fizyki - powinna być w technical/ |
| `THERMODYNAMIC_VALIDATION.md` | Walidacja termodynamiczna - powinna być w technical/ |
| `SCIENTIFIC_INTEGRITY_VERIFICATION.md` | Weryfikacja integralności naukowej - powinna być w technical/ |
| `SKAD_METRICS_BIERZE_CLUSTERS.md` | Dokumentacja metryk - powinna być w technical/ |
| `SPATIAL_HASHING_SUCCESS.md` | Dokumentacja spatial hashing - powinna być w technical/ |

**Podsumowanie A.4:** 6 plików → `docs/technical/`

### A.5 → `docs/analysis/` (już istnieje)

| Plik | Uzasadnienie |
|------|--------------|
| `ANALIZA_WYNIKOW_TESTU.md` | Analiza wyników testu - powinna być w analysis/ |
| `CHECK_REAL_RESULTS.md` | Check wyników - powinien być w analysis/ |
| `PERFORMANCE_ANALYSIS.md` | Analiza wydajności - powinna być w analysis/ |
| `PERFORMANCE_DIAGNOSIS_FINAL.md` | Diagnoza wydajności - powinna być w analysis/ |
| `PHASE2_RESULTS_ASSESSMENT.md` | Ocena wyników Phase 2 - powinna być w analysis/ |
| `TEST_RESULTS_OCT13.md` | Wyniki testów - powinny być w analysis/ |

**Podsumowanie A.5:** 6 plików → `docs/analysis/`

### A.6 → `docs/guides/` (już istnieje)

| Plik | Uzasadnienie |
|------|--------------|
| `QUICK_START.md` | Quick start guide - powinien być w guides/ |
| `QUICK_START_PHASE2.md` | Quick start Phase 2 - powinien być w guides/ |
| `ENVIRONMENT_SETUP.md` | Setup środowiska - powinien być w guides/ |
| `INSTALLATION.md` | Instrukcja instalacji - powinna być w guides/ |
| `DIAGNOSTICS_QUICKSTART.md` | Quickstart diagnostyki - powinien być w guides/ |
| `DIAGNOSTICS.md` | Dokumentacja diagnostyki - powinna być w guides/ |
| `PHASE2_USAGE_GUIDE.md` | Przewodnik Phase 2 - powinien być w guides/ |
| `PHASE3_ANALYSIS_GUIDE.md` | Przewodnik analizy Phase 3 - powinien być w guides/ |
| `RUN_LOCAL_RTX5070.md` | Przewodnik uruchomienia lokalnego - powinien być w guides/ |
| `CLOUD_DEPLOYMENT_GUIDE.md` | Przewodnik deploymentu w chmurze - powinien być w guides/ |

**Podsumowanie A.6:** 10 plików → `docs/guides/`

---

## B. PLIKI DO PRZENIESIENIA DO NOWYCH PODKATALOGÓW

### B.1 → `docs/plans/` (NOWY)

| Plik | Uzasadnienie |
|------|--------------|
| `AGGRESSIVE_OPTIMIZATION_PLAN.md` | Plan optymalizacji - powinien być w plans/ |
| `BOND_ENHANCEMENT_PLAN.md` | Plan ulepszenia wiązań - powinien być w plans/ |
| `CLUSTER_ENHANCEMENT_PLAN.md` | Plan ulepszenia klastrów - powinien być w plans/ |
| `PHASE2_NEXT_STEPS_PL.md` | Następne kroki Phase 2 - powinien być w plans/ |
| `Live2_v1_plan.md` | Plan Live2 v1 - powinien być w plans/ |
| `live2-roadmap.md` | Roadmap Live2 - powinien być w plans/ |
| `ROADMAP_COMPARISON_2025-10-13.md` | Porównanie roadmap - powinno być w plans/ |
| `ROADMAP_UPDATE_2025-10-13.md` | Aktualizacja roadmap - powinna być w plans/ |
| `Live 2-plan walidacji naukowej.md` | Plan walidacji naukowej - powinien być w plans/ |

**Podsumowanie B.1:** 9 plików → `docs/plans/` (NOWY)

### B.2 → `docs/optimization/` (NOWY)

| Plik | Uzasadnienie |
|------|--------------|
| `FPS_OPTIMIZATION.md` | Optymalizacja FPS - powinna być w optimization/ |
| `HYBRID_GPU_CPU_ANALYSIS.md` | Analiza hybrid GPU/CPU - powinna być w optimization/ |
| `HYBRID_GPU_CPU_GUIDE.md` | Przewodnik hybrid GPU/CPU - powinien być w optimization/ |
| `HYBRID_MODE_SUMMARY.md` | Podsumowanie trybu hybrid - powinno być w optimization/ |
| `OPTIMIZACJA_PRODUKCJI.md` | Optymalizacja produkcji - powinna być w optimization/ |
| `OPTIMIZATION_SUMMARY.md` | Podsumowanie optymalizacji - powinno być w optimization/ |
| `OPTIMIZATION_SUMMARY_FINAL.md` | Finalne podsumowanie optymalizacji - powinno być w optimization/ |
| `OPTIMIZE_SIMULATION.md` | Optymalizacja symulacji - powinna być w optimization/ |
| `PERFORMANCE_OPTIMIZATION_14CORES.md` | Optymalizacja wydajności 14 core - powinna być w optimization/ |
| `PERFORMANCE_TUNING.md` | Tuning wydajności - powinien być w optimization/ |
| `PODSUMOWANIE_OPTYMALIZACJI.md` | Podsumowanie optymalizacji - powinno być w optimization/ |
| `PRODUCTION_OPTIMIZATION.md` | Optymalizacja produkcji - powinna być w optimization/ |

**Podsumowanie B.2:** 12 plików → `docs/optimization/` (NOWY)

### B.3 → `docs/fixes/` (NOWY)

| Plik | Uzasadnienie |
|------|--------------|
| `FIX_PATH_PROBLEM.md` | Fix problemu ścieżki - powinien być w fixes/ |
| `FIX_PHYSICS_DATABASE.md` | Fix bazy fizyki - powinien być w fixes/ |
| `FIX_AWS_SCRIPTS.md` | Fix skryptów AWS - powinien być w fixes/ |
| `ROZWIAZANIE_KROK_PO_KROKU.md` | Rozwiązanie krok po kroku - powinno być w fixes/ |
| `PODSUMOWANIE_PROBLEM_KLASTROW.md` | Podsumowanie problemu klastrów - powinno być w fixes/ |

**Uwaga:** Część plików z "FIX" może być już w troubleshooting/ - wymaga weryfikacji czy nie ma duplikatów.

**Podsumowanie B.3:** 5 plików → `docs/fixes/` (NOWY)

### B.4 → `docs/archive/` (NOWY - plany archiwizacji)

| Plik | Uzasadnienie |
|------|--------------|
| `ARCHIVE_PLAN.md` | Plan archiwizacji (wykonany) - powinien być w archive/ |
| `ARCHIVE_PLAN_FINAL.md` | Finalny plan archiwizacji (wykonany) - powinien być w archive/ |
| `DUPLICATES_CANONICAL_PLAN.md` | Plan duplikatów (wykonany) - powinien być w archive/ |
| `MD_FILES_MIGRATION_PLAN.md` | Plan migracji plików MD (wykonany) - powinien być w archive/ |
| `REPO_STRUCTURE_ANALYSIS.md` | Analiza struktury repo - powinna być w archive/ |
| `STRUCTURE_DEEP_ANALYSIS.md` | Głęboka analiza struktury - powinna być w archive/ |

**Uzasadnienie:** Plany archiwizacji są już wykonane i służą jako dokumentacja historyczna.

**Podsumowanie B.4:** 6 plików → `docs/archive/` (NOWY)

---

## C. PLIKI DO ZOSTAWIENIA W `docs/` (GŁÓWNE DOKUMENTY)

### C.1 Dokumenty główne (zostają)

| Plik | Uzasadnienie |
|------|--------------|
| `README.md` | Główny README (jeśli istnieje) - zostaje |
| `INDEX.md` | Indeks dokumentacji - zostaje |
| `NAVIGATION_GUIDE.md` | Przewodnik nawigacji - zostaje (główny dokument) |
| `ARCHIVE_POLICY.md` | Polityka archiwizacji - zostaje (główny dokument) |
| `VALIDATION_ROADMAP.md` | Roadmap walidacji - zostaje (główny dokument) |
| `PROJECT_INDEX.md` | Indeks projektu - zostaje |
| `protocol.md` | Protokół - zostaje |

**Podsumowanie C.1:** 7 plików - **ZOSTAJĄ w docs/**

### C.2 Dokumenty Phase 2 (zostają lub do phase2b/)

| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `PHASE2_DEMO_COMPLETE.md` | Demo Phase 2 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `PHASE2_EXPERIMENTS.md` | Eksperymenty Phase 2 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `PHASE2_INFRASTRUCTURE_READY.md` | Infrastruktura Phase 2 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `PHASE2_INTEGRATION_SUCCESS.md` | Sukces integracji Phase 2 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `PHASE1_COMPLETION_SUMMARY.md` | Podsumowanie Phase 1 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `START_REAL_SIMULATIONS.md` | Start symulacji - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |

**Podsumowanie C.2:** 6 plików - **DO POTWIERDZENIA** (zostają lub do phase2b/)

---

## D. PLIKI DO PRZENIESIENIA DO `archive/old_docs/` (STARE/NIEAKTUALNE)

### D.1 Stare podsumowania i raporty (DO POTWIERDZENIA)

| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `ORGANIZATION_SUMMARY.md` | Stare podsumowanie organizacji - może być zastąpione | ⚠️ DO POTWIERDZENIA |
| `GIT_PUSH_READY.md` | Stary status git push - może być nieaktualny | ⚠️ DO POTWIERDZENIA |
| `REVIEW_TASKS.md` | Stary review zadań - może być nieaktualny | ⚠️ DO POTWIERDZENIA |
| `CI_CD_GUIDE.md` | Stary przewodnik CI/CD - może być zastąpiony przez nowszy | ⚠️ DO POTWIERDZENIA |
| `CI_CHEATSHEET.md` | Stary cheatsheet CI/CD - może być zastąpiony | ⚠️ DO POTWIERDZENIA |

**Uzasadnienie:** Wymaga weryfikacji czy te dokumenty są aktualne czy zastąpione przez nowsze wersje.

**Podsumowanie D.1:** 5 plików - **DO POTWIERDZENIA** → `archive/old_docs/`

### D.2 Dokumenty matcher (DO POTWIERDZENIA)

| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `LIVE 2-matcher.md` | Dokumentacja matcher - może być zastąpiona przez MATCHER_V2.md | ⚠️ DO POTWIERDZENIA |
| `MATCHER_V2.md` | Dokumentacja matcher v2 - może być aktualna | ⚠️ DO POTWIERDZENIA |
| `README_MATCHER.md` | README matcher - może być zastąpione | ⚠️ DO POTWIERDZENIA |

**Uzasadnienie:** Wymaga weryfikacji która wersja jest aktualna.

**Podsumowanie D.2:** 3 pliki - **DO POTWIERDZENIA**

### D.3 Dokumenty quantum/expansion (DO POTWIERDZENIA)

| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `LIVE2_QUANTUM_AI_EXPANSION.md` | Dokumentacja quantum AI expansion - może być nieaktualna | ⚠️ DO POTWIERDZENIA |
| `PATENT_ANALYSIS.md` | Analiza patentu - może być w patents/ | ⚠️ DO POTWIERDZENIA |
| `SCIENTIFIC_OVERVIEW.md` | Przegląd naukowy - może być w technical/ lub patents/ | ⚠️ DO POTWIERDZENIA |

**Podsumowanie D.3:** 3 pliki - **DO POTWIERDZENIA**

### D.4 Dokumenty reorganizacji (DO POTWIERDZENIA)

| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `DOKUMENTACJA_REORGANIZACJA.md` | Dokument reorganizacji - może być w archive/ | ⚠️ DO POTWIERDZENIA |

**Podsumowanie D.4:** 1 plik - **DO POTWIERDZENIA**

---

## E. PLIKI DO PRZENIESIENIA DO `docs/phase2b/` (DO POTWIERDZENIA)

| Plik | Uzasadnienie | Status |
|------|--------------|--------|
| `PHASE2_DEMO_COMPLETE.md` | Demo Phase 2 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `PHASE2_EXPERIMENTS.md` | Eksperymenty Phase 2 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `PHASE2_INFRASTRUCTURE_READY.md` | Infrastruktura Phase 2 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `PHASE2_INTEGRATION_SUCCESS.md` | Sukces integracji Phase 2 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `PHASE1_COMPLETION_SUMMARY.md` | Podsumowanie Phase 1 - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |
| `START_REAL_SIMULATIONS.md` | Start symulacji - może być w phase2b/ | ⚠️ DO POTWIERDZENIA |

**Podsumowanie E:** 6 plików - **DO POTWIERDZENIA** → `docs/phase2b/`

---

## 📊 PODSUMOWANIE PLANU

### Nowe podkatalogi do utworzenia:
1. `docs/plans/` - 9 plików
2. `docs/optimization/` - 12 plików
3. `docs/fixes/` - 5 plików
4. `docs/archive/` - 6 plików

### Przeniesienia do istniejących podkatalogów:
- `docs/sessions/` - 11 plików
- `docs/aws_test/` - 7 plików
- `docs/troubleshooting/` - 10 plików
- `docs/technical/` - 6 plików
- `docs/analysis/` - 6 plików
- `docs/guides/` - 10 plików

### Pliki zostające w docs/:
- 7 plików głównych (INDEX, NAVIGATION_GUIDE, etc.)

### Pliki do potwierdzenia:
- 6 plików Phase 2 → `docs/phase2b/` lub zostają
- 12 plików → `archive/old_docs/` (stare/nieaktualne)

### Statystyki:
- ✅ **PEWNE do przeniesienia:** 61 plików
- ⚠️ **DO POTWIERDZENIA:** 18 plików
- ✅ **ZOSTAJĄ w docs/:** 7 plików

---

## 🔧 SZCZEGÓŁOWY PLAN PRZENIESIENIA

### Krok 1: Utworzenie nowych podkatalogów

```bash
mkdir docs/plans
mkdir docs/optimization
mkdir docs/fixes
mkdir docs/archive
```

### Krok 2: Przeniesienie do istniejących podkatalogów (PEWNE)

**A.1 → docs/sessions/ (11 plików)**
- SESSION_SUMMARY_OCT13.md
- SESSION_FINAL_OCT13_EVENING.md
- SESSION_SUMMARY_2025-10-13.md
- SESSION_SUMMARY_2025-10-16_AWS_PIPELINE.md
- SESSION_SUMMARY_NOV8_2025.md
- FINAL_SESSION_SUMMARY_OCT13.md
- README_SESSION_OCT13.md
- FINAL_PROGRESS_OCT13.md
- TODAYS_FINAL_SUMMARY.md
- WEEK2_DAY1_SUMMARY.md
- WEEK4_COMPLETION.md

**A.2 → docs/aws_test/ (7 plików)**
- AWS_RESULTS_PIPELINE.md
- AWS_RUN_STATUS.md
- CHECK_AWS_STATUS.md
- DEBUG_AWS_PHASE2B.md
- FINAL_AWS_INSTRUCTIONS.md
- FIX_AWS_SCRIPTS.md
- INSTRUKCJA_AWS_UPD.md

**A.3 → docs/troubleshooting/ (10 plików)**
- CLUSTER_DETECTION_ISSUE.md
- CRASH_REPORT.md
- FAQ_KLASTRY.md
- GPU_MEMORY_ISSUE.md
- FIX_PATH_PROBLEM.md
- FIX_PHYSICS_DATABASE.md
- WEBSOCKET_FIXES.md
- MATCHER_BUTTON_FIX.md
- RUNTIME_TIMER_FIX.md
- MEMORY_PERFORMANCE_FIX.md

**A.4 → docs/technical/ (6 plików)**
- PHYSICS_DATABASE.md
- PHYSICS_DB_INTEGRATION.md
- THERMODYNAMIC_VALIDATION.md
- SCIENTIFIC_INTEGRITY_VERIFICATION.md
- SKAD_METRICS_BIERZE_CLUSTERS.md
- SPATIAL_HASHING_SUCCESS.md

**A.5 → docs/analysis/ (6 plików)**
- ANALIZA_WYNIKOW_TESTU.md
- CHECK_REAL_RESULTS.md
- PERFORMANCE_ANALYSIS.md
- PERFORMANCE_DIAGNOSIS_FINAL.md
- PHASE2_RESULTS_ASSESSMENT.md
- TEST_RESULTS_OCT13.md

**A.6 → docs/guides/ (10 plików)**
- QUICK_START.md
- QUICK_START_PHASE2.md
- ENVIRONMENT_SETUP.md
- INSTALLATION.md
- DIAGNOSTICS_QUICKSTART.md
- DIAGNOSTICS.md
- PHASE2_USAGE_GUIDE.md
- PHASE3_ANALYSIS_GUIDE.md
- RUN_LOCAL_RTX5070.md
- CLOUD_DEPLOYMENT_GUIDE.md

### Krok 3: Przeniesienie do nowych podkatalogów (PEWNE)

**B.1 → docs/plans/ (9 plików)**
- AGGRESSIVE_OPTIMIZATION_PLAN.md
- BOND_ENHANCEMENT_PLAN.md
- CLUSTER_ENHANCEMENT_PLAN.md
- PHASE2_NEXT_STEPS_PL.md
- Live2_v1_plan.md
- live2-roadmap.md
- ROADMAP_COMPARISON_2025-10-13.md
- ROADMAP_UPDATE_2025-10-13.md
- Live 2-plan walidacji naukowej.md

**B.2 → docs/optimization/ (12 plików)**
- FPS_OPTIMIZATION.md
- HYBRID_GPU_CPU_ANALYSIS.md
- HYBRID_GPU_CPU_GUIDE.md
- HYBRID_MODE_SUMMARY.md
- OPTIMIZACJA_PRODUKCJI.md
- OPTIMIZATION_SUMMARY.md
- OPTIMIZATION_SUMMARY_FINAL.md
- OPTIMIZE_SIMULATION.md
- PERFORMANCE_OPTIMIZATION_14CORES.md
- PERFORMANCE_TUNING.md
- PODSUMOWANIE_OPTYMALIZACJI.md
- PRODUCTION_OPTIMIZATION.md

**B.3 → docs/fixes/ (5 plików)**
- FIX_PATH_PROBLEM.md
- FIX_PHYSICS_DATABASE.md
- FIX_AWS_SCRIPTS.md
- ROZWIAZANIE_KROK_PO_KROKU.md
- PODSUMOWANIE_PROBLEM_KLASTROW.md

**B.4 → docs/archive/ (6 plików)**
- ARCHIVE_PLAN.md
- ARCHIVE_PLAN_FINAL.md
- DUPLICATES_CANONICAL_PLAN.md
- MD_FILES_MIGRATION_PLAN.md
- REPO_STRUCTURE_ANALYSIS.md
- STRUCTURE_DEEP_ANALYSIS.md

---

## ⚠️ WAŻNE UWAGI

1. **Duplikaty:** Niektóre pliki mogą być wymienione w wielu kategoriach (np. FIX_AWS_SCRIPTS.md w aws_test/ i fixes/) - wymaga weryfikacji przed przeniesieniem
2. **Read-Only Zones:** Pliki związane z Phase 2B mogą być w read-only zone - wymaga weryfikacji
3. **Git mv:** Wszystkie operacje będą wykonane przez `git mv` aby zachować historię
4. **Weryfikacja:** Pliki oznaczone jako "DO POTWIERDZENIA" wymagają weryfikacji przed przeniesieniem

---

## ❓ CZY ZATWIERDZASZ PLAN?

**Plan zawiera:**
- ✅ **61 plików PEWNYCH** do przeniesienia
- ⚠️ **18 plików DO POTWIERDZENIA**
- ✅ **7 plików** zostaje w docs/ (główne dokumenty)
- ✅ **4 nowe podkatalogi** do utworzenia

**Proponowane działanie:**
1. **Utworzenie nowych podkatalogów**
2. **Zatwierdzenie 61 pewnych plików** → wykonanie `git mv`
3. **Weryfikacja plików DO POTWIERDZENIA** → decyzja o każdym z osobna

**Czy zatwierdzasz plan?**

