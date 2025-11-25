---
date: 2025-11-23
label: log
---

# Archive Log - Live 2.0

**Log wszystkich operacji archiwizacji zgodnie z `docs/ARCHIVE_POLICY.md`**

---

## [2025-11-23] Archiwizacja zatwierdzonych plików (sekcja A.1-A.4)

### A.1 Skrypty diagnostyczne/debug

**[2025-11-23] Archived: check_real_clusters.py**  
Reason: One-off debugging script, used once during cluster debugging.  
Moved to: archive/one_off_scripts/check_real_clusters.py

**[2025-11-23] Archived: diagnose_chemistry.py**  
Reason: Diagnostic script for chemistry analysis, one-time use.  
Moved to: archive/one_off_scripts/diagnose_chemistry.py

**[2025-11-23] Archived: diagnose_round1.sh**  
Reason: Diagnostic script for specific round (round1), one-time use.  
Moved to: archive/one_off_scripts/diagnose_round1.sh

**[2025-11-23] Archived: fix_catalog_timeline.py**  
Reason: One-time fix for catalog timeline, likely already applied.  
Moved to: archive/one_off_scripts/fix_catalog_timeline.py

**[2025-11-23] Archived: force_cluster_detection.py**  
Reason: Debug/test script for cluster detection, one-time use.  
Moved to: archive/one_off_scripts/force_cluster_detection.py

**[2025-11-23] Archived: QUICK_RUN_PHASE2.py**  
Reason: Quick test script, one-time use.  
Moved to: archive/one_off_scripts/QUICK_RUN_PHASE2.py

### A.2 Skrypty AWS emergency/fix

**[2025-11-23] Archived: aws_start_missing_9.sh**  
Reason: One-time fix for specific problem (missing run 9).  
Moved to: archive/one_off_scripts/aws_start_missing_9.sh

**[2025-11-23] Archived: CHECK_AWS_RESULTS.sh**  
Reason: One-time check script for AWS results, likely replaced by aws_test/scripts/check_*.  
Moved to: archive/one_off_scripts/CHECK_AWS_RESULTS.sh

**[2025-11-23] Archived: copy_fix_to_aws.ps1**  
Reason: One-time fix script for copying to AWS.  
Moved to: archive/one_off_scripts/copy_fix_to_aws.ps1

**[2025-11-23] Archived: copy_to_aws.ps1**  
Reason: One-time copy script, likely replaced by better solution.  
Moved to: archive/one_off_scripts/copy_to_aws.ps1

### A.3 Skrypty benchmark/test

**[2025-11-23] Archived: analyze_benchmark_results.ps1**  
Reason: Benchmark results analysis script, one-time use.  
Moved to: archive/one_off_scripts/analyze_benchmark_results.ps1

**[2025-11-23] Archived: run_benchmark.ps1**  
Reason: Benchmark script, one-time tests.  
Moved to: archive/one_off_scripts/run_benchmark.ps1

**[2025-11-23] Archived: run_cpu_benchmark.ps1**  
Reason: CPU benchmark script, one-time tests.  
Moved to: archive/one_off_scripts/run_cpu_benchmark.ps1

**[2025-11-23] Archived: run_hybrid_test.ps1**  
Reason: Hybrid mode test script, one-time use.  
Moved to: archive/one_off_scripts/run_hybrid_test.ps1

### A.4 Skrypty fix/cleanup

**[2025-11-23] Archived: cleanup_processes.ps1**  
Reason: Process cleanup script, one-time use.  
Moved to: archive/one_off_scripts/cleanup_processes.ps1

**[2025-11-23] Archived: fix_taichi_version.ps1**  
Reason: Taichi version fix script, one-time use (if already applied).  
Moved to: archive/one_off_scripts/fix_taichi_version.ps1

**[2025-11-23] Archived: fix_taichi_version.sh**  
Reason: Duplicate of fix_taichi_version.ps1, bash version.  
Moved to: archive/one_off_scripts/fix_taichi_version.sh

---

---

## [2025-11-23] Przeniesienie plików .md do docs/ (MD_FILES_MIGRATION_PLAN.md)

### B.1 Dokumentacja AWS → docs/aws_test/

**[2025-11-23] Moved: aws_test/scripts/CPU_THREADS_GUIDE.md**  
Reason: AWS documentation should be in docs/aws_test/ according to NAVIGATION_GUIDE.  
Moved to: docs/aws_test/CPU_THREADS_GUIDE.md

**[2025-11-23] Moved: aws_test/scripts/MONITORING_SCRIPTS_COMPARISON.md**  
Reason: AWS documentation should be in docs/aws_test/ according to NAVIGATION_GUIDE.  
Moved to: docs/aws_test/MONITORING_SCRIPTS_COMPARISON.md

### B.2 Dokumentacja techniczna → docs/technical/

**[2025-11-23] Moved: backend/sim/io/schema.md**  
Reason: Technical documentation should be in docs/technical/ according to NAVIGATION_GUIDE.  
Moved to: docs/technical/backend_sim_io_schema.md

**[2025-11-23] Moved: backend/tests/tests.md**  
Reason: Technical documentation should be in docs/technical/ according to NAVIGATION_GUIDE.  
Moved to: docs/technical/backend_tests.md

**[2025-11-23] Moved: backend/tests/TEST_SUMMARY.md**  
Reason: Technical documentation should be in docs/technical/ according to NAVIGATION_GUIDE.  
Moved to: docs/technical/backend_tests_summary.md

### B.4 Analizy wyników → docs/analysis/

**[2025-11-23] Moved: analysis/phase2b_miller_urey/ANALYSIS_FINDINGS.md**  
Reason: Analysis results should be in docs/analysis/ according to NAVIGATION_GUIDE.  
Moved to: docs/analysis/phase2b_miller_urey_ANALYSIS_FINDINGS.md

**[2025-11-23] Moved: analysis/phase2b_miller_urey/PAPER_SUMMARY.md**  
Reason: Analysis results should be in docs/analysis/ according to NAVIGATION_GUIDE.  
Moved to: docs/analysis/phase2b_miller_urey_PAPER_SUMMARY.md

### B.3, B.5 Dokumentacja CI/CD i scripts → docs/

**[2025-11-23] Moved: .github/CI_CHEATSHEET.md**  
Reason: CI/CD documentation moved to docs/ for better organization.  
Moved to: docs/CI_CHEATSHEET.md

**[2025-11-23] Moved: scripts/ANALYSIS_QUICK_REF.md**  
Reason: Scripts documentation moved to docs/ for better organization.  
Moved to: docs/scripts_ANALYSIS_QUICK_REF.md

### C.1 Raporty wyników → archive/old_docs/

**[2025-11-23] Archived: results/FINAL_ETA_REPORT.md**  
Reason: Results report, likely replaced by newer reports in docs/.  
Moved to: archive/old_docs/results_FINAL_ETA_REPORT.md

**[2025-11-23] Archived: results/FINAL_OPTIMIZATION_REPORT.md**  
Reason: Results report, likely replaced by newer reports in docs/.  
Moved to: archive/old_docs/results_FINAL_OPTIMIZATION_REPORT.md

**[2025-11-23] Archived: results/OPTIMIZATION_SUMMARY.md**  
Reason: Results report, likely replaced by newer reports in docs/.  
Moved to: archive/old_docs/results_OPTIMIZATION_SUMMARY.md

**[2025-11-23] Archived: results/spatial_hash_test/PERFORMANCE_REPORT.md**  
Reason: Results report, likely replaced by newer reports in docs/.  
Moved to: archive/old_docs/results_spatial_hash_PERFORMANCE_REPORT.md

**[2025-11-23] Archived: results/test_formamide_10k/TEST_REPORT.md**  
Reason: Results report, likely replaced by newer reports in docs/.  
Moved to: archive/old_docs/results_test_formamide_10k_TEST_REPORT.md

---

## [2025-11-23] Archiwizacja duplikatów (DUPLICATES_CANONICAL_PLAN.md)

### 1. Duplikaty check_status*.ps1

**[2025-11-23] Archived: backend/tests/check_status.ps1**  
Reason: Older version (2025-09-27), replaced by tests/check_status.ps1 (2025-10-06).  
Moved to: archive/deprecated/backend_tests_check_status.ps1

**[2025-11-23] Archived: backend/tests/check_status2.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/backend_tests_check_status2.ps1

**[2025-11-23] Archived: backend/tests/check_status3.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/backend_tests_check_status3.ps1

**[2025-11-23] Archived: backend/tests/check_status4.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/backend_tests_check_status4.ps1

**[2025-11-23] Archived: backend/tests/check_status5.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/backend_tests_check_status5.ps1

**[2025-11-23] Archived: backend/tests/check_status6.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/backend_tests_check_status6.ps1

**[2025-11-23] Archived: backend/tests/check_status7.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/backend_tests_check_status7.ps1

**[2025-11-23] Archived: backend/tests/check_status8.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/backend_tests_check_status8.ps1

**[2025-11-23] Archived: tests/check_status2.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/tests_check_status2.ps1

**[2025-11-23] Archived: tests/check_status3.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/tests_check_status3.ps1

**[2025-11-23] Archived: tests/check_status4.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/tests_check_status4.ps1

**[2025-11-23] Archived: tests/check_status5.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/tests_check_status5.ps1

**[2025-11-23] Archived: tests/check_status6.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/tests_check_status6.ps1

**[2025-11-23] Archived: tests/check_status7.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/tests_check_status7.ps1

**[2025-11-23] Archived: tests/check_status8.ps1**  
Reason: Iterative development version, check_status.ps1 is canonical.  
Moved to: archive/one_off_scripts/tests_check_status8.ps1

**Canonical version:** `tests/check_status.ps1` (remains in place)

### 3. Duplikaty start_backend*.ps1

**[2025-11-23] Archived: start_backend_simple.ps1**  
Reason: "Simple" version, replaced by start_backend.ps1 (canonical production version).  
Moved to: archive/deprecated/start_backend_simple.ps1

**Canonical version:** `start_backend.ps1` (remains in place)

---

## Podsumowanie

**Data:** 2025-11-23  
**Liczba zarchiwizowanych plików:** 17 (sekcja A.1-A.4) + 15 (duplikaty) + 5 (raporty) = 37 plików  
**Liczba przeniesionych plików do docs/:** 9 plików  
**Kategorie:**
- archive/one_off_scripts/: 17 + 14 = 31 plików
- archive/deprecated/: 2 pliki
- archive/old_docs/: 5 plików
- docs/aws_test/: 2 pliki
- docs/technical/: 3 pliki
- docs/analysis/: 2 pliki
- docs/: 2 pliki

**Plany źródłowe:**
- docs/ARCHIVE_PLAN_FINAL.md (sekcja A.1-A.4)
- docs/MD_FILES_MIGRATION_PLAN.md
- docs/DUPLICATES_CANONICAL_PLAN.md

**Operacje wykonane przez:** git mv (zachowanie historii Git)

---

## [2025-11-23] Organizacja dokumentów w docs/ (DOCS_ORGANIZATION_PLAN.md)

### Utworzone nowe podkatalogi:
- `docs/plans/` - plany i roadmapy
- `docs/optimization/` - optymalizacje i performance
- `docs/fixes/` - fixy i rozwiązania problemów
- `docs/archive/` - plany archiwizacji (wykonane)

### Przeniesienia do istniejących podkatalogów:

**docs/sessions/ (11 plików):**
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

**docs/aws_test/ (7 plików):**
- AWS_RESULTS_PIPELINE.md
- AWS_RUN_STATUS.md
- CHECK_AWS_STATUS.md
- DEBUG_AWS_PHASE2B.md
- FINAL_AWS_INSTRUCTIONS.md
- FIX_AWS_SCRIPTS.md
- INSTRUKCJA_AWS_UPD.md

**docs/troubleshooting/ (10 plików):**
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

**docs/technical/ (6 plików):**
- PHYSICS_DATABASE.md
- PHYSICS_DB_INTEGRATION.md
- THERMODYNAMIC_VALIDATION.md
- SCIENTIFIC_INTEGRITY_VERIFICATION.md
- SKAD_METRICS_BIERZE_CLUSTERS.md
- SPATIAL_HASHING_SUCCESS.md

**docs/analysis/ (6 plików):**
- ANALIZA_WYNIKOW_TESTU.md
- CHECK_REAL_RESULTS.md
- PERFORMANCE_ANALYSIS.md
- PERFORMANCE_DIAGNOSIS_FINAL.md
- PHASE2_RESULTS_ASSESSMENT.md
- TEST_RESULTS_OCT13.md

**docs/guides/ (10 plików):**
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

### Przeniesienia do nowych podkatalogów:

**docs/plans/ (9 plików):**
- AGGRESSIVE_OPTIMIZATION_PLAN.md
- BOND_ENHANCEMENT_PLAN.md
- CLUSTER_ENHANCEMENT_PLAN.md
- PHASE2_NEXT_STEPS_PL.md
- Live2_v1_plan.md
- live2-roadmap.md
- ROADMAP_COMPARISON_2025-10-13.md
- ROADMAP_UPDATE_2025-10-13.md
- Live 2-plan walidacji naukowej.md

**docs/optimization/ (12 plików):**
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

**docs/fixes/ (2 pliki):**
- ROZWIAZANIE_KROK_PO_KROKU.md
- PODSUMOWANIE_PROBLEM_KLASTROW.md

**docs/archive/ (7 plików):**
- ARCHIVE_PLAN.md
- ARCHIVE_PLAN_FINAL.md
- DUPLICATES_CANONICAL_PLAN.md
- MD_FILES_MIGRATION_PLAN.md
- REPO_STRUCTURE_ANALYSIS.md
- STRUCTURE_DEEP_ANALYSIS.md
- DOCS_ORGANIZATION_PLAN.md

**Podsumowanie organizacji docs/:**
- **Utworzone podkatalogi:** 4 nowe
- **Przeniesione pliki:** 61 plików
- **Kategorie:** sessions (11), aws_test (7), troubleshooting (10), technical (6), analysis (6), guides (10), plans (9), optimization (12), fixes (2), archive (7)

---

*Zgodnie z `docs/ARCHIVE_POLICY.md` sekcja 5, krok 4*

