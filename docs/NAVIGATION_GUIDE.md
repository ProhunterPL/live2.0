---
date: 2025-11-23
label: guide
---

# Live 2.0 – Navigation Guide  

**Centralny przewodnik po repozytorium | Phase 2B / publikacja**

Ten dokument prowadzi badaczy, developerów i agentów AI przez **strukturę Live 2.0**, wskazując:

- gdzie jest faktyczna logika symulacji,
- gdzie są skrypty do Phase 2B,
- gdzie znajdują się wyniki + jak z nich korzystać,
- które foldery są **święte** (read-only),
- gdzie trzymać eksperymenty, stare wersje i jednorazowe pliki,
- jak szukać rzeczy według tematu (chemia, fizyka, pipeline'y, wyniki, AWS).

To jest oficjalny "Single Source of Truth" do nawigacji po repo.

---

# 📁 1. HIGH-LEVEL OVERVIEW

Live 2.0 składa się z czterech głównych warstw:

## 1.1 Simulation Engine (CORE)  
**`backend/sim/**`**  
→ wykorzystywana przez Phase 2B, publikację i wszystkie eksperymenty naukowe.

**Kluczowe pliki:**
- `backend/sim/core/stepper.py` - główna pętla symulacji (1642 linie)
- `backend/sim/core/binding.py` - formowanie/zerwanie wiązań (Taichi kernels)
- `backend/sim/core/particles.py` - stan cząstek i fizyka
- `backend/sim/core/grid.py` - spatial hashing dla wydajności
- `backend/sim/core/catalog.py` - wykrywanie i śledzenie molekuł
- `backend/sim/chemistry/physics_db.py` - tablica okresowa, energie wiązań
- `backend/sim/chemistry/reactions.py` - wykrywanie i walidacja reakcji
- `backend/sim/config.py` - konfiguracja symulacji (SimulationConfig)
- `backend/sim/molecule_extractor.py` - post-processing: ekstrakcja molekuł ze snapshotów

## 1.2 Phase 2B Orchestration  
**`scripts/run_phase2_full.py`, `aws_test/*`, skrypty AWS queue**  
→ wykonywanie 500K-step symulacji, snapshoty, checkpointy, restart queue.

**Kluczowe pliki:**
- `scripts/run_phase2_full.py` - główny runner dla Phase 2 (Phase2FullRunner)
- `aws_test/configs/*_SUPER_FAST.yaml` - konfiguracje produkcyjne
- `aws_test/configs/phase2_*.yaml` - definicje scenariuszy (miller_urey, hydrothermal, formamide)
- `aws_test/scripts/auto_queue_restart.sh` - system auto-restart (max 4 równoległe)
- `aws_test/scripts/check_actual_progress.py` - monitorowanie postępu
- `aws_test/scripts/kill_stuck_phase2b.sh` - awaryjne czyszczenie

## 1.3 Data & Results  
**`results/**`, `phase2b_*_results/**`, `run_*/snapshots/**`**  
→ surowe dane, wyniki, JSON-y, struktury plików generowane automatycznie.

**Struktura wyników:**
```
results/
└── phase2b_additional/
    └── miller_urey_extended/
        └── run_X/              # Każdy run zawiera:
            ├── results.json     # Metryki symulacji
            ├── molecules.json   # Wykryte molekuły
            ├── snapshots/       # Snapshoty co 50K kroków
            └── checkpoints/     # Checkpointy co 100K kroków (markery)
```

## 1.4 Documentation  
**`docs/**`**  
→ cała wiedza naukowa, techniczna, analizy, raporty, notatki.

**Struktura dokumentacji:**
```
docs/
├── INDEX.md                          ⭐ START HERE - Główny indeks
├── NAVIGATION_GUIDE.md              📍 Ten plik
├── ARCHIVE_POLICY.md                📋 Polityka archiwizacji
│
├── 📅 sessions/                      Sesje robocze (chronologicznie)
│   ├── SESSION_SUMMARY_OCT13.md
│   ├── SESSION_SUMMARY_2025-10-13.md
│   └── [inne podsumowania sesji]
│
├── 🔧 technical/                     Dokumentacja techniczna
│   ├── PHYSICS_DATABASE.md
│   ├── THERMODYNAMIC_VALIDATION.md
│   └── [inne dokumenty techniczne]
│
├── 📖 guides/                        Przewodniki użytkownika
│   ├── QUICK_START.md
│   ├── ENVIRONMENT_SETUP.md
│   └── [inne przewodniki]
│
├── 🐛 troubleshooting/                Problemy i rozwiązania
│   ├── CLUSTER_DETECTION_ISSUE.md
│   ├── GPU_MEMORY_ISSUE.md
│   └── [inne fixy]
│
├── 📊 analysis/                      Analizy wyników
│   ├── PERFORMANCE_ANALYSIS.md
│   ├── PHASE2_RESULTS_ASSESSMENT.md
│   └── [inne analizy]
│
├── ☁️ aws_test/                      Dokumentacja AWS
│   ├── AUTO_RESTART_GUIDE.md
│   ├── CLUSTER_DEADLOCK_FIX.md
│   └── [inne dokumenty AWS]
│
├── 📋 plans/                         Plany i roadmapy
│   ├── live2-roadmap.md
│   ├── AGGRESSIVE_OPTIMIZATION_PLAN.md
│   └── [inne plany]
│
├── ⚡ optimization/                  Optymalizacje i performance
│   ├── FPS_OPTIMIZATION.md
│   ├── HYBRID_GPU_CPU_GUIDE.md
│   └── [inne optymalizacje]
│
├── 🔧 fixes/                         Fixy i rozwiązania problemów
│   ├── ROZWIAZANIE_KROK_PO_KROKU.md
│   └── PODSUMOWANIE_PROBLEM_KLASTROW.md
│
├── 📦 archive/                       Plany archiwizacji (wykonane)
│   ├── ARCHIVE_PLAN_FINAL.md
│   ├── DUPLICATES_CANONICAL_PLAN.md
│   └── [inne plany archiwizacji]
│
├── phase2b/                          Dokumentacja Phase 2B
│   └── PHASE2B_STATUS.md
│
└── 🔬 [inne katalogi]                Pozostała dokumentacja
    ├── SCIENTIFIC_OVERVIEW.md
    ├── VALIDATION_ROADMAP.md
    └── ...
```

---

# 🟥 2. READ-ONLY ZONES | (Do not modify during Phase 2B)

Te katalogi są **krytyczne naukowo** i ich zmiana może:
- unieważnić Phase 2B,
- zniszczyć możliwość odtworzenia wyników,
- uszkodzić pipeline publikacji.

**Dlatego traktujemy je jako read-only:**

## 2.1 🔬 CORE (fizyka, chemia, główny symulator)

**`backend/sim/core/**`**  
- `stepper.py` - główna pętla symulacji
- `binding.py` - formowanie/zerwanie wiązań (Taichi kernels)
- `particles.py` - stan cząstek i fizyka
- `grid.py` - spatial hashing
- `catalog.py` - wykrywanie molekuł

**`backend/sim/chemistry/**`**  
- `physics_db.py` - tablica okresowa, energie wiązań
- `reactions.py` - wykrywanie i walidacja reakcji

**`backend/sim/config.py`**  
- `SimulationConfig` - wszystkie parametry symulacji

**`backend/sim/molecule_extractor.py`**  
- Post-processing: ekstrakcja molekuł ze snapshotów

## 2.2 📊 Phase 2B Orchestration & Analysis

**`scripts/run_phase2_full.py`**  
- `Phase2FullRunner` - orkiestracja pełnych symulacji

**`scripts/phase2.py`**  
- Dodatkowe skrypty Phase 2

**`scripts/*phase2*.py`**  
- Wszystkie skrypty związane z Phase 2

**`aws_test/scripts/**`**  
- Queue, monitoring, kill scripts
- `auto_queue_restart.sh` - system auto-restart
- `check_actual_progress.py` - monitorowanie postępu
- `kill_stuck_phase2b.sh` - awaryjne czyszczenie

**`aws_test/configs/phase2_*.yaml`**  
- Definicje scenariuszy (miller_urey, hydrothermal, formamide)

**`aws_test/configs/*SUPER_FAST*.yaml`**  
- Konfiguracje produkcyjne (muszą mieć `cluster_check_interval: 999999999`)

## 2.3 📚 Scientific Documentation

**`docs/phase2b/**`**  
- Dokumentacja Phase 2B

**`docs/technical/**`**  
- Dokumentacja techniczna

**`docs/VALIDATION_ROADMAP.md`**  
- Roadmap walidacji

**`docs/SCIENTIFIC_VALIDATION_ANALYSIS.md`**  
- Analiza walidacji naukowej

## 2.4 📦 Completed Results (Raw Data)

**`results/**`**  
- Wszystkie wyniki symulacji

**`phase2b_*_results/**`**  
- Wyniki Phase 2B

**Wszystkie katalogi `run_*/` zawierające:**
- `results.json`
- `molecules.json`
- `snapshots/` - snapshoty co 50K kroków
- `checkpoints/` - checkpointy co 100K kroków (markery)

**`run_*/snapshots/**`**  
- Wszystkie snapshoty z ukończonych biegów

**`run_*/checkpoints/**`**  
- Wszystkie checkpointy z ukończonych biegów

### ⚠️ Zasady dla Read-Only Zones:

1. **Nie usuwać, nie zmieniać nazw, nie przenosić** plików bez wyraźnej zgody użytkownika
2. **Przy zmianach kodu:**
   - Zawsze pokazać jasny `PLAN` najpierw
   - Wyjaśnić wpływ naukowy i zaproponować testy/backtesty
3. **Nigdy nie wykonywać** masowego "czyszczenia" lub refaktoryzacji w tych ścieżkach jako część automatycznej organizacji

---

# 🗂️ 3. ARCHIVE POLICY | Struktura Archiwum

Live 2.0 to długotrwały projekt naukowy. **Nigdy nie usuwamy na stałe** rzeczy podczas czyszczenia. Zamiast tego używamy **strukturyzowanego archiwum**.

## 3.1 Struktura Archiwum

Wszystkie pliki i foldery jednorazowe, przestarzałe lub zduplikowane powinny być przeniesione (przez `git mv`) do:

- **`archive/one_off_scripts/`** – skrypty jednorazowe / debugujące
- **`archive/old_docs/`** – przestarzała dokumentacja, zastąpiona przez nowsze w `docs/`
- **`archive/experiments/`** – eksperymentalne pipeline'y, prototypy, kod scratch
- **`archive/deprecated/`** – kod/konfiguracje zastąpione finalną wersją
- **`archive/tmp_results/`** – tymczasowe foldery wyników, które nie są już potrzebne do analizy
- **`archive/ARCHIVE_LOG.md`** – log wszystkich operacji archiwizacji

**Struktura archive/:**
```
archive/
├── ARCHIVE_LOG.md              # Log wszystkich operacji archiwizacji
├── one_off_scripts/            # Skrypty jednorazowe (32 pliki)
├── old_docs/                   # Stara dokumentacja (6 plików)
├── experiments/                # Eksperymenty i prototypy
├── deprecated/                 # Zastąpione wersje (3 pliki)
└── tmp_results/                # Tymczasowe wyniki
```

## 3.2 Zasady Czyszczenia

Gdy proszony o "czyszczenie", "organizację" lub "usunięcie duplikatów":

1. **NIGDY nie usuwać plików na stałe**
   - Użyj `git mv` aby przenieść je do `archive/...`
   - Zachowaj oryginalną względną strukturę gdy możliwe (np. `scripts/foo.py` → `archive/one_off_scripts/scripts/foo.py`)

2. **ZAWSZE zaproponować plan najpierw:**
   - Wylistować wszystkie pliki, które uważasz za:
     - zduplikowane,
     - jednorazowe / debug,
     - do przeniesienia do archiwum
   - Poczekać na potwierdzenie użytkownika przed zastosowaniem planu

3. **Preferować oznaczanie jako legacy zamiast nadpisywania:**
   - Gdy znajdziesz dwa podobne skrypty:
     - Wybierz JEDEN jako wersję "kanoniczną"
     - Przenieś pozostałe do `archive/` i dodaj krótką notatkę w pliku kanonicznym

4. **NIE przenosić:**
   - niczego pod `backend/sim/core/**` lub `backend/sim/chemistry/**`,
   - niczego pod `docs/phase2b/**` i `docs/technical/**`,
   - żadnych ukończonych katalogów wyników (`run_*/` z `results.json`)

---

# 📘 3. WHERE TO FIND THE MOST IMPORTANT THINGS

## 3.1 Core Simulation Code

> **Jeśli chcesz wiedzieć, jak działa symulator → zacznij tutaj**

| Obszar | Plik / folder |
|--------|----------------|
| Główny loop | `backend/sim/core/stepper.py` |
| Siatka (spatial hashing) | `backend/sim/core/grid.py` |
| Fizyka cząstek | `backend/sim/core/particles.py` |
| Tworzenie/zerwanie wiązań | `backend/sim/core/binding.py` |
| Wykrywanie molekuł | `backend/sim/core/catalog.py` |
| Konfiguracja symulacji | `backend/sim/config.py` |

## 3.2 Phase 2B Simulations (500K steps)

To jest kluczowy pipeline do publikacji:

| Funkcja | Lokalizacja |
|---------|-------------|
| Uruchamianie pełnej symulacji | `scripts/run_phase2_full.py` |
| Kolejkowanie, restartowanie AWS | `aws_test/scripts/auto_queue_restart.sh` |
| Monitoring | `aws_test/scripts/check_actual_progress.py` |
| Analiza wyników Phase 2B | `scripts/analyze_phase2b_complete.py` |
| Post-processing molekuł | `backend/sim/molecule_extractor.py` |

## 3.3 AWS Configs (najważniejsze parametry)

**`aws_test/configs/*SUPER_FAST.yaml`**  
**`aws_test/configs/phase2_*.yaml`**

**Parametr, który musi być poprawny:**
```yaml
cluster_check_interval: 999999999
```

## 3.4 Results & Snapshots

**Struktura wyników:**
```
results/<scenario>/run_X/
├── results.json      # Metryki symulacji
├── molecules.json    # Wykryte molekuły
├── snapshots/        # Snapshoty co 50K kroków
└── checkpoints/      # Checkpointy co 100K kroków (markery)
```

**Uwaga:** snapshoty i checkpointy mogą mieć dziesiątki tysięcy plików → nie otwieraj w AI bez potrzeby.

## 3.5 Documentation Structure

| Cel dokumentacji | Folder |
|------------------|--------|
| AWS deployment | `docs/aws_test/` |
| Lokalne uruchamianie | `docs/local/` |
| Phase 2B | `docs/phase2b/` |
| Technical (chemia, fizyka, architektura) | `docs/technical/` |
| Walidacja naukowa | `docs/VALIDATION_ROADMAP.md` |
| Sesje GPT/notatki | `docs/sessions/` |
| Analizy wyników | `docs/analysis/` |
| Przewodniki użytkownika | `docs/guides/` |
| Problemy i fixy | `docs/troubleshooting/` |
| Plany i roadmapy | `docs/plans/` |
| Optymalizacje | `docs/optimization/` |
| Fixy i rozwiązania | `docs/fixes/` |
| Plany archiwizacji (wykonane) | `docs/archive/` |

---

# 🔍 4. JAK ZNALEŹĆ TO CZEGO SZUKASZ?

## 4.1 Według Tematu

### 🔬 Fizyczne jądro symulacji
- `backend/sim/core/particles.py` - fizyka cząstek
- `backend/sim/core/stepper.py` - główna pętla symulacji

### ⚛️ Wiązania i chemia
- `backend/sim/chemistry/physics_db.py` - reguły wiązań, tablica okresowa
- `backend/sim/chemistry/reactions.py` - reakcje chemiczne
- `backend/sim/core/binding.py` - formowanie/zerwanie wiązań

### 🧪 Molekuły i detekcja
- `backend/sim/core/catalog.py` - wykrywanie molekuł
- `backend/sim/molecule_extractor.py` - ekstrakcja molekuł ze snapshotów

### 📈 Analizy i narzędzia
- `scripts/analyze_phase2b_complete.py` - analiza Phase 2B
- `scripts/analyze_results.py` - analiza wyników

### 🏭 AWS Execution
- `aws_test/scripts/**` - skrypty AWS (queue, monitoring)
- `aws_test/configs/**` - konfiguracje AWS

### 🧪 Chemia i Fizyka (szczegóły)
- **Reguły wiązań**: `backend/sim/chemistry/physics_db.py`
- **Reakcje**: `backend/sim/chemistry/reactions.py`
- **Fizyka cząstek**: `backend/sim/core/particles.py`
- **Formowanie wiązań**: `backend/sim/core/binding.py`

### ⚙️ Pipeline Symulacji
- **Główna pętla**: `backend/sim/core/stepper.py` → metoda `step()`
- **Runner Phase 2**: `scripts/run_phase2_full.py` → klasa `Phase2FullRunner`
- **Konfiguracja**: `backend/sim/config.py` → klasa `SimulationConfig`

### 📊 Wyniki i Analiza
- **Wyniki**: `results/phase2b_additional/`
- **Ekstrakcja molekuł**: `backend/sim/molecule_extractor.py`
- **Analiza**: `scripts/analyze_results.py`

### ☁️ AWS Deployment
- **Konfiguracje**: `aws_test/configs/*_SUPER_FAST.yaml`
- **Queue system**: `aws_test/scripts/auto_queue_restart.sh`
- **Monitoring**: `aws_test/scripts/check_actual_progress.py`
- **Dokumentacja AWS**: `docs/technical/aws/`

### 📚 Dokumentacja
- **Główny indeks**: `docs/INDEX.md`
- **Quick start**: `docs/guides/QUICK_START.md`
- **Parametry**: `docs/technical/parameters/`
- **Matcher**: `docs/technical/matcher/`
- **Status Phase 2B**: `docs/phase2b/PHASE2B_STATUS.md`

## 4.2 Według Akcji

### 📱 "Chcę szybko zacząć"
→ `docs/guides/QUICK_START.md`

### 🔍 "Szukam konkretnego dokumentu"
→ `docs/INDEX.md` - Pełny indeks z opisami

### 📅 "Co było zmienione ostatnio?"
→ `docs/sessions/` - Najnowsza sesja

### 🔧 "Chcę zmienić parametry"
→ `docs/technical/parameters/`

### 🧪 "Jak działa matcher?"
→ `docs/technical/matcher/`

### ☁️ "Jak wdrożyć na AWS?"
→ `docs/technical/aws/`

### 🐛 "Mam problem"
→ `docs/CRASH_REPORT.md`  
→ `docs/PERFORMANCE_DIAGNOSIS_FINAL.md`

### 📊 "Jakie są plany rozwoju?"
→ `docs/live2-roadmap.md`  
→ `docs/VALIDATION_ROADMAP.md`

---

# ⚠️ 5. ZNANE PROBLEMY & ROZWIĄZANIA

## 5.1 Cluster Detection Deadlock (NAPRAWIONE)
**Problem**: O(N²) cluster detection powoduje nieskończoną pętlę w złożonych sieciach  
**Lokalizacja**: `backend/sim/core/stepper.py` linia 510  
**Naprawa**: Czyta `cluster_check_interval` z konfiguracji zamiast hardcoded 1200  
**Fix konfiguracji**: Ustaw `cluster_check_interval: 999999999` we wszystkich SUPER_FAST configs  
**Wpływ**: Chemia nadal dokładna, tylko metryki cluster wyłączone

## 5.2 Pusty Molecule Catalog (Post-Processing Fix)
**Problem**: Catalog nie aktualizowany podczas symulacji → pusty `molecules.json`  
**Lokalizacja**: `backend/sim/core/catalog.py`  
**Rozwiązanie**: Użyj `molecule_extractor.py` aby ekstrahować ze snapshotów post-symulacja  
**Skrypt**: `scripts/fix_run1_molecules.py`

## 5.3 Checkpoints Nie Wznawiają
**Problem**: Checkpointy to tylko markery, nie pełny stan  
**Lokalizacja**: `scripts/run_phase2_full.py` linia 434  
**Workaround**: Restart od początku (system queue to obsługuje)

---

# 🚀 6. COMMON OPERATIONS

## 6.1 Uruchomienie Pojedynczej Symulacji (Lokalnie)
```bash
python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \
    --output results/test_run \
    --seed 100 \
    --steps 500000 \
    --force-cpu
```

## 6.2 Uruchomienie Phase 2B na AWS (Wiele Równoległych)
```bash
# Start auto-restart queue system (4 parallel max)
cd ~/live2.0
nohup bash aws_test/scripts/auto_queue_restart.sh > logs/auto_restart_main.log 2>&1 &

# Monitor progress
python3 aws_test/scripts/check_actual_progress.py
```

## 6.3 Ekstrakcja Molekuł z Ukończonego Run
```python
from backend.sim.molecule_extractor import extract_molecules_from_results
results = extract_molecules_from_results("results/phase2b_additional/miller_urey_extended/run_1")
```

## 6.4 Analiza Wyników
```bash
python3 aws_test/scripts/analyze_additional_results.py
```

---

# 📋 7. KONWENCJE NAZEWNICTWA

## 7.1 Typy Plików Dokumentacji
- `README.md` - Wprowadzenie do katalogu
- `*_GUIDE.md` - Przewodnik użytkownika
- `*_FIX.md` - Naprawa/poprawka
- `*_ANALYSIS.md` - Analiza techniczna
- `SESSION_*.md` - Podsumowanie sesji
- `PHASE*_*.md` - Dokumentacja fazy projektu

## 7.2 Priorytety Dokumentacji
- ⭐ **START HERE** - Zacznij tutaj
- 📍 **Important** - Ważne dokumenty
- 🔍 **Reference** - Dokumenty referencyjne

## 7.3 Konfiguracje
- `phase2_<scenario>_<mode>.yaml` (np. `phase2_miller_urey_extended_SUPER_FAST.yaml`)

## 7.4 Wyniki
- `run_<number>/` (np. `run_1/`, `run_2/`)

---

# 📂 8. WHERE TO PUT NEW CODE (BEST PRACTICES)

## 8.1 Nowe skrypty (analiza, pomocnicze)
→ **`scripts/`**

## 8.2 Nowe dokumenty (przepisy, analizy, status)
→ **`docs/<odpowiednia-kategoria>/`**

**Dostępne kategorie:**
- `docs/sessions/` - podsumowania sesji
- `docs/aws_test/` - dokumentacja AWS
- `docs/technical/` - dokumentacja techniczna
- `docs/analysis/` - analizy wyników
- `docs/guides/` - przewodniki użytkownika
- `docs/troubleshooting/` - problemy i fixy
- `docs/plans/` - plany i roadmapy
- `docs/optimization/` - optymalizacje i performance
- `docs/fixes/` - fixy i rozwiązania problemów
- `docs/phase2b/` - dokumentacja Phase 2B
- `docs/local/` - dokumentacja lokalna

**Zasada:**  
**Nigdy nie tworzymy dokumentacji poza `docs/`.**

## 8.3 Eksperymenty / prototypy (nie-do-końca pewne)
→ **`archive/experiments/`**

## 8.4 Stare wersje / deprecated
→ **`archive/old_docs/`**, **`archive/deprecated/`**, **`archive/one_off_scripts/`**

---

# 🧹 9. CLEANUP RULES (SAFE MODE)

1. **Nic nie usuwamy** podczas Phase 2B → tylko przenosimy.

2. Używamy **`git mv`** aby zachować historię.

3. Przed sprzątaniem agent musi pokazać plan:
   - lista plików,
   - lista folderów,
   - ocena duplikatów,
   - propozycje przeniesienia.

4. Zasada: **CORE / PHASE 2B / RESULTS** = zawsze read-only.

5. Wszystkie "jednorazówki", debug-scripts, playgroundy → do `archive/`.

---

# 🧠 10. AI INTERACTION RULES (For Cursor & LLMs)

## 10.1 Kiedy agent może zmieniać kod:
- naprawa błędów,
- lokalny refactor jednego pliku,
- przygotowanie małych narzędzi analitycznych.

## 10.2 Kiedy AI *musi* pytać użytkownika:
- zmiany w:
  - `backend/sim/core/**`,
  - `backend/sim/chemistry/**`,
  - Phase 2B pipeline,
  - YAML SUPER_FAST configs,
  - dokumentacji naukowej.

## 10.3 Kiedy AI nie powinno czytać danych:
- >2 MB wyników lub snapshotów,
- `results/**`,
- `phase2b_*_results/**`.

---

# 🧪 11. PHASE 2B PIPELINE (SHORT SUMMARY)

1. **Start run:**  
   `scripts/run_phase2_full.py`

2. **AWS monitoring:**  
   `aws_test/scripts/check_actual_progress.py`

3. **Queue management:**  
   `aws_test/scripts/auto_queue_restart.sh`

4. **Outputs:**  
   `results/<scenario>/run_X/`

5. **Post-processing:**  
   `molecule_extractor.py` → `molecules.json`

6. **Final analysis for publication:**  
   `scripts/analyze_phase2b_complete.py`  
   `docs/phase2b/**`

---

# 🧹 12. LOW-PRIORITY / NO-CONTEXT PATHS

Podczas wyszukiwania kontekstu, **unikaj** lub traktuj jako bardzo niski priorytet:

- **Duże surowe dane i drzewa wyników:**
  - `results/**`
  - `phase2b_*_results/**`
  - wszystkie `run_*/snapshots/**`
  - wszystkie `run_*/checkpoints/**`

- **Logi i metryki:**
  - `**/*.log`
  - `logs/**`
  - `**/metrics_*.csv`
  - `diagnostics/**`

- **Cache i pliki tymczasowe:**
  - `**/__pycache__/**`
  - `.pytest_cache/**`
  - `.mypy_cache/**`
  - `.ruff_cache/**`

**Użyj kodu i dokumentacji jako głównego źródła prawdy:**
- `backend/sim/**`
- `aws_test/configs/**`
- `scripts/**`
- `docs/**`

Czytaj pliki danych tylko gdy:
- użytkownik wyraźnie prosi o analizę konkretnego run lub datasetu, lub
- potrzebujesz przykładu formatu rekordu.

---

# 📌 13. QUICK LOOKUP

**Znajdź: fizykę cząstek** → `backend/sim/core/particles.py`  

**Znajdź: wiązania** → `backend/sim/core/binding.py`  

**Znajdź: wykrywanie klastrów/molekuł** → `backend/sim/core/catalog.py`  

**Znajdź: konfiguracje** → `backend/sim/config.py`  

**Znajdź: pipeline Phase2B** → `scripts/run_phase2_full.py`  

**Znajdź: raporty naukowe** → `docs/phase2b/`  

**Znajdź: AWS** → `aws_test/scripts/`, `aws_test/configs/`

---

# 🔗 14. SZYBKIE LINKI

| Kategoria | Link | Opis |
|-----------|------|------|
| 🚀 Start | `docs/INDEX.md` | Główny punkt wejścia |
| 📅 Najnowsze | `docs/sessions/` | Ostatnie zmiany |
| 🔧 Parametry | `docs/technical/parameters/` | Parametry naukowe |
| 🧪 Matcher | `docs/technical/matcher/` | Identyfikacja molekuł |
| ☁️ AWS | `docs/technical/aws/` | Cloud deployment |
| 📖 Quick Start | `docs/guides/QUICK_START.md` | Szybki start |
| 🔬 Core | `backend/sim/core/stepper.py` | Główna pętla symulacji |
| 📊 Phase 2B | `scripts/run_phase2_full.py` | Runner Phase 2B |

---

# 💡 15. WSKAZÓWKI

## 9.1 Dla Nowych Użytkowników
1. Zacznij od `docs/INDEX.md`
2. Przeczytaj `docs/guides/QUICK_START.md`
3. Eksploruj `docs/sessions/` dla historii

## 9.2 Dla Powracających
1. Sprawdź najnowszą sesję w `docs/sessions/`
2. Zobacz "Historia Aktualizacji" w `docs/INDEX.md`
3. Przejrzyj zmiany w `docs/technical/`

## 9.3 Dla Rozwijających
1. Czytaj `docs/sessions/` dla kontekstu decyzji
2. Sprawdzaj `docs/technical/` dla szczegółów
3. Aktualizuj odpowiednie README przy zmianach
4. **Zawsze sprawdzaj read-only zones przed modyfikacjami**

---

# 🏁 16. FINAL NOTES

- Live 2.0 jest projektem **naukowym**, nie tylko programistycznym.

- Kluczem jest **zachowanie odtwarzalności eksperymentów**.

- Ta nawigacja + `.cursorrules` + `archive/` dają:
  - bezpieczne sprzątanie,
  - stabilny workflow dla AI,
  - 100% kontrolę nad tym, co jest current i co jest legacy.

---

*Masz pytania? Sprawdź `docs/INDEX.md` lub otwórz issue na GitHub.*

**Last Updated:** 2025-11-23  

**Maintainer:** Michał Klawikowski  

**Purpose:** Centralna mapa Live 2.0 dla ludzi i agentów AI
