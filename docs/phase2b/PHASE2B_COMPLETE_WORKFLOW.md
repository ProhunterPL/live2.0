# 🔄 Kompletny Workflow Analizy Phase 2B

## 📋 Przegląd Workflow

Po zakończeniu symulacji na AWS, wykonaj następujące kroki w kolejności:

### Krok 1: Napraw Molekuły (Jeśli Katalog Pusty)

**Problem**: Symulacje mogą mieć pusty katalog (molekuły nie były rejestrowane podczas symulacji).

**Rozwiązanie**: Wyekstraktuj molekuły z snapshotów.

```bash
# Napraw wszystkie runy w wszystkich scenariuszach
python scripts/fix_run1_molecules.py --all

# LUB napraw tylko konkretny scenariusz
python scripts/fix_run1_molecules.py --scenario miller_urey_extended

# LUB napraw tylko jeden run
python scripts/fix_run1_molecules.py --run results/phase2b_additional/miller_urey_extended/run_1
```

**Co to robi:**
- Ekstraktuje molekuły z snapshotów (używa bonds i clusters)
- Aktualizuje `results.json` z wyekstraktowanymi molekułami
- Aktualizuje `molecules.json`
- Działa dla wszystkich runów automatycznie

**Czas**: ~1-2 minuty na run (zależy od liczby snapshotów)

---

### Krok 2: Batch Analysis (Podstawowa Analiza)

**Skrypt**: `scripts/analyze_phase2_batch.py`

```bash
# Analizuj wszystkie runy rekurencyjnie
python scripts/analyze_phase2_batch.py \
    --input results/phase2b_additional \
    --output analysis/phase2b_complete \
    --recursive \
    --use-matcher
```

**Co to robi:**
- Znajduje wszystkie `results.json` w katalogu
- **Czyta molekuły z results.json** (które zostały naprawione w Kroku 1)
- Tworzy raporty analizy dla każdego runu
- Agreguje wyniki przez scenariusze
- Opcjonalnie używa MatcherV2 do identyfikacji PubChem

**UWAGA**: Ten skrypt **NIE** ekstraktuje z snapshotów - używa molekuł z `results.json`. 
**Zawsze najpierw uruchom Krok 1**, aby upewnić się, że `results.json` zawiera molekuły!

**Output:**
- `analysis/phase2b_complete/batch_results.json` - wszystkie wyniki
- `analysis/phase2b_complete/batch_report.txt` - raport tekstowy
- `analysis/phase2b_complete/{scenario}_summary.json` - podsumowania per scenariusz

---

### Krok 3: Kompletna Analiza (Autocatalysis + Complexity)

**Skrypt**: `scripts/analyze_phase2b_complete.py`

```bash
# Pełna analiza z autocatalysis i complexity metrics
python scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data
```

**Co to robi:**
- Analizuje każdy scenariusz (10 runów)
- Wykrywa cykle autocatalytic
- Oblicza metryki złożoności
- Porównuje scenariusze
- Generuje dane dla figurek
- Generuje tabele LaTeX
- Tworzy snippetki LaTeX dla publikacji

**Output:**
- `paper/results_data/{scenario}_analysis.json` - analiza per scenariusz
- `paper/results_data/scenario_comparison.json` - porównanie
- `paper/results_data/summary_table.csv` - tabela podsumowująca
- `paper/results_data/figure_data.json` - dane dla figurek
- `paper/results_data/latex_snippets.txt` - snippetki LaTeX

---

### Krok 4: Generowanie Figurek i Tabel

**Skrypty** (jeśli istnieją):
- `scripts/generate_all_figures.py` - generuje wszystkie figureki
- `scripts/generate_all_tables.py` - generuje wszystkie tabele

```bash
# Generuj figureki
python scripts/generate_all_figures.py \
    --data paper/results_data \
    --output paper/figures

# Generuj tabele
python scripts/generate_all_tables.py \
    --data paper/results_data \
    --output paper/tables
```

---

## 🔄 Kompletny Pipeline (Jedna Komenda)

Możesz stworzyć master script, który wykonuje wszystko:

```bash
# 1. Napraw molekuły
python scripts/fix_run1_molecules.py --all

# 2. Batch analysis
python scripts/analyze_phase2_batch.py \
    --input results/phase2b_additional \
    --output analysis/phase2b_complete \
    --recursive \
    --use-matcher

# 3. Kompletna analiza
python scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data

# 4. Generuj figureki i tabele (jeśli skrypty istnieją)
python scripts/generate_all_figures.py --data paper/results_data --output paper/figures
python scripts/generate_all_tables.py --data paper/results_data --output paper/tables
```

---

## ⚠️ Ważne Uwagi

### 1. Kolejność Jest Ważna

**Zawsze najpierw Krok 1** (napraw molekuły), potem Krok 2 i 3.

### 2. Sprawdź Czy Molekuły Są W results.json

Przed Krok 2, sprawdź czy `results.json` zawiera molekuły:

```bash
python -c "import json; d=json.load(open('results/phase2b_additional/miller_urey_extended/run_1/results.json')); print('Molecules:', len(d.get('molecules_detected', [])))"
```

Jeśli 0, uruchom Krok 1.

### 3. analyze_phase2_batch.py vs analyze_phase2b_complete.py

- **`analyze_phase2_batch.py`**: Podstawowa ekstrakcja i agregacja
- **`analyze_phase2b_complete.py`**: Zaawansowana analiza (autocatalysis, complexity)

**Użyj obu** - najpierw batch, potem complete.

---

## 📊 Oczekiwane Wyniki

Po wykonaniu wszystkich kroków powinieneś mieć:

1. ✅ **Wszystkie results.json z molekułami** (po Kroku 1)
2. ✅ **Batch analysis report** (po Kroku 2)
3. ✅ **Complete analysis data** (po Kroku 3)
4. ✅ **Figureki i tabele** (po Kroku 4)

---

## 🚀 Quick Start (Po Zakończeniu AWS)

```bash
# 1. Napraw wszystkie runy
python scripts/fix_run1_molecules.py --all

# 2. Sprawdź czy działa
python -c "import json; d=json.load(open('results/phase2b_additional/miller_urey_extended/run_1/results.json')); print('Molecules:', len(d.get('molecules_detected', [])))"

# 3. Batch analysis
python scripts/analyze_phase2_batch.py --input results/phase2b_additional --output analysis/phase2b_complete --recursive

# 4. Complete analysis
python scripts/analyze_phase2b_complete.py --input results/phase2b_additional --output paper/results_data
```

---

## 📝 Checklist

- [ ] Krok 1: Napraw molekuły (`fix_run1_molecules.py --all`)
- [ ] Weryfikacja: Sprawdź czy molekuły są w results.json
- [ ] Krok 2: Batch analysis (`analyze_phase2_batch.py`)
- [ ] Krok 3: Complete analysis (`analyze_phase2b_complete.py`)
- [ ] Krok 4: Generuj figureki i tabele (jeśli skrypty istnieją)
- [ ] Sprawdź outputy w `paper/results_data/`

---

**Workflow jest gotowy!** 🎉

