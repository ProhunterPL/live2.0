---
date: 2025-11-28
label: plan
---

# Truth-Filter: System Walidacji Wyników dla Publikacji

## 🎯 Cel

Truth-Filter to kompleksowy system walidacji wyników symulacji przed publikacją, zapewniający że tylko wiarygodne, naukowo poprawne wyniki trafiają do manuskryptu.

## 📋 Podsumowanie

Truth-Filter integruje istniejące komponenty walidacji (termodynamika, chemia, dopasowania PubChem) w jeden spójny pipeline, który:
- Filtruje molekuły pod kątem wiarygodności chemicznej
- Weryfikuje zgodność z literaturą
- Sprawdza jakość dopasowań PubChem
- Waliduje termodynamikę symulacji
- Generuje raporty walidacji dla każdego runa

**Status**: Projekt do implementacji

---

## 🏗️ Architektura

### Moduły Filtrowania

```
TruthFilter
├── MoleculeFilter          # Filtrowanie molekuł (real vs clusters)
├── ThermodynamicValidator  # Walidacja termodynamiki (już istnieje)
├── LiteratureValidator     # Zgodność z literaturą (benchmark_reactions)
├── MatchConfidenceFilter   # Wiarygodność dopasowań PubChem
├── SimulationQualityFilter # Jakość symulacji (completion, stability)
└── TruthReportGenerator    # Generowanie raportów
```

### Poziomy Walidacji

1. **STRICT** (dla publikacji)
   - Wszystkie filtry muszą przejść
   - Confidence > 0.8 dla dopasowań
   - Zgodność z literaturą ±20%
   - Termodynamika: drift < 0.5%

2. **MEDIUM** (dla analizy)
   - Większość filtrów musi przejść
   - Confidence > 0.6
   - Zgodność z literaturą ±30%
   - Termodynamika: drift < 1%

3. **LENIENT** (dla eksploracji)
   - Podstawowe filtry
   - Confidence > 0.4
   - Zgodność z literaturą ±50%
   - Termodynamika: drift < 2%

---

## 🔧 Komponenty

### 1. MoleculeFilter

**Lokalizacja**: `backend/validation/molecule_filter.py`

**Funkcje**:
- Filtrowanie real molecules vs clusters (używa `filter_real_molecules.py`)
- Walidacja struktury (valence, charge, bond orders)
- Sprawdzanie spójności (formula, atom count)

**Kryteria**:
- Size >= 2
- Bonds >= 0.4 * (size - 1)
- Bond density >= 0.3
- Valence OK
- Charge balance OK

### 2. ThermodynamicValidator

**Lokalizacja**: `backend/sim/core/thermodynamics.py` (już istnieje)

**Użycie**:
- Sprawdzanie logów symulacji
- Weryfikacja energy drift
- Sprawdzanie momentum conservation
- Maxwell-Boltzmann distribution

**Kryteria STRICT**:
- Energy drift < 0.5%
- Momentum drift < 0.1%
- Temperature stability ±5%

### 3. LiteratureValidator

**Lokalizacja**: `backend/validation/literature_validator.py`

**Funkcje**:
- Porównanie z benchmark reactions
- Weryfikacja expected products
- Sprawdzanie yield ranges
- Validation against known reactions

**Kryteria STRICT**:
- ≥50% expected products detected
- Yield within ±20% of literature
- At least 2 benchmark molecules

### 4. MatchConfidenceFilter

**Lokalizacja**: `matcher/confidence.py` (już istnieje)

**Użycie**:
- Filtrowanie dopasowań PubChem
- Sprawdzanie confidence scores
- Walidacja chemical plausibility

**Kryteria STRICT**:
- Confidence > 0.8
- Reliability = HIGH
- Validation status = PASS
- No valence violations

### 5. SimulationQualityFilter

**Lokalizacja**: `backend/validation/simulation_quality.py`

**Funkcje**:
- Sprawdzanie completion rate
- Weryfikacja stability
- Sprawdzanie performance
- Detection of crashes/errors

**Kryteria STRICT**:
- Completion: 100% (all steps)
- No crashes
- Memory stable
- Performance acceptable

### 6. TruthReportGenerator

**Lokalizacja**: `backend/validation/truth_report.py`

**Funkcje**:
- Generowanie raportów walidacji
- Statystyki filtrowania
- Warnings i recommendations
- JSON + Markdown output

---

## 📊 Pipeline Przetwarzania

```
Input: results/phase2b_additional/*/run_X/
│
├─ 1. Load Results
│  ├─ results.json
│  ├─ molecules.json
│  ├─ snapshots/
│  └─ simulation.log
│
├─ 2. Simulation Quality Check
│  ├─ Completion rate
│  ├─ Stability
│  └─ Performance
│
├─ 3. Thermodynamic Validation
│  ├─ Energy conservation
│  ├─ Momentum conservation
│  └─ Temperature stability
│
├─ 4. Molecule Filtering
│  ├─ Real molecules vs clusters
│  ├─ Structure validation
│  └─ Chemical plausibility
│
├─ 5. Literature Validation
│  ├─ Expected products
│  ├─ Yield comparison
│  └─ Benchmark reactions
│
├─ 6. Match Confidence Filter
│  ├─ PubChem matches
│  ├─ Confidence scores
│  └─ Reliability checks
│
└─ 7. Generate Report
   ├─ truth_report.json
   ├─ truth_report.md
   └─ filtered_results.json
```

---

## 🎯 Interfejs API

### TruthFilter Class

```python
from backend.validation.truth_filter import TruthFilter

# Initialize
filter = TruthFilter(
    validation_level="STRICT",  # STRICT, MEDIUM, LENIENT
    output_dir="validation_reports"
)

# Filter single run
result = filter.filter_run(
    run_path="results/phase2b_additional/miller_urey_extended/run_1"
)

# Filter batch
results = filter.filter_batch(
    results_dir="results/phase2b_additional",
    scenario="miller_urey_extended"
)

# Get report
report = filter.generate_report(results)
```

### Output Structure

```json
{
  "run_id": "miller_urey_extended/run_1",
  "validation_level": "STRICT",
  "overall_status": "PASS" | "WARNING" | "FAIL",
  "filters": {
    "simulation_quality": {
      "status": "PASS",
      "completion_rate": 1.0,
      "stability": "stable",
      "warnings": []
    },
    "thermodynamics": {
      "status": "PASS",
      "energy_drift": 0.003,
      "momentum_drift": 0.0001,
      "warnings": []
    },
    "molecules": {
      "status": "PASS",
      "total_molecules": 45,
      "real_molecules": 32,
      "clusters_removed": 13,
      "warnings": []
    },
    "literature": {
      "status": "WARNING",
      "expected_products_detected": 0.4,
      "yield_match_rate": 0.6,
      "warnings": ["Low expected products rate"]
    },
    "match_confidence": {
      "status": "PASS",
      "high_confidence_matches": 15,
      "medium_confidence_matches": 10,
      "low_confidence_matches": 7,
      "warnings": []
    }
  },
  "filtered_results": {
    "molecules": [...],  # Only validated molecules
    "reactions": [...],  # Only validated reactions
    "matches": [...]     # Only high-confidence matches
  }
}
```

---

## 📁 Struktura Plików

```
backend/validation/
├── __init__.py
├── truth_filter.py          # Main TruthFilter class
├── molecule_filter.py        # Molecule filtering
├── literature_validator.py   # Literature validation
├── simulation_quality.py     # Simulation quality checks
└── truth_report.py           # Report generation

scripts/
├── run_truth_filter.py       # CLI script
└── batch_truth_filter.py     # Batch processing

docs/validation/
├── TRUTH_FILTER_GUIDE.md     # User guide
└── TRUTH_FILTER_EXAMPLES.md  # Usage examples
```

---

## 🚀 Plan Implementacji

### Faza 1: Core Infrastructure (Agent-Architect)
- [ ] Projekt interfejsu API
- [ ] Struktura modułów
- [ ] Integracja z istniejącymi komponentami

### Faza 2: Molecule Filtering (Agent-Implementer)
- [ ] MoleculeFilter class
- [ ] Integracja z filter_real_molecules.py
- [ ] Structure validation

### Faza 3: Literature Validation (Agent-Implementer)
- [ ] LiteratureValidator class
- [ ] Integracja z benchmark_reactions.py
- [ ] Yield comparison logic

### Faza 4: Simulation Quality (Agent-Implementer)
- [ ] SimulationQualityFilter class
- [ ] Log parsing
- [ ] Quality metrics

### Faza 5: Report Generation (Agent-Implementer)
- [ ] TruthReportGenerator class
- [ ] JSON + Markdown output
- [ ] Statistics aggregation

### Faza 6: Integration & Testing (Agent-Reviewer)
- [ ] Integration tests
- [ ] Test on Phase 2B results
- [ ] Documentation

---

## 📝 Przykłady Użycia

### Przykład 1: Filter Single Run

```python
from backend.validation.truth_filter import TruthFilter

filter = TruthFilter(validation_level="STRICT")
result = filter.filter_run("results/phase2b_additional/miller_urey_extended/run_1")

if result.overall_status == "PASS":
    print("✅ Run passed validation")
    # Use filtered_results for publication
else:
    print(f"⚠️ Run has warnings: {result.warnings}")
```

### Przykład 2: Filter All Phase 2B Runs

```python
filter = TruthFilter(validation_level="STRICT")
results = filter.filter_batch("results/phase2b_additional")

# Get summary
summary = filter.generate_summary(results)
print(f"Passed: {summary.passed_runs}")
print(f"Warnings: {summary.warning_runs}")
print(f"Failed: {summary.failed_runs}")
```

### Przykład 3: CLI Usage

```bash
# Filter single run
python scripts/run_truth_filter.py \
    --input results/phase2b_additional/miller_urey_extended/run_1 \
    --level STRICT \
    --output validation_reports/run_1

# Filter all runs
python scripts/batch_truth_filter.py \
    --input results/phase2b_additional \
    --level STRICT \
    --output validation_reports/phase2b
```

---

## 🔍 Integracja z Istniejącymi Komponentami

### Używa:
- `scripts/filter_real_molecules.py` → MoleculeFilter
- `backend/sim/core/thermodynamics.py` → ThermodynamicValidator
- `backend/sim/core/benchmark_reactions.py` → LiteratureValidator
- `matcher/confidence.py` → MatchConfidenceFilter

### Rozszerza:
- Dodaje pipeline przetwarzania
- Dodaje poziomy walidacji
- Dodaje raportowanie
- Dodaje batch processing

---

## ✅ Success Criteria

1. **Funkcjonalność**:
   - Wszystkie filtry działają
   - Pipeline przetwarza wyniki
   - Raporty są generowane

2. **Jakość**:
   - Testy przechodzą
   - Dokumentacja kompletna
   - Przykłady działają

3. **Integracja**:
   - Działa z Phase 2B results
   - Kompatybilny z istniejącymi skryptami
   - Nie psuje istniejącego kodu

---

## 📚 Dokumentacja

- **User Guide**: `docs/validation/TRUTH_FILTER_GUIDE.md`
- **API Reference**: Docstrings w kodzie
- **Examples**: `docs/validation/TRUTH_FILTER_EXAMPLES.md`

---

**Status**: Ready for implementation
**Priority**: High (needed for publication)
**Estimated Time**: 2-3 days

