---
date: 2025-11-28
label: guide
---

# Truth-Filter: User Guide

## 🎯 Wprowadzenie

Truth-Filter to system walidacji wyników symulacji przed publikacją. Zapewnia, że tylko wiarygodne, naukowo poprawne wyniki trafiają do manuskryptu.

## 📋 Funkcje

Truth-Filter wykonuje następujące walidacje:

1. **Simulation Quality** - Sprawdza completion rate, stability, performance
2. **Thermodynamics** - Weryfikuje energy/momentum drift
3. **Molecule Filtering** - Filtruje real molecules vs clusters
4. **Literature Validation** - Porównuje z benchmark reactions
5. **Match Confidence** - Sprawdza wiarygodność dopasowań PubChem

## 🚀 Quick Start

### Filter Single Run

```bash
python scripts/run_truth_filter.py \
    --input results/phase2b_additional/miller_urey_extended/run_1 \
    --level STRICT \
    --output validation_reports/run_1
```

### Filter Batch

```bash
python scripts/batch_truth_filter.py \
    --input results/phase2b_additional \
    --scenario miller_urey_extended \
    --level STRICT \
    --output validation_reports/phase2b
```

## 📊 Poziomy Walidacji

### STRICT (dla publikacji)
- Wszystkie filtry muszą przejść
- Confidence > 0.8 dla dopasowań
- Zgodność z literaturą ±20%
- Termodynamika: drift < 0.5%

### MEDIUM (dla analizy)
- Większość filtrów musi przejść
- Confidence > 0.6
- Zgodność z literaturą ±30%
- Termodynamika: drift < 1%

### LENIENT (dla eksploracji)
- Podstawowe filtry
- Confidence > 0.4
- Zgodność z literaturą ±50%
- Termodynamika: drift < 2%

## 🔧 Użycie Programatyczne

### Python API

```python
from backend.validation import TruthFilter, ValidationLevel

# Initialize
filter = TruthFilter(
    validation_level=ValidationLevel.STRICT,
    output_dir="validation_reports"
)

# Filter single run
result = filter.filter_run("results/phase2b_additional/miller_urey_extended/run_1")

# Check status
if result.is_pass():
    print("✅ Run passed validation")
    # Use filtered_results for publication
else:
    print(f"⚠️ Run has warnings: {result.get_warnings()}")

# Generate report
output_files = filter.generate_report(result, format="both")

# Filter batch
results = filter.filter_batch(
    "results/phase2b_additional",
    scenario="miller_urey_extended"
)

# Generate summary
summary = filter.generate_summary(results)
```

## 📁 Struktura Wyników

### TruthResult

```python
{
    'run_id': 'run_1',
    'run_path': 'results/.../run_1',
    'validation_level': 'STRICT',
    'overall_status': 'PASS',
    'filters': {
        'simulation_quality': {...},
        'thermodynamics': {...},
        'molecule_filter': {...},
        'literature_validator': {...},
        'match_confidence': {...}
    },
    'filtered_results': {
        'molecules': [...],  # Only validated molecules
        'reactions': [...],
        'matches': [...]
    },
    'summary': {
        'total_molecules': 100,
        'filtered_molecules': 75,
        'retention_rate': 0.75
    }
}
```

## 📄 Raporty

Truth-Filter generuje raporty w dwóch formatach:

### JSON Report (`truth_report_*.json`)
- Strukturalne dane
- Łatwe do przetwarzania programatycznego
- Zawiera wszystkie szczegóły walidacji

### Markdown Report (`truth_report_*.md`)
- Human-readable format
- Łatwe do przeglądania
- Zawiera statusy, warnings, errors

## ⚠️ Troubleshooting

### Problem: "Benchmark database not found"
**Rozwiązanie**: Upewnij się, że plik `data/benchmark_reactions.json` istnieje.

### Problem: "No molecules found"
**Rozwiązanie**: Sprawdź czy `molecules.json` lub `results.json` istnieje w run directory.

### Problem: "Match confidence evaluator not available"
**Rozwiązanie**: To nie jest błąd - match confidence filtering jest opcjonalne. Walidacja będzie kontynuowana bez tego filtra.

## 🔍 Przykłady

Zobacz `docs/validation/TRUTH_FILTER_EXAMPLES.md` dla więcej przykładów użycia.

---

**Status**: ✅ Ready to use
**Version**: 1.0.0

