---
date: 2025-11-28
label: guide
---

# Truth-Filter: Przykłady Użycia

## Przykład 1: Filter Single Run (STRICT)

```bash
python scripts/run_truth_filter.py \
    --input results/phase2b_additional/miller_urey_extended/run_1 \
    --level STRICT \
    --output validation_reports/run_1
```

**Output:**
```
TRUTH-FILTER VALIDATION RESULT
======================================================================
Run: run_1
Validation Level: STRICT
Overall Status: PASS

Filter Results:
  ✅ simulation_quality: PASS
  ✅ thermodynamics: PASS
  ✅ molecule_filter: PASS
  ✅ literature_validator: PASS
  ✅ match_confidence: PASS

Summary:
  Total Molecules: 100
  Filtered Molecules: 75
  Retention Rate: 75.0%

Report Files:
  json: validation_reports/run_1/truth_report_run_1.json
  markdown: validation_reports/run_1/truth_report_run_1.md
======================================================================
```

## Przykład 2: Filter Batch (MEDIUM)

```bash
python scripts/batch_truth_filter.py \
    --input results/phase2b_additional \
    --scenario miller_urey_extended \
    --level MEDIUM \
    --output validation_reports/phase2b
```

**Output:**
```
TRUTH-FILTER BATCH VALIDATION SUMMARY
======================================================================
Total Runs: 18
Passed: 15 (83.3%)
Warnings: 3
Failed: 0

Filter Statistics:
  simulation_quality:
    Pass: 18
    Warning: 0
    Fail: 0
  molecule_filter:
    Pass: 15
    Warning: 3
    Fail: 0
  ...

Molecule Statistics:
  Total Original: 1800
  Total Filtered: 1350
  Retention Rate: 75.0%

Summary saved to: validation_reports/phase2b/truth_filter_summary.md
======================================================================
```

## Przykład 3: Python API

```python
from backend.validation import TruthFilter, ValidationLevel

# Initialize
filter = TruthFilter(
    validation_level=ValidationLevel.STRICT,
    output_dir="validation_reports"
)

# Filter single run
result = filter.filter_run(
    "results/phase2b_additional/miller_urey_extended/run_1"
)

# Check status
if result.is_pass():
    print("✅ Run passed validation")
    # Use filtered_results for publication
    validated_molecules = result.filtered_results['molecules']
else:
    print(f"⚠️ Run has warnings: {result.get_warnings()}")
    print(f"❌ Errors: {result.get_errors()}")

# Generate report
output_files = filter.generate_report(result, format="both")
print(f"Reports saved to: {output_files}")

# Filter batch
results = filter.filter_batch(
    "results/phase2b_additional",
    scenario="miller_urey_extended"
)

# Generate summary
summary = filter.generate_summary(results)
print(f"Pass rate: {summary['pass_rate']:.1%}")
```

## Przykład 4: Custom Validation Level

```python
from backend.validation import TruthFilter, ValidationLevel

# LENIENT for exploration
filter = TruthFilter(validation_level=ValidationLevel.LENIENT)

result = filter.filter_run("results/test_run")

# More permissive thresholds
# - Confidence > 0.4
# - Literature match ±50%
# - Energy drift < 2%
```

## Przykład 5: Accessing Filter Results

```python
result = filter.filter_run("results/phase2b_additional/.../run_1")

# Access individual filter results
sim_quality = result.filters['simulation_quality']
print(f"Completion rate: {sim_quality.details['completion_rate']:.1%}")

mol_filter = result.filters['molecule_filter']
print(f"Retention rate: {mol_filter.details['retention_rate']:.1%}")

lit_validator = result.filters['literature_validator']
print(f"Detection rate: {lit_validator.details['detection_rate']:.1%}")
```

## Przykład 6: Batch Processing with Summary

```python
from backend.validation import TruthFilter, ValidationLevel

filter = TruthFilter(validation_level=ValidationLevel.MEDIUM)

# Filter all scenarios
scenarios = ['miller_urey_extended', 'hydrothermal_extended', 'formamide_extended']
all_results = []

for scenario in scenarios:
    results = filter.filter_batch(
        "results/phase2b_additional",
        scenario=scenario
    )
    all_results.extend(results)

# Generate overall summary
summary = filter.generate_summary(all_results)
print(f"Overall pass rate: {summary['pass_rate']:.1%}")
print(f"Total runs: {summary['total_runs']}")
```

---

**Status**: ✅ Ready to use

