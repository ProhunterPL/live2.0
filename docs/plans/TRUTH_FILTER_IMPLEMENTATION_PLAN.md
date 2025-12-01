---
date: 2025-11-28
label: plan
---

# Truth-Filter: Plan Implementacji

## 📋 Podsumowanie

Truth-Filter to system walidacji wyników przed publikacją. Integruje istniejące komponenty (termodynamika, chemia, dopasowania) w jeden pipeline z trzema poziomami walidacji (STRICT/MEDIUM/LENIENT).

**Cel**: Zapewnić, że tylko wiarygodne wyniki trafiają do publikacji.

**Czas**: 2-3 dni implementacji

---

## 🎯 Plan Działania

### Krok 1: Architektura i Interfejs (Agent-Architect)
**Czas**: 2-3h

**Zadania**:
1. Zaprojektować strukturę modułów `backend/validation/`
2. Zdefiniować interfejs API `TruthFilter` class
3. Zaprojektować strukturę danych (TruthResult, FilterResult)
4. Określić integrację z istniejącymi komponentami

**Output**:
- `backend/validation/__init__.py`
- `backend/validation/truth_filter.py` (szkielet klasy)
- `backend/validation/types.py` (dataclasses)

**Komendy dla Agent-Architect**:
```
Zaprojektuj moduł backend/validation/ z klasą TruthFilter:
- Interfejs: filter_run(run_path) -> TruthResult
- Interfejs: filter_batch(results_dir) -> List[TruthResult]
- Poziomy walidacji: STRICT, MEDIUM, LENIENT
- Integracja z: filter_real_molecules, thermodynamics, benchmark_reactions, matcher/confidence
- Output: JSON + Markdown reports
```

---

### Krok 2: Molecule Filter (Agent-Implementer)
**Czas**: 3-4h

**Zadania**:
1. Stworzyć `MoleculeFilter` class
2. Zintegrować z `scripts/filter_real_molecules.py`
3. Dodać structure validation (valence, charge, bonds)
4. Dodać chemical plausibility checks

**Output**:
- `backend/validation/molecule_filter.py`

**Komendy dla Agent-Implementer**:
```
Zaimplementuj backend/validation/molecule_filter.py:
- Klasa MoleculeFilter z metodą filter(molecules) -> filtered_molecules
- Użyj logiki z scripts/filter_real_molecules.py (is_real_molecule)
- Dodaj walidację: valence_check, charge_balance, bond_orders
- Zwróć: {filtered_molecules, clusters_removed, warnings}
- Testy: test_molecule_filter.py
```

---

### Krok 3: Literature Validator (Agent-Implementer)
**Czas**: 3-4h

**Zadania**:
1. Stworzyć `LiteratureValidator` class
2. Zintegrować z `backend/sim/core/benchmark_reactions.py`
3. Dodać yield comparison logic
4. Dodać expected products detection

**Output**:
- `backend/validation/literature_validator.py`

**Komendy dla Agent-Implementer**:
```
Zaimplementuj backend/validation/literature_validator.py:
- Klasa LiteratureValidator z metodą validate(molecules, scenario) -> validation_result
- Użyj BenchmarkReactionDatabase z benchmark_reactions.py
- Sprawdź expected products dla scenariusza (miller_urey, hydrothermal, formamide)
- Porównaj yield z literaturą (tolerance zależny od poziomu walidacji)
- Zwróć: {expected_products_detected, yield_match_rate, warnings}
- Testy: test_literature_validator.py
```

---

### Krok 4: Simulation Quality Filter (Agent-Implementer)
**Czas**: 2-3h

**Zadania**:
1. Stworzyć `SimulationQualityFilter` class
2. Dodać log parsing (completion, stability, performance)
3. Dodać quality metrics calculation
4. Dodać error detection

**Output**:
- `backend/validation/simulation_quality.py`

**Komendy dla Agent-Implementer**:
```
Zaimplementuj backend/validation/simulation_quality.py:
- Klasa SimulationQualityFilter z metodą check(run_path) -> quality_result
- Parsuj simulation.log: completion rate, crashes, errors
- Sprawdź results.json: final_step, max_steps, stability metrics
- Zwróć: {completion_rate, stability, performance, warnings}
- Testy: test_simulation_quality.py
```

---

### Krok 5: Thermodynamic Validator Integration (Agent-Implementer)
**Czas**: 2-3h

**Zadania**:
1. Zintegrować `ThermodynamicValidator` z Truth-Filter
2. Dodać log parsing dla thermodynamic metrics
3. Dodać threshold checking (zależny od poziomu)
4. Dodać warnings generation

**Output**:
- Rozszerzenie `backend/validation/truth_filter.py`

**Komendy dla Agent-Implementer**:
```
Zintegruj ThermodynamicValidator z Truth-Filter:
- W backend/validation/truth_filter.py dodaj metodę _validate_thermodynamics()
- Parsuj simulation.log dla energy/momentum drift
- Użyj thresholdów zależnych od poziomu walidacji:
  - STRICT: energy < 0.5%, momentum < 0.1%
  - MEDIUM: energy < 1%, momentum < 0.5%
  - LENIENT: energy < 2%, momentum < 1%
- Zwróć: {energy_drift, momentum_drift, temperature_stability, warnings}
```

---

### Krok 6: Match Confidence Filter Integration (Agent-Implementer)
**Czas**: 2-3h

**Zadania**:
1. Zintegrować `MatchConfidenceEvaluator` z Truth-Filter
2. Dodać filtering based on confidence scores
3. Dodać reliability checks
4. Dodać warnings for low-confidence matches

**Output**:
- Rozszerzenie `backend/validation/truth_filter.py`

**Komendy dla Agent-Implementer**:
```
Zintegruj MatchConfidenceEvaluator z Truth-Filter:
- W backend/validation/truth_filter.py dodaj metodę _filter_matches()
- Użyj MatchConfidenceEvaluator z matcher/confidence.py
- Filtruj dopasowania według poziomu:
  - STRICT: confidence > 0.8, reliability = HIGH
  - MEDIUM: confidence > 0.6, reliability >= MEDIUM
  - LENIENT: confidence > 0.4
- Zwróć: {high_confidence, medium_confidence, low_confidence, warnings}
```

---

### Krok 7: Report Generator (Agent-Implementer)
**Czas**: 3-4h

**Zadania**:
1. Stworzyć `TruthReportGenerator` class
2. Dodać JSON output
3. Dodać Markdown output
4. Dodać statistics aggregation
5. Dodać summary generation

**Output**:
- `backend/validation/truth_report.py`

**Komendy dla Agent-Implementer**:
```
Zaimplementuj backend/validation/truth_report.py:
- Klasa TruthReportGenerator z metodami:
  - generate_report(result) -> JSON + Markdown
  - generate_summary(results) -> Summary statistics
- Format JSON: truth_report.json (struktura z TRUTH_FILTER_DESIGN.md)
- Format Markdown: truth_report.md (human-readable)
- Statystyki: passed/warning/failed runs, retention rates, etc.
- Testy: test_truth_report.py
```

---

### Krok 8: Main TruthFilter Integration (Agent-Implementer)
**Czas**: 4-5h

**Zadania**:
1. Zaimplementować główną klasę `TruthFilter`
2. Zintegrować wszystkie filtry w pipeline
3. Dodać batch processing
4. Dodać error handling
5. Dodać logging

**Output**:
- `backend/validation/truth_filter.py` (pełna implementacja)

**Komendy dla Agent-Implementer**:
```
Zaimplementuj backend/validation/truth_filter.py:
- Klasa TruthFilter z metodami:
  - filter_run(run_path) -> TruthResult
  - filter_batch(results_dir, scenario) -> List[TruthResult]
  - generate_report(results) -> Report
- Pipeline: simulation_quality → thermodynamics → molecules → literature → matches
- Poziomy walidacji: STRICT, MEDIUM, LENIENT (configurable)
- Error handling: graceful failures, warnings
- Logging: użyj logging.getLogger(__name__)
- Testy: test_truth_filter.py
```

---

### Krok 9: CLI Scripts (Agent-Implementer)
**Czas**: 2-3h

**Zadania**:
1. Stworzyć `scripts/run_truth_filter.py` (single run)
2. Stworzyć `scripts/batch_truth_filter.py` (batch processing)
3. Dodać argument parsing
4. Dodać progress reporting

**Output**:
- `scripts/run_truth_filter.py`
- `scripts/batch_truth_filter.py`

**Komendy dla Agent-Implementer**:
```
Zaimplementuj CLI scripts:
- scripts/run_truth_filter.py:
  - --input: path to run directory
  - --level: STRICT|MEDIUM|LENIENT
  - --output: output directory
  - --format: json|markdown|both
- scripts/batch_truth_filter.py:
  - --input: results directory
  - --scenario: optional scenario filter
  - --level: STRICT|MEDIUM|LENIENT
  - --output: output directory
  - --parallel: optional parallel processing
```

---

### Krok 10: Testing & Documentation (Agent-Reviewer)
**Czas**: 4-5h

**Zadania**:
1. Napisać testy jednostkowe
2. Napisać testy integracyjne
3. Przetestować na Phase 2B results
4. Napisać dokumentację użytkownika
5. Dodać przykłady użycia

**Output**:
- `tests/test_truth_filter.py`
- `docs/validation/TRUTH_FILTER_GUIDE.md`
- `docs/validation/TRUTH_FILTER_EXAMPLES.md`

**Komendy dla Agent-Reviewer**:
```
Przetestuj i udokumentuj Truth-Filter:
- Testy jednostkowe: wszystkie moduły
- Testy integracyjne: pełny pipeline na przykładowych danych
- Test na Phase 2B: przetestuj na 2-3 runach z results/phase2b_additional/
- Dokumentacja:
  - docs/validation/TRUTH_FILTER_GUIDE.md (user guide)
  - docs/validation/TRUTH_FILTER_EXAMPLES.md (usage examples)
- Sprawdź: czy wszystkie komponenty działają razem
```

---

## 📊 Harmonogram

| Krok | Agent | Czas | Priorytet |
|------|-------|------|-----------|
| 1. Architektura | Architect | 2-3h | Wysoki |
| 2. Molecule Filter | Implementer | 3-4h | Wysoki |
| 3. Literature Validator | Implementer | 3-4h | Wysoki |
| 4. Simulation Quality | Implementer | 2-3h | Średni |
| 5. Thermodynamics | Implementer | 2-3h | Średni |
| 6. Match Confidence | Implementer | 2-3h | Średni |
| 7. Report Generator | Implementer | 3-4h | Wysoki |
| 8. Main Integration | Implementer | 4-5h | Wysoki |
| 9. CLI Scripts | Implementer | 2-3h | Średni |
| 10. Testing & Docs | Reviewer | 4-5h | Wysoki |

**Total**: ~25-35 godzin (3-4 dni robocze)

---

## 🎯 Success Criteria

1. ✅ Wszystkie filtry działają
2. ✅ Pipeline przetwarza wyniki
3. ✅ Raporty są generowane (JSON + Markdown)
4. ✅ Testy przechodzą
5. ✅ Dokumentacja kompletna
6. ✅ Działa na Phase 2B results

---

## 🚀 Quick Start (Po Implementacji)

```bash
# Filter single run
python scripts/run_truth_filter.py \
    --input results/phase2b_additional/miller_urey_extended/run_1 \
    --level STRICT \
    --output validation_reports/run_1

# Filter all Phase 2B runs
python scripts/batch_truth_filter.py \
    --input results/phase2b_additional \
    --level STRICT \
    --output validation_reports/phase2b
```

---

**Status**: Ready for implementation
**Next Step**: Agent-Architect - zaprojektuj strukturę modułów

