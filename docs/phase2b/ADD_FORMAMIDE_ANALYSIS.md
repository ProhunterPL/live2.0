---
date: 2025-11-27
label: guide
---

# Dodawanie Formamide Extended do analizy Phase 2B

## Szybki start

Po zakończeniu testów formamide_extended, uruchom:

### PowerShell (Windows):
```powershell
.\scripts\add_formamide_analysis.ps1
```

### Python (cross-platform):
```bash
python scripts/add_formamide_to_analysis.py
```

### Pełna analiza (ręcznie):
```bash
python scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data
```

## Co się stanie?

1. **Automatyczne wykrywanie**: Skrypt automatycznie wykryje wszystkie dostępne runy w `results/phase2b_additional/formamide_extended/`

2. **Pełna analiza**: Uruchomi kompletną analizę Phase 2B, która obejmuje:
   - Miller Urey Extended (już przeanalizowane)
   - Hydrothermal Extended (już przeanalizowane)
   - **Formamide Extended** (nowe)

3. **Aktualizacja wyników**: Zaktualizuje wszystkie pliki w `paper/results_data/`:
   - `summary_table.csv` / `summary_table.tex`
   - `formamide_extended_analysis.json`
   - `scenario_comparison.json`
   - `figure_data.json`
   - `latex_snippets.txt`

## Wymagania

- Katalog `results/phase2b_additional/formamide_extended/` musi istnieć
- Każdy run musi mieć plik `results.json` (symulacja zakończona)
- Opcjonalnie: `reaction_network.json` dla każdego runa (lepsze wyniki)

## Sprawdzenie gotowości

Przed uruchomieniem możesz sprawdzić ile runów jest gotowych:

```powershell
# PowerShell
$runs = Get-ChildItem "results/phase2b_additional/formamide_extended" -Directory -Filter "run_*"
$completed = ($runs | Where-Object { Test-Path (Join-Path $_.FullName "results.json") }).Count
Write-Host "Completed runs: $completed"
```

```bash
# Bash
find results/phase2b_additional/formamide_extended -name "results.json" | wc -l
```

## Wyniki

Po zakończeniu analizy, formamide_extended będzie widoczny w:

- **Summary table**: `paper/results_data/summary_table.csv`
- **Scenario comparison**: `paper/results_data/scenario_comparison.json`
- **Individual analysis**: `paper/results_data/formamide_extended_analysis.json`

## Uwagi

- Analiza automatycznie wykrywa wszystkie dostępne runy (nie trzeba podawać liczby)
- Jeśli niektóre runy nie są jeszcze gotowe, zostaną pominięte (można uruchomić ponownie później)
- Analiza może zająć kilka godzin dla dużej liczby runów (zależy od liczby cykli w sieciach)

## Troubleshooting

### "Scenario directory not found"
- Sprawdź czy katalog `results/phase2b_additional/formamide_extended/` istnieje
- Sprawdź czy ścieżka jest poprawna

### "No completed runs found"
- Sprawdź czy runy mają plik `results.json`
- Sprawdź czy symulacje rzeczywiście się zakończyły

### Timeout podczas wykrywania cykli
- To normalne dla gęstych sieci
- Analiza kontynuuje z częściowymi wynikami
- Można zwiększyć timeout w `autocatalysis_detector.py` (domyślnie 600s)

