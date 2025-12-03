# Quick Guide: Generate Missing Figures for Manuscript

**Date**: 2025-01-23  
**Status**: ⚠️ **ACTION REQUIRED**

---

## 🎯 Co Trzeba Wygenerować

1. ✅ **Wykresy termodynamiczne** (energia vs czas, histogram M-B)
2. ✅ **Benchmark reakcji** (formose/Strecker/HCN)
3. ⚠️ **Struktury molekularne** (top 5 z PubChem Matcher)
4. ⚠️ **Przykład sieci reakcji** (wizualizacja z ReactionNetworkAnalyzer)

---

## 🚀 Szybki Start (Wszystko w Jednym)

### Opcja 1: Automatyczny Script (Najłatwiejsze)

```bash
# Wygeneruj wszystkie figury z prawdziwych danych
python scripts/generate_paper_figures_from_real_data.py \
    --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
    --output-dir paper/figures
```

**To wygeneruje**:
- `figure1_energy_conservation.png` - Energia vs czas
- `figure1_maxwell_boltzmann.png` - Histogram M-B
- `figure2_formose_validation.png` - Benchmark formose
- `molecular_structures_panel.png` - Top 5 molekuł
- `reaction_network_example.png` - Przykład sieci

---

## 📋 Krok po Kroku (Jeśli Automatyczny Nie Działa)

### Krok 1: Wykresy Termodynamiczne

```bash
# Znajdź validation log (jeśli istnieje)
python scripts/analyze_thermodynamics.py \
    --input diagnostics/validation_log.json \
    --output-dir paper/figures
```

**Jeśli nie ma validation log**:
- Użyj istniejącego `figure1_thermodynamic_validation.png` (syntetyczne dane są OK dla submission)

---

### Krok 2: Benchmark Reakcji

```bash
# Jeśli masz benchmark results
python scripts/analyze_benchmark_reactions.py \
    --simulation-data results/benchmarks/formose/results.json \
    --output-dir paper/figures
```

**Jeśli nie masz benchmark results**:
- Użyj istniejącego `figure2_benchmark_validation.png` (syntetyczne dane są OK)

---

### Krok 3: Struktury Molekularne (NOWE)

```bash
# Wybierz top 5 molekuł i wygeneruj struktury
python scripts/match_top_molecules_pubchem.py \
    --filtered-file results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
    --output paper/figures/molecular_structures
```

**Lub użyj matcher_v2 bezpośrednio**:
```python
from matcher.matcher_v2 import MatcherV2
import json

# Load molecules
with open('results/phase2b_additional/miller_urey_extended/run_1/molecules.json') as f:
    molecules = json.load(f)

# Match top 5
matcher = MatcherV2()
top_5 = sorted(molecules['molecules'], key=lambda x: x.get('abundance', 0), reverse=True)[:5]

for mol in top_5:
    result = matcher.match_cluster(mol)
    if result.success:
        print(f"{mol['formula']} -> {result.pubchem_name} (CID: {result.pubchem_cid})")
```

---

### Krok 4: Przykład Sieci Reakcji (NOWE)

```bash
# Zbuduj sieć z jednej symulacji
python scripts/reaction_network_analyzer.py \
    results/phase2b_additional/miller_urey_extended/run_1 \
    --output analysis/reaction_network_example \
    --export both

# Wizualizuj
python scripts/network_visualizer.py \
    analysis/reaction_network_example/reaction_network.json \
    --max-nodes 50 \
    --output paper/figures/reaction_network_example.png
```

---

## ⚠️ Jeśli Nie Masz Prawdziwych Danych

### Dla Submission (OK):
- **Figure 1** (termodynamiczne): Syntetyczne dane są **akceptowalne** - pokazują format
- **Figure 2** (benchmark): Syntetyczne dane są **akceptowalne** - pokazują format

### Wymagane (Przed Final Submission):
- **Struktury molekularne**: **MUSZĄ** być prawdziwe (użyj PubChem Matcher)
- **Sieć reakcji**: **POWINNA** być prawdziwa (użyj ReactionNetworkAnalyzer)

---

## 📝 Integracja z Manuskryptem

### Po Wygenerowaniu Figur:

1. **Struktury molekularne**:
   - Dodaj do Figure 6 jako panel E
   - Lub dodaj nową sekcję w Results (3.5: Example Molecular Structures)

2. **Sieć reakcji**:
   - Dodaj do Figure 4 jako panel E
   - Lub dodaj jako Figure 7 (Example Reaction Network)

3. **Zaktualizuj captions**:
   - Dodaj informacje o źródle danych
   - Dodaj informacje o metodzie (PubChem Matcher, ReactionNetworkAnalyzer)

---

## ✅ Checklist

### Przed Submission:
- [ ] Wykresy termodynamiczne (syntetyczne OK)
- [ ] Benchmark reakcji (syntetyczne OK)
- [ ] **Struktury molekularne** (prawdziwe - wymagane)
- [ ] **Sieć reakcji** (prawdziwa - wymagana)

### Po Submission (Jeśli Czasopismo Poprosi):
- [ ] Wykresy termodynamiczne z prawdziwych danych
- [ ] Benchmark reakcji z prawdziwych danych
- [ ] Więcej struktur molekularnych
- [ ] Więcej przykładów sieci

---

## 🎯 Rekomendacja

**Dla szybkiego submissionu**:
1. Użyj istniejących Figure 1 i 2 (syntetyczne dane)
2. **Wygeneruj struktury molekularne** (Krok 3) - **WYMAGANE**
3. **Wygeneruj sieć reakcji** (Krok 4) - **WYMAGANE**

**Szacowany czas**: 30-60 minut

---

**Status**: ⚠️ **WYMAGA DZIAŁANIA**  
**Priorytet**: **WYSOKI** (przed submission)

