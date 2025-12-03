# Plan Dodania Brakujących Elementów do Artykułu

**Date**: 2025-01-23  
**Status**: ⚠️ **ACTION REQUIRED**

---

## 📋 Wymagane Elementy

### 1. ✅ Wykresy Termodynamiczne (1-2)
- **Status**: Figure 1 już istnieje, ale używa syntetycznych danych
- **Potrzebne**: Prawdziwe dane z symulacji
- **Elementy**:
  - Energia całkowita vs czas
  - Histogram Maxwell-Boltzmann

### 2. ✅ Benchmark Reakcji (1)
- **Status**: Figure 2 już istnieje, ale używa syntetycznych danych
- **Potrzebne**: Prawdziwe dane z benchmark symulacji
- **Opcje**: Formose / Strecker / HCN polymerization

### 3. ⚠️ Wykryte Struktury Molekularne (kilka)
- **Status**: Brak w manuskrypcie
- **Narzędzie**: PubChem Matcher (matcher_v2)
- **Potrzebne**: Wybrać top 3-5 molekuł i pokazać struktury

### 4. ⚠️ Przykład Sieci Reakcji
- **Status**: Figure 4 istnieje, ale może potrzebować konkretnego przykładu
- **Narzędzie**: ReactionNetworkAnalyzer
- **Potrzebne**: Wizualizacja konkretnej sieci z jednej symulacji

---

## 🎯 Plan Działania

### Krok 1: Wykresy Termodynamiczne z Prawdziwych Danych

**Plik**: `scripts/analyze_thermodynamics.py` (już istnieje)

**Akcja**:
```bash
# Znajdź wyniki symulacji z validation log
# Uruchom analizę termodynamiczną
python scripts/analyze_thermodynamics.py \
    --input diagnostics/validation_log.json \
    --output-dir paper/figures
```

**Output**:
- `fig1_energy_conservation.png` - Energia vs czas
- `fig2_maxwell_boltzmann.png` - Histogram M-B

**Integracja z manuskryptem**:
- Zastąp `figure1_thermodynamic_validation.png` prawdziwymi wykresami
- Upewnij się, że dane są z rzeczywistej symulacji (nie syntetyczne)

---

### Krok 2: Benchmark Reakcji z Prawdziwych Danych

**Plik**: `scripts/analyze_benchmark_reactions.py` (już istnieje)

**Akcja**:
```bash
# Uruchom benchmark symulacje (jeśli nie były uruchomione)
python scripts/run_benchmark_simulations.py \
    --scenario formose \
    --output results/benchmarks/formose

# Analizuj wyniki
python scripts/analyze_benchmark_reactions.py \
    --simulation-data results/benchmarks/formose/results.json \
    --output-dir paper/figures
```

**Output**:
- `figure3_formose_validation.png` - Formose reaction validation

**Integracja z manuskryptem**:
- Zastąp `figure2_benchmark_validation.png` prawdziwymi danymi
- Lub dodaj jako dodatkowy panel do Figure 2

---

### Krok 3: Wykryte Struktury Molekularne

**Narzędzie**: PubChem Matcher (matcher_v2)

**Akcja**:
```bash
# 1. Wybierz top molekuły z wyników Phase 2B
python scripts/match_top_molecules_pubchem.py \
    --filtered-file results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
    --output paper/figures/molecular_structures

# 2. Wygeneruj struktury dla top 5 molekuł
python matcher/matcher_v2.py \
    --input results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
    --top-n 5 \
    --output paper/figures/molecular_structures
```

**Output**:
- Struktury molekularne (PNG) dla top 5 molekuł
- Porównanie z PubChem matches

**Integracja z manuskryptem**:
- Dodać nową sekcję w Results (3.5: Example Molecular Structures)
- Lub dodać do Figure 6 (Novel Molecules) jako panel E

---

### Krok 4: Przykład Sieci Reakcji

**Narzędzie**: ReactionNetworkAnalyzer

**Akcja**:
```bash
# 1. Zbuduj sieć reakcji z jednej symulacji
python scripts/reaction_network_analyzer.py \
    results/phase2b_additional/miller_urey_extended/run_1 \
    --output analysis/reaction_network_example \
    --export both

# 2. Wizualizuj sieć
python scripts/network_visualizer.py \
    analysis/reaction_network_example/reaction_network.json \
    --max-nodes 50 \
    --output paper/figures/reaction_network_example.png
```

**Output**:
- `reaction_network_example.png` - Wizualizacja sieci
- `reaction_network.json` - Dane sieci

**Integracja z manuskryptem**:
- Dodać do Figure 4 jako panel E (Example Network)
- Lub dodać jako nową figurę (Figure 7)

---

## 📝 Szczegółowy Plan Implementacji

### Task 1: Prawdziwe Wykresy Termodynamiczne

**Status**: ⚠️ Wymaga prawdziwych danych z symulacji

**Kroki**:
1. Znajdź wyniki symulacji z validation log
2. Uruchom `scripts/analyze_thermodynamics.py`
3. Zastąp syntetyczne wykresy prawdziwymi
4. Zaktualizuj caption w manuskrypcie

**Pliki do modyfikacji**:
- `paper/manuscript_draft.tex` (Figure 1 caption)
- `paper/figures/figure1_thermodynamic_validation.png` (zastąp)

---

### Task 2: Prawdziwe Benchmark Reakcji

**Status**: ⚠️ Wymaga uruchomienia benchmark symulacji

**Kroki**:
1. Uruchom benchmark symulacje (formose/Strecker/HCN)
2. Analizuj wyniki
3. Wygeneruj wykresy
4. Zastąp syntetyczne dane prawdziwymi

**Pliki do modyfikacji**:
- `paper/manuscript_draft.tex` (Figure 2 caption)
- `paper/figures/figure2_benchmark_validation.png` (zastąp)

---

### Task 3: Struktury Molekularne

**Status**: ⚠️ Nowy element - wymaga implementacji

**Kroki**:
1. Wybierz top 5 molekuł z Phase 2B results
2. Uruchom PubChem Matcher dla każdej
3. Wygeneruj struktury (PNG)
4. Utwórz panel z strukturami
5. Dodaj do manuskryptu

**Nowe pliki**:
- `paper/figures/molecular_structures_panel.png` (nowy)
- `paper/manuscript_draft.tex` (nowa sekcja lub panel)

**Script do utworzenia**:
- `scripts/generate_molecular_structures_panel.py` (nowy)

---

### Task 4: Przykład Sieci Reakcji

**Status**: ⚠️ Wymaga wygenerowania z prawdziwych danych

**Kroki**:
1. Wybierz jedną symulację (np. miller_urey_extended/run_1)
2. Zbuduj sieć reakcji
3. Wizualizuj (top 50 molekuł)
4. Dodaj do manuskryptu

**Pliki do modyfikacji**:
- `paper/figures/figure4_reaction_networks.png` (zastąp lub dodaj panel)
- `paper/manuscript_draft.tex` (Figure 4 caption)

---

## 🔧 Narzędzia Gotowe

### ✅ Dostępne Skrypty

1. **Thermodynamic Analysis**:
   - `scripts/analyze_thermodynamics.py` ✅
   - Generuje wykresy energii i M-B

2. **Benchmark Reactions**:
   - `scripts/analyze_benchmark_reactions.py` ✅
   - `scripts/run_benchmark_simulations.py` ✅

3. **PubChem Matcher**:
   - `matcher/matcher_v2.py` ✅
   - `scripts/match_top_molecules_pubchem.py` ✅

4. **Reaction Network**:
   - `scripts/reaction_network_analyzer.py` ✅
   - `scripts/network_visualizer.py` ✅

---

## 📊 Priorytety

### Wysoki Priorytet (Przed Submission):
1. ✅ **Wykresy termodynamiczne** - Użyj prawdziwych danych
2. ✅ **Benchmark reakcji** - Użyj prawdziwych danych (jeśli dostępne)

### Średni Priorytet (Można dodać później):
3. ⚠️ **Struktury molekularne** - Nowy element, wymaga implementacji
4. ⚠️ **Przykład sieci** - Można użyć istniejącego Figure 4

---

## 🎯 Szybki Start

### Najszybsza Opcja (Użyj Istniejących):

1. **Thermodynamic**: Użyj `scripts/analyze_thermodynamics.py` z prawdziwymi danymi
2. **Benchmark**: Użyj `scripts/analyze_benchmark_reactions.py` z prawdziwymi danymi
3. **Structures**: Użyj `matcher/matcher_v2.py` dla top 5 molekuł
4. **Network**: Użyj `scripts/reaction_network_analyzer.py` + `network_visualizer.py`

---

## 📝 Checklist

### Przed Submission:
- [ ] Wykresy termodynamiczne z prawdziwych danych
- [ ] Benchmark reakcji z prawdziwych danych (jeśli dostępne)
- [ ] Struktury molekularne (top 3-5)
- [ ] Przykład sieci reakcji

### Po Submission (Jeśli Czasopismo Poprosi):
- [ ] Dodatkowe struktury molekularne
- [ ] Więcej przykładów sieci
- [ ] Extended benchmark validation

---

**Status**: ⚠️ **WYMAGA DZIAŁANIA**  
**Priorytet**: **WYSOKI** (przed submission)  
**Szacowany Czas**: 2-4 godziny

