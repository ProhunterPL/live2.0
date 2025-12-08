---
date: 2025-12-04
label: plan
---

# Paper 2: Plan Szczegółowej Analizy Autocatalytic Cycles

**Status**: 📋 Planning  
**Timeline**: Styczeń-Luty 2026  
**Data Source**: Phase 2B (43 runs, 769,315 cycles)

---

## 📊 Dostępne Dane

### Phase 2B Results
- ✅ **43 runs completed**: 18 Miller-Urey, 17 Hydrothermal, 8 Formamide
- ✅ **769,315 autocatalytic cycles** detected across all scenarios
- ✅ **Reaction networks** analyzed (molecules as nodes, reactions as edges)
- ✅ **All snapshots** saved (every 50K steps)
- ✅ **All checkpoints** saved (every 100K steps)

### Lokalizacja Danych
- `results/phase2b_additional/miller_urey_extended/run_1/` through `run_18/`
- `results/phase2b_additional/hydrothermal_extended/run_1/` through `run_17/`
- `results/phase2b_additional/formamide_extended/run_1/` through `run_8/`

### Pliki z Analizą
- `results/phase2b_additional/hydrothermal_extended/run_X/autocatalytic_cycles.json` - Wykryte cykle
- `results/phase2b_additional/hydrothermal_extended/run_X/complexity_metrics.json` - Metryki złożoności
- `results/phase2b_additional/phase2b_analysis_results.json` - Agregowane wyniki

---

## 🎯 Cele Analizy

### 1. Klasyfikacja Cykli (Direct, Indirect, Hypercycles)
- **Cel**: Podzielić wszystkie 769,315 cykli na kategorie
- **Metoda**: Analiza struktury cykli (reakcje, uczestnicy)
- **Output**: `analysis/phase2b/autocatalytic_cycles_classified.json`

### 2. Amplification Factors
- **Cel**: Obliczyć amplification factors dla każdego cyklu
- **Metoda**: Analiza temporalna (d[A]/dt przed vs po aktywacji cyklu)
- **Output**: `analysis/phase2b/amplification_factors.json`

### 3. Network Topology Analysis
- **Cel**: Analiza topologii sieci reakcji
- **Metodyki**:
  - Degree distribution (in-degree, out-degree, total)
  - Hub identification (top 10% by degree centrality)
  - Centrality metrics (betweenness, closeness, eigenvector)
  - Network motifs (3-node, 4-node)
- **Output**: `analysis/phase2b/network_topology_metrics.json`

### 4. Scenario Comparison
- **Cel**: Statystyczne porównanie między scenariuszami
- **Metodyki**:
  - ANOVA dla amplification factors
  - Chi-square dla typów cykli
  - T-tests dla topology metrics
- **Output**: `analysis/phase2b/scenario_comparison.json`

---

## 📝 Plan Implementacji

### Krok 1: Przygotowanie Danych (Tydzień 1)

**Zadania**:
- [ ] Zweryfikować dostępność wszystkich plików `autocatalytic_cycles.json`
- [ ] Załadować wszystkie cykle do jednej struktury danych
- [ ] Zweryfikować kompletność danych (wszystkie 43 runy)

**Skrypt**: `scripts/paper2/load_phase2b_cycles.py`
```python
# Pseudokod
cycles = []
for scenario in ["miller_urey", "hydrothermal", "formamide"]:
    for run_id in range(1, max_runs[scenario] + 1):
        cycles_file = f"results/phase2b_additional/{scenario}_extended/run_{run_id}/autocatalytic_cycles.json"
        if exists(cycles_file):
            cycles.extend(load_json(cycles_file))
```

**Output**: `analysis/phase2b/all_cycles_raw.json`

---

### Krok 2: Klasyfikacja Cykli (Tydzień 1-2)

**Zadania**:
- [ ] Zaimplementować klasyfikację direct autocatalysis (A + B → 2A)
- [ ] Zaimplementować klasyfikację indirect cycles (A → B → C → A)
- [ ] Zaimplementować klasyfikację hypercycles (mutual catalysis)
- [ ] Zweryfikować klasyfikację (manual check top 100)

**Skrypt**: `scripts/paper2/classify_cycles.py`
```python
# Pseudokod
def classify_cycle(cycle):
    if is_direct_autocatalysis(cycle):
        return "direct"
    elif is_hypercycle(cycle):
        return "hypercycle"
    else:
        return "indirect"

classified = [classify_cycle(c) for c in cycles]
```

**Output**: `analysis/phase2b/autocatalytic_cycles_classified.json`

**Struktura Output**:
```json
{
  "total_cycles": 769315,
  "direct": {
    "count": 123456,
    "percentage": 16.0,
    "cycles": [...]
  },
  "indirect": {
    "count": 567890,
    "percentage": 73.8,
    "cycles": [...]
  },
  "hypercycles": {
    "count": 77969,
    "percentage": 10.1,
    "cycles": [...]
  }
}
```

---

### Krok 3: Amplification Factors (Tydzień 2)

**Zadania**:
- [ ] Zaimplementować obliczanie amplification factors
- [ ] Analiza temporalna (d[A]/dt przed vs po)
- [ ] Agregacja per scenario

**Skrypt**: `scripts/paper2/calculate_amplification_factors.py`
```python
# Pseudokod
def calculate_amplification_factor(cycle, snapshots):
    # Znajdź uczestników cyklu
    participants = get_cycle_participants(cycle)
    
    # Oblicz d[A]/dt przed aktywacją cyklu
    rate_before = calculate_rate(participants, snapshots, before_cycle_start)
    
    # Oblicz d[A]/dt po aktywacji cyklu
    rate_after = calculate_rate(participants, snapshots, after_cycle_start)
    
    # Amplification factor
    return rate_after / rate_before if rate_before > 0 else 0
```

**Output**: `analysis/phase2b/amplification_factors.json`

**Struktura Output**:
```json
{
  "miller_urey": {
    "mean": 2.34,
    "std": 0.56,
    "range": [1.05, 8.92],
    "factors": [...]
  },
  "hydrothermal": {
    "mean": 3.12,
    "std": 0.78,
    "range": [1.02, 12.45],
    "factors": [...]
  },
  "formamide": {
    "mean": 2.89,
    "std": 0.67,
    "range": [1.08, 9.34],
    "factors": [...]
  }
}
```

---

### Krok 4: Network Topology Analysis (Tydzień 3)

**Zadania**:
- [ ] Zbudować reaction networks (molecules as nodes, reactions as edges)
- [ ] Obliczyć degree distribution
- [ ] Zidentyfikować hub molecules (top 10% by degree)
- [ ] Obliczyć centrality metrics (betweenness, closeness, eigenvector)
- [ ] Analiza network motifs (3-node, 4-node)

**Skrypt**: `scripts/paper2/network_topology_analysis.py`
```python
# Pseudokod
import networkx as nx

# Zbuduj graf
G = build_reaction_network(all_reactions)

# Degree distribution
degree_dist = dict(G.degree())

# Hub molecules (top 10%)
hubs = sorted(degree_dist.items(), key=lambda x: x[1], reverse=True)[:int(len(G) * 0.1)]

# Centrality metrics
betweenness = nx.betweenness_centrality(G)
closeness = nx.closeness_centrality(G)
eigenvector = nx.eigenvector_centrality(G)

# Network motifs (użyj Mfinder lub NetworkX)
motifs = analyze_motifs(G, size=3)  # 3-node motifs
motifs_4 = analyze_motifs(G, size=4)  # 4-node motifs
```

**Output**: `analysis/phase2b/network_topology_metrics.json`

---

### Krok 5: Scenario Comparison (Tydzień 3-4)

**Zadania**:
- [ ] ANOVA dla amplification factors
- [ ] Chi-square dla typów cykli
- [ ] T-tests dla topology metrics
- [ ] Hub molecule overlap analysis

**Skrypt**: `scripts/paper2/scenario_comparison.py`
```python
# Pseudokod
from scipy import stats

# ANOVA dla amplification factors
f_stat, p_value = stats.f_oneway(
    miller_urey_factors,
    hydrothermal_factors,
    formamide_factors
)

# Chi-square dla typów cykli
chi2, p_value = stats.chi2_contingency(cycle_type_contingency_table)

# T-tests dla topology metrics
t_stat, p_value = stats.ttest_ind(miller_urey_metrics, hydrothermal_metrics)
```

**Output**: `analysis/phase2b/scenario_comparison.json`

---

## 📊 Generowanie Figur

### Figure 1: Cycle Type Distribution
**Skrypt**: `scripts/paper2/generate_figure1_cycle_types.py`
- Panel A: Pie chart (direct, indirect, hypercycles)
- Panel B: Examples of each type (reaction diagrams)
- Panel C: Cycle length distribution (histogram)

### Figure 2: Amplification Factors
**Skrypt**: `scripts/paper2/generate_figure2_amplification.py`
- Panel A: Violin plots by scenario
- Panel B: Amplification vs cycle length (scatter)
- Panel C: Temporal evolution of amplification (time series)

### Figure 3: Network Topology
**Skrypt**: `scripts/paper2/generate_figure3_topology.py`
- Panel A: Degree distribution (log-log plot, all scenarios)
- Panel B: Hub molecules (network visualization, top 20)
- Panel C: Network motifs (frequency comparison)

### Figure 4: Scenario Comparison
**Skrypt**: `scripts/paper2/generate_figure4_comparison.py`
- Panel A: Cycle frequency comparison (bar chart)
- Panel B: Amplification factors comparison (box plots)
- Panel C: Hub molecule overlap (Venn diagram)
- Panel D: Topology metrics comparison (heatmap)

---

## ✅ Checklist Analizy

### Przygotowanie (Tydzień 1)
- [ ] Zweryfikować dostępność wszystkich danych
- [ ] Utworzyć strukturę katalogów `analysis/phase2b/`
- [ ] Napisać skrypty do ładowania danych

### Klasyfikacja (Tydzień 1-2)
- [ ] Zaimplementować klasyfikację cykli
- [ ] Zweryfikować klasyfikację (manual check)
- [ ] Wygenerować statystyki

### Amplification Factors (Tydzień 2)
- [ ] Zaimplementować obliczanie amplification factors
- [ ] Agregacja per scenario
- [ ] Statystyki (mean, std, range)

### Network Topology (Tydzień 3)
- [ ] Zbudować reaction networks
- [ ] Obliczyć wszystkie metryki
- [ ] Zidentyfikować hub molecules

### Scenario Comparison (Tydzień 3-4)
- [ ] Wykonać wszystkie testy statystyczne
- [ ] Wygenerować porównania
- [ ] Przygotować tabele

### Generowanie Figur (Tydzień 4)
- [ ] Figure 1: Cycle types
- [ ] Figure 2: Amplification factors
- [ ] Figure 3: Network topology
- [ ] Figure 4: Scenario comparison

---

## 📅 Timeline

### Styczeń 2026 (Tydzień 1-4)
- **Tydzień 1**: Przygotowanie danych + Klasyfikacja cykli
- **Tydzień 2**: Amplification factors
- **Tydzień 3**: Network topology analysis
- **Tydzień 4**: Scenario comparison + Generowanie figur

### Luty 2026 (Tydzień 5-8)
- **Tydzień 5-6**: Writing Methods & Results
- **Tydzień 7-8**: Writing Introduction & Discussion

---

**Last Updated**: 2025-12-04  
**Status**: Plan gotowy, gotowy do rozpoczęcia analizy (styczeń 2026)

