---
date: 2025-11-25
label: guide
---

# Wykrywanie cykli autokatalitycznych - Przewodnik

## Problem

Analiza Phase 2B hydrothermal wykazała **0 wykrytych cykli autokatalitycznych**, mimo że system wykazuje oznaki samoorganizacji (self-organization index = 0.21).

## Przyczyna

Detektor autokatalityczny wymaga pliku `reaction_network.json`, który **nie został wygenerowany podczas symulacji**. Plik ten zawiera:
- Sieć reakcji chemicznych (które cząsteczki reagują z którymi)
- Graf reakcji (NetworkX DiGraph)
- Historię obfitości cząsteczek w czasie

## Rozwiązanie

### Opcja 1: Generowanie sieci reakcji z snapshotów (REKOMENDOWANE)

Mamy dostęp do snapshotów z danymi o wiązaniach. Możemy wygenerować sieć reakcji przez analizę temporalną:

```bash
# Dla pojedynczego runu
python scripts/build_reaction_network_from_snapshots.py \
    --run results/phase2b_additional/hydrothermal_extended/run_1

# Dla wszystkich runów (skrypt batch)
python scripts/build_reaction_networks_batch.py \
    --scenario hydrothermal_extended
```

**Jak to działa:**
1. Analizuje snapshoty temporalnie (step_00050000.json, step_00100000.json, ...)
2. Wykrywa zmiany w cząsteczkach między snapshotami
3. Inferuje reakcje: jeśli cząsteczka A znika, a B pojawia się → A → B
4. Buduje graf reakcji
5. Generuje `reaction_network.json`

### Opcja 2: Analiza temporalna obfitości cząsteczek

Alternatywnie, możemy analizować zmiany obfitości cząsteczek w czasie bezpośrednio:

```python
# Przykład analizy temporalnej
from backend.sim.analysis.autocatalysis_detector import AutocatalysisDetector
import networkx as nx

# Buduj graf z danych temporalnych
reaction_graph = nx.DiGraph()
abundance_history = {}

# Analizuj snapshoty
for snapshot in snapshots:
    molecules = extract_molecules(snapshot)
    for mol in molecules:
        if mol not in abundance_history:
            abundance_history[mol] = []
        abundance_history[mol].append(molecules[mol])
    
    # Inferuj reakcje z zmian
    # ...

# Wykryj cykle
detector = AutocatalysisDetector()
cycles = detector.detect_cycles_in_network(
    reaction_graph, 
    abundance_history, 
    molecule_names
)
```

### Opcja 3: Użycie danych z molecules.json

Mamy wyekstrahowane cząsteczki w `molecules.json`. Możemy:
1. Analizować zmiany w `first_seen` / `last_seen`
2. Budować sieć na podstawie podobieństwa strukturalnego
3. Inferować reakcje z sekwencji pojawiania się

## Implementacja

### Krok 1: Wygeneruj sieci reakcji dla wszystkich runów

```bash
# Stwórz skrypt batch
python scripts/build_reaction_networks_batch.py \
    --scenario hydrothermal_extended \
    --base-dir results/phase2b_additional
```

### Krok 2: Uruchom detektor autokatalityczny

```bash
# Dla pojedynczego runu
python scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data
```

Lub bezpośrednio:

```python
from backend.sim.analysis.autocatalysis_detector import analyze_scenario_autocatalysis

results = analyze_scenario_autocatalysis(
    "results/phase2b_additional",
    "hydrothermal_extended",
    list(range(1, 18))
)
```

### Krok 3: Weryfikacja wyników

Sprawdź wygenerowane pliki:
- `run_X/reaction_network.json` - sieć reakcji
- `run_X/autocatalytic_cycles.json` - wykryte cykle
- `paper/results_data/hydrothermal_extended_analysis.json` - zaktualizowana analiza

## Uwagi techniczne

### Ograniczenia obecnej metody

1. **Uproszczona inferencja reakcji**: Obecna metoda inferuje reakcje jako "A znika, B pojawia się → A → B", co jest uproszczeniem. W rzeczywistości:
   - Reakcje mogą być wielocząsteczkowe (A + B → C + D)
   - Mogą być katalizowane przez inne cząsteczki
   - Mogą być odwracalne

2. **Brak informacji o mechanizmach**: Nie mamy informacji o:
   - Energii aktywacji
   - Kinetyce reakcji
   - Katalizatorach

3. **Temporal resolution**: Snapshoty są co 50K kroków, więc możemy przegapić szybkie reakcje

### Ulepszenia (przyszłe)

1. **Lepsza inferencja reakcji**:
   - Analiza zmian strukturalnych (które wiązania się tworzą/rozpadają)
   - Matching z bazą reakcji znanych z literatury
   - Machine learning do przewidywania reakcji

2. **Analiza temporalna z wyższą rozdzielczością**:
   - Analiza checkpointów (co 100K kroków)
   - Interpolacja między snapshotami

3. **Integracja z symulacją**:
   - Zapisywanie reakcji podczas symulacji
   - Tracking zmian w czasie rzeczywistym

## Przykładowe wyniki

Po wygenerowaniu sieci reakcji, oczekiwane wyniki:

```json
{
  "autocatalysis": {
    "total_cycles": 5,
    "cycles_per_run_mean": 0.29,
    "cycle_types": {
      "direct": 2,
      "indirect": 3
    }
  }
}
```

## Referencje

- `backend/sim/analysis/autocatalysis_detector.py` - główny detektor
- `scripts/build_reaction_network_from_snapshots.py` - generator sieci
- `scripts/analyze_phase2b_complete.py` - pipeline analizy

---

**Status**: W trakcie implementacji  
**Ostatnia aktualizacja**: 2025-11-25

