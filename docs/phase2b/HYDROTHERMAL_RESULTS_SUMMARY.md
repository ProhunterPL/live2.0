---
date: 2025-11-25
label: summary
---

# Podsumowanie wyników Phase 2B: Hydrothermal Extended

## 📊 Przegląd

**Scenariusz**: Hydrothermal Extended  
**Liczba runów**: 17  
**Status**: Wszystkie runy zakończone pomyślnie  
**Data analizy**: 2025-11-25

---

## 🔬 Kluczowe wyniki

### Różnorodność cząsteczek

- **Średnia liczba unikalnych cząsteczek na run**: 59.5 ± 7.8
- **Zakres**: 47-78 cząsteczek na run
- **Łączna liczba cząsteczek** (wszystkie runy): ~1,000+ unikalnych formuł

### Metryki złożoności

| Metryka | Wartość średnia | Odchylenie standardowe |
|---------|----------------|------------------------|
| **Shannon Entropy** | 2.76 | ± 0.12 |
| **Species Richness** | 59.5 | ± 7.8 |
| **Evenness** | 0.68 | ± 0.03 |
| **Self-organization Index** | 0.21 | ± 0.01 |

### Interpretacja

- **Wysoka różnorodność**: Shannon entropy ~2.76 wskazuje na znaczną różnorodność cząsteczek
- **Dobra równomierność**: Evenness ~0.68 sugeruje względnie równomierny rozkład cząsteczek
- **Samoorganizacja**: Self-organization index ~0.21 wskazuje na umiarkowaną samoorganizację systemu

---

## 📈 Rozkład cząsteczek per run

| Run | Liczba cząsteczek | Run | Liczba cząsteczek |
|-----|-------------------|-----|-------------------|
| run_1 | 67 | run_10 | 58 |
| run_2 | 57 | run_11 | 69 |
| run_3 | 71 | run_12 | 78 |
| run_4 | 58 | run_13 | 60 |
| run_5 | 47 | run_14 | 54 |
| run_6 | 63 | run_15 | 52 |
| run_7 | 59 | run_16 | 51 |
| run_8 | 59 | run_17 | 52 |
| run_9 | 57 | | |

**Statystyki**:
- Minimum: 47 cząsteczek (run_5)
- Maximum: 78 cząsteczek (run_12)
- Mediana: ~58 cząsteczek

---

## 🔄 Cykle autokatalityczne

- **Wykryte cykle**: 0
- **Cycles per run**: 0.0 ± 0.0

### Problem

Brak wykrytych cykli autokatalitycznych wynika z **braku plików `reaction_network.json`**, które są wymagane przez detektor autokatalityczny. Pliki te nie zostały wygenerowane podczas symulacji.

### Rozwiązanie

**Opcja 1: Generowanie sieci reakcji z snapshotów (REKOMENDOWANE)**

Wygeneruj sieć reakcji przez analizę temporalną snapshotów:

```bash
# Dla pojedynczego runu
python scripts/build_reaction_network_from_snapshots.py \
    --run results/phase2b_additional/hydrothermal_extended/run_1

# Dla wszystkich runów (po implementacji batch script)
python scripts/build_reaction_networks_batch.py \
    --scenario hydrothermal_extended
```

**Jak to działa:**
1. Analizuje snapshoty temporalnie (step_00050000.json, step_00100000.json, ...)
2. Wykrywa zmiany w cząsteczkach między snapshotami
3. Inferuje reakcje: jeśli cząsteczka A znika, a B pojawia się → A → B
4. Buduje graf reakcji
5. Generuje `reaction_network.json`

**Po wygenerowaniu sieci:**
```bash
# Uruchom ponownie analizę autokatalityczną
python scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data
```

**Szczegóły**: Zobacz `docs/phase2b/AUTOCATALYSIS_DETECTION_GUIDE.md` dla pełnego przewodnika.

### Uwagi

- Obecna metoda inferencji reakcji jest uproszczona (A znika, B pojawia się → A → B)
- W rzeczywistości reakcje mogą być wielocząsteczkowe i katalizowane
- Snapshoty co 50K kroków mogą przegapić szybkie reakcje
- System wykazuje oznaki samoorganizacji (index = 0.21), więc cykle mogą istnieć, ale wymagają lepszej detekcji

---

## 📁 Wygenerowane pliki

### Analiza
- `paper/results_data/hydrothermal_extended_analysis.json` - szczegółowa analiza
- `paper/results_data/summary_table.csv` - podsumowanie statystyczne
- `paper/results_data/scenario_comparison.json` - porównanie scenariuszy

### Wykresy
- `paper/figures/figure3_molecular_diversity.png` - różnorodność cząsteczek
- `paper/figures/figure4_reaction_networks.png` - sieci reakcji
- `paper/figures/figure5_autocatalytic_cycles.png` - cykle autokatalityczne
- `paper/figures/figure6_novel_molecules.png` - nowe cząsteczki

### Tabele
- `paper/tables/table5_hub_molecules.csv` - cząsteczki hub
- `paper/tables/tableS1_parameters.tex` - parametry symulacji

---

## 🎯 Wnioski

1. **Sukces ekstrakcji**: Wszystkie 17 runów zostało pomyślnie przetworzonych z wyekstrahowanymi cząsteczkami

2. **Stabilna różnorodność**: Średnia ~60 cząsteczek na run z niskim odchyleniem standardowym (7.8) wskazuje na stabilność procesu

3. **Wysoka jakość danych**: Metryki złożoności wskazują na dobrze zorganizowany system chemiczny

4. **Potencjał do dalszej analizy**: 
   - Analiza temporalna zmian cząsteczek w czasie
   - Identyfikacja reakcji między cząsteczkami
   - Porównanie z innymi scenariuszami (Miller-Urey, Formamide)

---

## 📝 Następne kroki

1. **Analiza temporalna**: Prześledzenie ewolucji cząsteczek w czasie (snapshots)
2. **Identyfikacja reakcji**: Budowa sieci reakcji z danych temporalnych
3. **PubChem matching**: Identyfikacja znanych cząsteczek w bazie PubChem
4. **Porównanie scenariuszy**: Analiza różnic między hydrothermal, Miller-Urey i Formamide

---

**Wygenerowano**: 2025-11-25  
**Skrypt analizy**: `scripts/analyze_phase2b_complete.py`  
**Ekstrakcja cząsteczek**: `scripts/fix_run1_molecules.py`

