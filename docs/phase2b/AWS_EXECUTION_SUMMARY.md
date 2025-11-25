---
date: 2025-11-25
label: status
---

# Podsumowanie wykonania na AWS

## ✅ Co zostało zrobione

### 1. Generowanie sieci reakcji
- **Status**: ✅ ZAKOŃCZONE
- **Wynik**: 17/17 runów przetworzonych pomyślnie
- **Czas**: ~20 sekund (z równoległością 4 workers)
- **Wygenerowane pliki**: `reaction_network.json` dla każdego runu

### 2. Analiza autokatalityczna
- **Status**: ⏳ W TRAKCIE / ZAKOŃCZONE
- **Skrypt**: `analyze_phase2b_complete.py`
- **Czas**: ~1-2 minuty

## 📊 Wyniki generowania sieci reakcji

| Run | Molecules | Reactions | Edges |
|-----|-----------|-----------|-------|
| run_1 | 27 | 230 | ~230 |
| run_2 | 43 | 768 | ~768 |
| run_3 | 45 | 429 | ~429 |
| run_4 | 39 | 608 | ~608 |
| run_5 | 50 | 248 | ~248 |
| run_6 | 35 | 396 | ~396 |
| run_7 | 36 | 443 | ~443 |
| run_8 | 39 | 512 | ~512 |
| run_9 | 38 | 465 | ~465 |
| run_10 | 41 | 618 | ~618 |
| run_11 | 35 | 442 | ~442 |
| run_12 | 24 | 234 | ~234 |
| run_13 | 35 | 463 | ~463 |
| run_14 | 39 | 486 | ~486 |
| run_15 | 33 | 456 | ~456 |
| run_16 | 44 | 560 | ~560 |
| run_17 | 32 | 362 | ~362 |

**Łącznie**: ~7,500 reakcji wygenerowanych z 170 snapshotów

## 🔄 Następne kroki

### Jeśli analiza jeszcze trwa:
- Poczekaj ~1-2 minuty
- Sprawdź status: `ps aux | grep analyze_phase2b`

### Jeśli analiza zakończona:
1. **Pobierz wyniki** (opcjonalnie):
   ```bash
   scp -i "C:\Users\klawi\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem" \
       -r ubuntu@63.178.224.65:~/live2.0/paper/results_data/* \
       paper/results_data/
   ```

2. **Sprawdź wykryte cykle**:
   ```bash
   ssh -i "..." ubuntu@63.178.224.65 \
       "grep -r 'total_cycles' ~/live2.0/paper/results_data/"
   ```

## ⚡ Szybkie sprawdzenie statusu

```bash
# Sprawdź czy wszystkie sieci wygenerowane
ssh -i "..." ubuntu@63.178.224.65 \
    "find ~/live2.0/results/phase2b_additional/hydrothermal_extended/run_*/reaction_network.json | wc -l"
# Powinno być: 17

# Sprawdź czy analiza zakończona
ssh -i "..." ubuntu@63.178.224.65 \
    "test -f ~/live2.0/paper/results_data/hydrothermal_extended_analysis.json && echo 'Done' || echo 'In progress'"
```

---

**Status**: Sieci reakcji wygenerowane ✅  
**Analiza**: W trakcie / Zakończona ⏳  
**Ostatnia aktualizacja**: 2025-11-25 09:08

