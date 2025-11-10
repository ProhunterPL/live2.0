# 📊 Status Wyników Phase 2B Additional

**Data analizy**: 2025-11-10  
**Katalog**: `results/phase2b_additional`

## ✅ Zakończone Symulacje

### run_1 ✅ COMPLETED
- **Status**: Zakończona pomyślnie
- **Pliki**: 
  - ✅ `results.json` - istnieje
  - ✅ `molecules.json` - istnieje
  - ✅ `summary.txt` - istnieje
- **Postęp**: 500,000/500,000 kroków (100%)
- **Checkpoints**: 4 (100K, 200K, 300K, 400K)
- **Snapshots**: 10 (co 50K kroków)
- **Seed**: 100
- **Czas zakończenia**: 2025-11-09 11:48:44

## 🔄 Symulacje w Trakcie (z logów)

### run_2 ⏸️ STOPPED
- **Status**: Zatrzymana (brak procesu)
- **Ostatni log**: Step 50,000 (snapshot)
- **Checkpoints**: Brak
- **Snapshots**: 1 (50K)
- **Postęp**: ~50,000/500,000 (10%)

### run_3 ⏸️ STOPPED
- **Status**: Zatrzymana (brak procesu)
- **Ostatni log**: Step 300,000 (checkpoint)
- **Checkpoints**: 3 (100K, 200K, 300K)
- **Snapshots**: 6 (do 300K)
- **Postęp**: ~300,000/500,000 (60%)

### run_4 ⏸️ STOPPED
- **Status**: Zatrzymana (brak procesu)
- **Ostatni log**: Step 100,000 (checkpoint)
- **Checkpoints**: 1 (100K)
- **Snapshots**: 3 (do 150K)
- **Postęp**: ~100,000/500,000 (20%)

### run_5 🔄 RUNNING (na serwerze)
- **Status**: Działa na serwerze (PID 28017)
- **Ostatni log w pliku**: Step 439,000 (2025-11-10 02:30:56)
- **Checkpoints**: 4 (100K, 200K, 300K, 400K)
- **Snapshots**: 9 (do 450K)
- **Postęp**: ~439,000/500,000 (87.8%)
- **Uwaga**: Logi są buforowane - rzeczywisty postęp może być wyższy

### run_6 🔄 RUNNING (na serwerze)
- **Status**: Działa na serwerze (PID 28426)
- **Ostatni log w pliku**: Step 78,000 (2025-11-09 23:03:15)
- **Checkpoints**: Brak
- **Snapshots**: 1 (50K)
- **Postęp**: ~78,000/500,000 (15.6%)
- **Uwaga**: Logi są buforowane - rzeczywisty postęp może być wyższy

### run_7 🔄 RUNNING (na serwerze)
- **Status**: Działa na serwerze (PID 28427)
- **Ostatni log w pliku**: Step 363,000 (2025-11-10 07:26:37)
- **Checkpoints**: 3 (100K, 200K, 300K)
- **Snapshots**: 7 (do 350K)
- **Postęp**: ~363,000/500,000 (72.6%)
- **Uwaga**: Logi są buforowane - rzeczywisty postęp może być wyższy

### run_8 🔄 RUNNING (na serwerze)
- **Status**: Działa na serwerze (PID 28428)
- **Ostatni log w pliku**: Step 104,000 (2025-11-10 00:08:34)
- **Checkpoints**: 1 (100K)
- **Snapshots**: 2 (50K, 100K)
- **Postęp**: ~104,000/500,000 (20.8%)
- **Uwaga**: Logi są buforowane - rzeczywisty postęp może być wyższy

## 📊 Podsumowanie

### Status ogólny (z phase2b_results.json):
- **Total runs**: 30
- **Completed**: 1/30 (3.3%)
- **Failed**: 3/30 (10.0%)
- **W trakcie**: 4 (run_5, run_6, run_7, run_8)
- **Zatrzymane**: 3 (run_2, run_3, run_4)

### Postęp według ostatnich logów:
| Run | Ostatni Krok | Postęp | Status |
|-----|--------------|--------|--------|
| run_1 | 500,000 | 100% | ✅ Zakończona |
| run_2 | ~50,000 | 10% | ⏸️ Zatrzymana |
| run_3 | ~300,000 | 60% | ⏸️ Zatrzymana |
| run_4 | ~100,000 | 20% | ⏸️ Zatrzymana |
| run_5 | 439,000 | 87.8% | 🔄 Działa |
| run_6 | 78,000 | 15.6% | 🔄 Działa |
| run_7 | 363,000 | 72.6% | 🔄 Działa |
| run_8 | 104,000 | 20.8% | 🔄 Działa |

## ⚠️ Uwagi

1. **Buforowanie logów**: Symulacje run_5, run_6, run_7, run_8 działają na serwerze, ale logi są buforowane. Rzeczywisty postęp może być wyższy niż pokazują logi.

2. **Zatrzymane symulacje**: run_2, run_3, run_4 zostały zatrzymane przed zakończeniem. Można je wznowić z checkpointów.

3. **Brakujące wyniki**: Tylko run_1 ma kompletne wyniki (results.json). Pozostałe są w trakcie lub zatrzymane.

4. **Checkpoints**: Większość symulacji ma checkpoints, które można użyć do wznowienia.

## 🔍 Następne Kroki

1. **Sprawdź na serwerze** czy run_5, run_6, run_7, run_8 są nadal aktywne
2. **Sprawdź czy run_5 i run_7** zakończyły się (mogły zakończyć się po pobraniu katalogu)
3. **Rozważ wznowienie** run_2, run_3, run_4 z checkpointów
4. **Poczekaj na zakończenie** run_6 i run_8

