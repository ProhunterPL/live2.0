# Phase 2B Fix - Podsumowanie

## Problem

1. **Zbyt wiele równoległych symulacji** - uruchomionych było 13 procesów jednocześnie (powinno być max 4)
2. **Zawieszone hydrothermal_extended symulacje** - wszystkie runs 1-10 zatrzymane od ~27 godzin
3. **Brak limitu równoległości** - ubijanie procesów powodowało włączanie kolejnych
4. **Hydrothermal wymagał fixa** - cluster_check_interval powinien być wyłączony (999999999)

## Rozwiązanie

### 1. Skrypty do zarządzania symulacjami

#### `kill_stuck_hydrothermal.sh`
Zatrzymuje zawieszone hydrothermal_extended symulacje (runs 1-10), które nie mają aktywności w logach przez >60 minut.

```bash
bash aws_test/scripts/kill_stuck_hydrothermal.sh
```

#### `limit_parallel_simulations.sh`
Ogranicza liczbę równoległych symulacji do maksymalnie 4 (domyślnie).

```bash
bash aws_test/scripts/limit_parallel_simulations.sh 4
```

#### `restart_from_checkpoint.sh`
Restartuje zawieszone symulacje od najnowszego checkpointu/snapshotu.

```bash
bash aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional hydrothermal_extended
```

**UWAGA**: Obecnie restartuje od początku (checkpoint loading nie jest jeszcze zaimplementowany w `run_phase2_full.py`), ale checkpointy są zachowane dla przyszłej implementacji.

#### `fix_phase2b_issues.sh`
Kompleksowy skrypt wykonujący wszystkie powyższe kroki:

```bash
bash aws_test/scripts/fix_phase2b_issues.sh
```

### 2. Fix dla Hydrothermal

Zastosowano fix z testów lokalnych:
- **Plik**: `aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml`
- **Zmiana**: `cluster_check_interval: 300` → `cluster_check_interval: 999999999` (wyłączone)
- **Efekt**: Zapobiega deadlockowi w cluster detection

### 3. Aktualizacja `run_phase2b_additional.py`

Dodano mechanizm sprawdzający liczbę działających symulacji przed uruchomieniem nowych:
- Sprawdza ile symulacji już działa
- Jeśli >= max_parallel, czeka 30 sekund i sprawdza ponownie
- Jeśli nadal za dużo, wyświetla ostrzeżenie i nie uruchamia nowych

## Instrukcje użycia

### Krok 1: Zatrzymaj zawieszone symulacje

```bash
cd ~/live2.0
bash aws_test/scripts/kill_stuck_hydrothermal.sh
```

### Krok 2: Ogranicz równoległe symulacje do 4

```bash
bash aws_test/scripts/limit_parallel_simulations.sh 4
```

### Krok 3: Sprawdź status

```bash
python3 aws_test/scripts/check_actual_progress.py
```

### Krok 4: (Opcjonalnie) Restartuj zawieszone symulacje

```bash
bash aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional hydrothermal_extended
```

**LUB** użyj kompleksowego skryptu:

```bash
bash aws_test/scripts/fix_phase2b_issues.sh
```

## Monitorowanie

### Sprawdź liczbę działających symulacji

```bash
ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l
```

Powinno pokazać <= 4.

### Sprawdź postęp

```bash
python3 aws_test/scripts/check_actual_progress.py
```

### Sprawdź czy hydrothermal fix jest zastosowany

```bash
grep "cluster_check_interval" aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml
```

Powinno pokazać: `cluster_check_interval: 999999999`

## Uwagi

1. **Checkpoint loading**: Obecnie `run_phase2_full.py` nie obsługuje restartu od checkpointu. Checkpointy są zapisywane, ale restart zawsze zaczyna od początku. To wymaga przyszłej implementacji.

2. **Automatyczne ograniczenie**: `run_phase2b_additional.py` teraz sprawdza liczbę działających symulacji, ale jeśli procesy są uruchamiane ręcznie (np. przez `nohup`), użyj `limit_parallel_simulations.sh`.

3. **Hydrothermal fix**: Fix jest już zastosowany w konfiguracji. Nowe symulacje hydrothermal będą używać poprawionej konfiguracji automatycznie.

## Następne kroki

1. Zatrzymaj zawieszone hydrothermal symulacje
2. Ogranicz równoległe symulacje do 4
3. Monitoruj postęp
4. Jeśli potrzebne, restartuj zawieszone symulacje (będą startować od początku, ale z poprawioną konfiguracją)

