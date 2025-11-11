# Opcje Monitorowania i Restartu

## Sytuacja Obecna

✅ **4 procesy działają** (miller_urey runs 5-8) - wszystkie aktywne  
⏳ **14 symulacji czeka na restart**:
- 10 hydrothermal_extended runs (1-10)
- 4 miller_urey_extended runs (2-4, 9)

## Opcje

### Opcja 1: Poczekaj aż runs 5-8 się zakończą (Zalecane)

Runs 5-8 są już w toku (krok ~24-26K z 500K). Gdy się zakończą, automatycznie uruchomią się kolejne z kolejki.

**Czas oczekiwania**: ~10-15 godzin (szacunkowo)

**Co zrobić**: Nic - po prostu monitoruj:
```bash
python3 aws_test/scripts/check_actual_progress.py
```

### Opcja 2: Automatyczny restart gdy sloty się zwolnią

Uruchom monitor, który automatycznie restartuje zawieszone symulacje gdy zwolnią się sloty:

```bash
# W tle (użyj screen lub nohup)
screen -S auto_restart
bash aws_test/scripts/auto_restart_when_ready.sh

# Lub z nohup
nohup bash aws_test/scripts/auto_restart_when_ready.sh > auto_restart.log 2>&1 &
```

Monitor będzie:
- Sprawdzał co 5 minut czy są wolne sloty
- Automatycznie restartował zawieszone symulacje
- Respektował limit 4 równoległych

**Zatrzymaj monitor**: `Ctrl+C` lub `kill <PID>`

### Opcja 3: Ręczny restart gdy sloty się zwolnią

Gdy runs 5-8 się zakończą (lub gdy chcesz je zatrzymać), uruchom:

```bash
# Restart hydrothermal
bash aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional hydrothermal_extended

# Restart miller_urey
bash aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional miller_urey_extended
```

### Opcja 4: Zatrzymaj niektóre runs 5-8 żeby zrobić miejsce

Jeśli priorytetem jest hydrothermal, możesz zatrzymać niektóre z runs 5-8:

```bash
# Zatrzymaj run_5 (PID 47784)
kill 47784

# Sprawdź dostępne sloty
ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l

# Restart hydrothermal
bash aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional hydrothermal_extended
```

**UWAGA**: To spowoduje utratę postępu w zatrzymanym runie (będzie musiał zacząć od nowa).

## Rekomendacja

**Opcja 1** (poczekaj) jest najlepsza, ponieważ:
- Runs 5-8 są już w toku i działają poprawnie
- Nie tracisz postępu
- Automatycznie uruchomią się kolejne gdy sloty się zwolnią

**Opcja 2** (automatyczny monitor) jest dobra jeśli chcesz być pewny że restart nastąpi natychmiast gdy sloty się zwolnią.

## Monitorowanie

Sprawdzaj postęp regularnie:

```bash
# Pełny status
python3 aws_test/scripts/check_actual_progress.py

# Szybki status
ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l

# Identyfikacja procesów
bash aws_test/scripts/identify_running_processes.sh
```

