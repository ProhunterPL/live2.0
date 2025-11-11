# 🚀 Quick Fix - Zastosuj Teraz

## Sytuacja

Z outputu widzę:
- ✅ **4 procesy działają** (PID 47784-47787) - to są miller_urey runs 5-8
- ❌ **Hydrothermal runs 1-10**: Wszystkie zawieszone (logi nieaktualne od 27+ godzin), procesy już nie działają
- ❌ **Miller-Urey runs 2-4, 9**: Zawieszone (logi nieaktualne od 3-27 godzin), procesy już nie działają
- ✅ **Miller-Urey runs 10-18**: Działają poprawnie

## Rozwiązanie - Krok po Kroku

### Krok 1: Zidentyfikuj działające procesy

```bash
cd ~/live2.0
bash aws_test/scripts/identify_running_processes.sh
```

To pokaże które PID odpowiadają którym runom.

### Krok 2: Zatrzymaj wszystkie zawieszone (jeśli jakieś jeszcze działają)

```bash
bash aws_test/scripts/kill_all_stuck.sh 60
```

To zatrzyma procesy z logami nieaktualnymi >60 minut.

### Krok 3: Ogranicz do 4 równoległych

```bash
bash aws_test/scripts/limit_parallel_simulations.sh 4
```

### Krok 4: Restartuj zawieszone hydrothermal runs

```bash
bash aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional hydrothermal_extended
```

**UWAGA**: To uruchomi wszystkie 10 hydrothermal runs, ale tylko 4 będą działać jednocześnie (reszta będzie czekać).

### Krok 5: (Opcjonalnie) Restartuj zawieszone miller_urey runs 2-4, 9

```bash
bash aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional miller_urey_extended
```

To zrestartuje tylko te, które nie mają działających procesów (runs 2-4, 9).

## LUB - Wszystko w Jednym

```bash
cd ~/live2.0
bash aws_test/scripts/complete_fix_and_restart.sh
```

To wykona wszystkie powyższe kroki automatycznie.

## Po Restarcie

Sprawdź status:

```bash
python3 aws_test/scripts/check_actual_progress.py
```

Powinieneś zobaczyć:
- ✅ Max 4 procesy działające jednocześnie
- ✅ Hydrothermal runs z nowymi logami (z poprawioną konfiguracją)
- ✅ Miller-Urey runs 2-4, 9 z nowymi logami

## Uwagi

1. **Checkpoint loading**: Obecnie restartuje od początku (checkpoint loading nie jest jeszcze zaimplementowany), ale używa poprawionej konfiguracji.

2. **Limit 4**: `run_phase2b_additional.py` teraz sprawdza liczbę działających symulacji, ale jeśli uruchamiasz ręcznie, użyj `limit_parallel_simulations.sh`.

3. **Hydrothermal fix**: Fix jest już w konfiguracji (`cluster_check_interval: 999999999`), więc nowe symulacje będą używać poprawionej wersji.

