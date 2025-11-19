# 🔄 Auto-Restart Guide - Phase 2B Miller-Urey

**Problem**: Masz 4 równoczesne runy (5-8) działające, 12 runów zatrzymanych (2-4, 10-18)

**Rozwiązanie**: Automatyczny restart w kolejce - gdy się 4 skończą, następne 4 startują

---

## 📊 Status Aktualny

| Status | Runy | Postęp | Uwagi |
|--------|------|--------|-------|
| ✅ Completed | run_1 | 500K (100%) | Gotowy |
| 🏃 Running | runs 5-8 | ~336K (67%) | ETA: 8h |
| ⏸️ Stopped | runs 2-4 | 0K (restart od początku) | W kolejce |
| ⏸️ Stopped | runs 10-18 | 0K (restart od początku) | W kolejce |

**WAŻNE**: Checkpointy NIE działają - restart będzie od początku (500K kroków każdy).

---

## 🚀 Jak Uruchomić Auto-Restart

### Metoda 1: Screen (Rekomendowana)

```bash
# SSH na AWS
ssh ubuntu@ip-172-31-0-42

cd ~/live2.0

# Uruchom w screen (możesz odłączyć SSH)
screen -S phase2b_queue

# Uruchom skrypt
chmod +x aws_test/scripts/auto_queue_restart.sh
bash aws_test/scripts/auto_queue_restart.sh

# Odłącz screen: Ctrl+A, potem D
# Podłącz ponownie: screen -r phase2b_queue
```

### Metoda 2: Nohup (Alternatywna)

```bash
cd ~/live2.0

chmod +x aws_test/scripts/auto_queue_restart.sh

# Uruchom w tle
nohup bash aws_test/scripts/auto_queue_restart.sh > logs/auto_restart.log 2>&1 &

# Monitoruj
tail -f logs/auto_restart.log
tail -f logs/phase2b_auto_restart.log
```

---

## 🔍 Co Robi Skrypt

### Główna Logika:

1. **Sprawdza co 5 minut** status symulacji
2. **Gdy miejsce się zwolni** (< 4 runy) → uruchamia następne z kolejki
3. **Kolejka wykonania**:
   - Obecnie: runs 5-8 (działają)
   - Następnie: runs 2, 3, 4, 10
   - Potem: runs 11, 12, 13, 14
   - Na końcu: runs 15, 16, 17, 18

4. **Kończy** gdy wszystkie 18 runów complete

### Przykładowy Output:

```
🔄 Auto Queue Restart - Miller-Urey Extended
Settings:
  Max parallel: 4
  Monitor interval: 300s (5 min)

[2025-11-12 10:41:36] Iteration 1 - Checking status...
[2025-11-12 10:41:36]   Running: 4, Queue: 12
[2025-11-12 10:41:36]   📊 Progress: 1 completed, 4 running, 12 queued
[2025-11-12 10:41:36]   ⏰ Sleeping for 300s...

# ~8 godzin później, runs 5-8 się kończą...

[2025-11-12 18:41:36] Iteration 96 - Checking status...
[2025-11-12 18:41:36]   Running: 0, Queue: 12
[2025-11-12 18:41:36]   🆓 Capacity available: 4 slots
[2025-11-12 18:41:36] 🚀 Starting run_2 (seed 101)...
[2025-11-12 18:41:36]    ✅ Started with PID 12345
[2025-11-12 18:41:39] 🚀 Starting run_3 (seed 102)...
[2025-11-12 18:41:39]    ✅ Started with PID 12346
[2025-11-12 18:41:42] 🚀 Starting run_4 (seed 103)...
[2025-11-12 18:41:42]    ✅ Started with PID 12347
[2025-11-12 18:41:45] 🚀 Starting run_10 (seed 109)...
[2025-11-12 18:41:45]    ✅ Started with PID 12348
[2025-11-12 18:41:45]   ✅ Started 4 new simulations
[2025-11-12 18:41:45]   📊 Progress: 5 completed, 4 running, 8 queued
```

---

## 📋 Monitorowanie

### Sprawdź Status:

```bash
# Quick check
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py

# Procesy
ps aux | grep "run_phase2_full.py" | grep -v grep

# Logi skryptu
tail -f ~/live2.0/logs/phase2b_auto_restart.log
```

### Sprawdź Kolejkę:

```bash
# Ile zostało w kolejce
cat ~/live2.0/logs/phase2b_queue.txt

# Ile linijek (= ile runów)
wc -l ~/live2.0/logs/phase2b_queue.txt
```

---

## ⏰ Timeline Przewidywany

| Czas | Event | Runs Completed | Running |
|------|-------|----------------|---------|
| **T+0h** (teraz) | Start auto-restart | 1 | 5-8 (4) |
| **T+8h** | Runs 5-8 done, start 2-4+10 | 5 | 2-4, 10 (4) |
| **T+16h** | Runs 2-4+10 done, start 11-14 | 9 | 11-14 (4) |
| **T+24h** | Runs 11-14 done, start 15-18 | 13 | 15-18 (4) |
| **T+32h** | Runs 15-18 done | **17** | 0 |

**Total: 32 godziny do ukończenia wszystkich 17 runów** (run_9 nie zostanie restartowany - był stuck)

---

## 🛠️ Zarządzanie Skryptem

### Zatrzymaj Auto-Restart:

**Jeśli w screen:**
```bash
screen -r phase2b_queue
# Naciśnij Ctrl+C
```

**Jeśli w nohup:**
```bash
pkill -f "auto_queue_restart.sh"
```

### Restart Skryptu:

```bash
# Zabij stary
pkill -f "auto_queue_restart.sh"

# Uruchom nowy
screen -S phase2b_queue
bash ~/live2.0/aws_test/scripts/auto_queue_restart.sh
```

### Modyfikuj Kolejkę:

Edytuj plik `~/live2.0/logs/phase2b_queue.txt`:

```bash
nano ~/live2.0/logs/phase2b_queue.txt
```

Format: `run_id:seed` (jedna linia = jeden run)

Przykład:
```
2:101
3:102
10:109
```

---

## 🚨 Troubleshooting

### Problem: Skrypt nie startuje nowych runów

**Sprawdź:**
```bash
# Czy skrypt działa?
ps aux | grep "auto_queue_restart.sh"

# Czy jest miejsce?
ps aux | grep "run_phase2_full.py" | wc -l  # Powinno być < 4

# Czy kolejka ma runy?
cat ~/live2.0/logs/phase2b_queue.txt
```

**Fix:**
- Jeśli brak procesu → restart skryptu
- Jeśli kolejka pusta → skrypt czeka na zakończenie aktualnych
- Jeśli 4 runy działają → czeka aż się zwolni miejsce

### Problem: Run się nie kończy (stuck)

**Sprawdź ostatnie logi:**
```bash
tail -50 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_X/simulation.log
```

**Fix:**
- Jeśli stuck → zabij proces: `pkill -f "run_X"`
- Skrypt automatycznie uruchomi następny z kolejki

### Problem: Za mało RAM/CPU

**Zmniejsz równoległość do 2:**

Edytuj `auto_queue_restart.sh` linię 21:
```bash
MAX_PARALLEL=2  # Było 4
```

Restart skryptu.

---

## ✅ Weryfikacja Poprawności

### Po 32h wszystkie runy powinny być gotowe:

```bash
# Sprawdź ile completed
ls -l ~/live2.0/results/phase2b_additional/miller_urey_extended/*/results.json | wc -l

# Powinno pokazać: 17 (runs 1-8, 10-18; run 9 skip)
```

### Każdy run powinien mieć:

```bash
# Dla każdego run_X:
ls -lh ~/live2.0/results/phase2b_additional/miller_urey_extended/run_X/

# Powinieneś zobaczyć:
# - results.json (~5-50 KB)
# - molecules.json (~50-500 KB)
# - snapshots/ (10 plików)
# - checkpoints/ (4-5 plików)
# - simulation.log (~300-500 KB)
```

---

## 📊 Co Dalej Po Zakończeniu

1. **Ekstraktuj molekuły z run_1** (katalog pusty):
   ```bash
   python3 scripts/fix_run1_molecules.py
   ```

2. **Analizuj wszystkie wyniki**:
   ```bash
   python3 aws_test/scripts/analyze_additional_results.py
   ```

3. **Generuj figurki do publikacji**:
   ```bash
   python3 scripts/generate_publication_figures.py
   ```

---

## 🎯 Podsumowanie

**1 linia do uruchomienia wszystkiego:**

```bash
screen -S phase2b_queue -dm bash ~/live2.0/aws_test/scripts/auto_queue_restart.sh && echo "✅ Auto-restart running in background (screen -r phase2b_queue to attach)"
```

**Monitoruj postęp:**

```bash
watch -n 300 "python3 ~/live2.0/aws_test/scripts/check_actual_progress.py | grep -A 3 'Progress:'"
```

**Odłącz SSH** - wszystko działa w tle w screen! 🚀

