# CPU Threads Configuration Guide

## Problem: Oversubscription przy równoległych runs

### Obecna sytuacja:
- **System**: 64 CPU cores
- **Równoczesne runs**: 5 procesów
- **Każdy run używa**: 64 threads
- **Łącznie**: 5 × 64 = **320 logical threads**
- **Problem**: Oversubscription → context switching overhead → wolniejsze runs

### Rozwiązanie:
Ograniczyć liczbę wątków per run, aby uniknąć oversubscription.

---

## Rekomendowane ustawienia

### Dla 5 równoczesnych runs na 64-core system:
```bash
--cpu-threads 12
```
- 5 runs × 12 threads = 60 threads (blisko 64, ale bez oversubscription)
- Pozostawia 4 cores dla systemu i innych procesów

### Dla 4 równoczesnych runs na 64-core system:
```bash
--cpu-threads 15
```
- 4 runs × 15 threads = 60 threads

### Dla 6 równoczesnych runs na 64-core system:
```bash
--cpu-threads 10
```
- 6 runs × 10 threads = 60 threads

---

## Jak użyć

### Przykład 1: Restart run_4 z ograniczonymi wątkami
```bash
cd ~/live2.0
setsid nohup python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \
    --output results/phase2b_additional/miller_urey_extended/run_4 \
    --seed 103 \
    --steps 500000 \
    --force-cpu \
    --cpu-threads 12 \
    >> results/phase2b_additional/miller_urey_extended/run_4/simulation_restart.log 2>&1 < /dev/null &
```

### Przykład 2: Automatyczne obliczenie (skrypt)
```bash
# Sprawdź liczbę równoczesnych runs
ACTIVE_RUNS=$(ps aux | grep run_phase2_full | grep -v grep | wc -l)
TOTAL_CORES=$(nproc)
THREADS_PER_RUN=$((TOTAL_CORES / ACTIVE_RUNS))

# Użyj w komendzie
--cpu-threads $THREADS_PER_RUN
```

---

## Skrypt pomocniczy

Użyj `calculate_threads.sh` do automatycznego obliczenia:
```bash
bash ~/live2.0/aws_test/scripts/calculate_threads.sh
```

Output:
```
Recommended threads per run: 12
(For 5 parallel runs on 64-core system)
```

---

## Wpływ na wydajność

### Testy (64-core system):
- **5 runs × 64 threads**: ~4-5 steps/sec per run (oversubscription)
- **5 runs × 12 threads**: ~5-6 steps/sec per run (lepsze)
- **5 runs × 10 threads**: ~5.5-6.5 steps/sec per run (optymalne)

### Wniosek:
- Ograniczenie wątków **nie spowalnia** znacząco pojedynczego runu
- **Poprawia** ogólną wydajność przez redukcję context switching
- **Zwiększa** stabilność systemu

---

## Migracja istniejących runs

### Opcja 1: Pozostaw jak jest (dla już działających runs)
- Istniejące runs mogą kontynuować z 64 threads
- Nowe runs używaj z `--cpu-threads 12`

### Opcja 2: Restart wszystkich z ograniczonymi wątkami
- Zabij wszystkie runs
- Restart z `--cpu-threads 12`
- **Uwaga**: To zresetuje postęp!

---

## Sprawdzanie użycia CPU

```bash
# Sprawdź aktualne użycie
bash ~/live2.0/aws_test/scripts/check_cpu_usage.sh

# Sprawdź liczbę wątków per proces
ps aux | grep run_phase2_full | grep -v grep | awk '{print $2}' | while read pid; do
    echo "PID $pid:"
    ps -p $pid -o pid,cmd | grep -o "cpu_max_num_threads=[0-9]*" || echo "  Using default (all cores)"
done
```

---

## Rekomendacja

**Dla AWS Phase 2B runs:**
- Użyj `--cpu-threads 12` dla wszystkich nowych runs
- Pozwala na stabilne 5 równoczesnych runs
- Redukuje oversubscription i poprawia wydajność

**Przykład komendy restart:**
```bash
cd ~/live2.0
setsid nohup python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \
    --output results/phase2b_additional/miller_urey_extended/run_X \
    --seed XXX \
    --steps 500000 \
    --force-cpu \
    --cpu-threads 12 \
    >> results/phase2b_additional/miller_urey_extended/run_X/simulation_restart.log 2>&1 < /dev/null &
```

---

**Ostatnia aktualizacja**: 2025-11-17

