# 🔄 Jak Wznowić Zatrzymane Symulacje Phase 2B

Symulacje się zatrzymały ~5 godzin temu. Oto jak je wznowić:

## 🔍 Krok 1: Sprawdź Dlaczego Się Zatrzymały

```bash
cd ~/live2.0
python3 aws_test/scripts/check_why_stopped.py
```

To pokaże:
- Szczegóły błędów w logach
- Status z phase2b_results.json
- Ostatnie linie logów

---

## 🔄 Krok 2: Wznowienie Symulacji

### Opcja A: Restart Zatrzymanych Runów (Zachowaj Postęp)

Symulacje Phase 2 nie mają wbudowanego resume, ale możesz:

1. **Sprawdź czy są checkpointy** (zapisane stany):
```bash
find ~/live2.0/results/phase2b_additional -name "checkpoint_*.pkl" -o -name "*.checkpoint"
```

2. **Jeśli NIE MA checkpointów** - musisz uruchomić od nowa:
```bash
cd ~/live2.0
python3 aws_test/scripts/run_phase2b_additional.py --scenario miller_urey_extended
```

### Opcja B: Uruchom Tylko Nieukończone Runy

Najpierw sprawdź które runy się nie ukończyły:
```bash
cd ~/live2.0/results/phase2b_additional
# Sprawdź które runy nie mają results.json
for dir in miller_urey_extended/run_*; do
    if [ ! -f "$dir/results.json" ]; then
        echo "Missing: $dir"
    fi
done
```

Potem uruchom tylko te runy ręcznie lub zmodyfikuj skrypt.

### Opcja C: Uruchom Wszystko Od Nowa (Najprostsze)

```bash
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run
```

**UWAGA**: To uruchomi wszystkie 30 symulacji od nowa, w tym te które już się częściowo wykonały.

---

## 💡 Zalecane Rozwiązanie

### Jeśli Symulacje Się Crashowały (błędy w logach):

1. **Sprawdź błędy**:
```bash
cd ~/live2.0
python3 aws_test/scripts/check_why_stopped.py
```

2. **Napraw problem** (np. brak pamięci, błąd w kodzie)

3. **Uruchom ponownie**:
```bash
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run
```

### Jeśli Symulacje Się Zatrzymały Bez Błędów (proces został zabity):

Możliwe przyczyny:
- **OOM (Out of Memory)** - sprawdź `dmesg | grep -i oom`
- **Timeout** - sprawdź czy były timeouty w logach
- **Proces został zabity ręcznie** - sprawdź historię komend

**Rozwiązanie**: Uruchom ponownie z większym limitem pamięci lub na większej instancji AWS.

---

## 🎯 Szybkie Wznowienie (Jeśli Nie Ma Checkpointów)

Ponieważ symulacje Phase 2 nie mają resume, najlepsze rozwiązanie to:

1. **Zapisz obecne logi** (dla analizy):
```bash
cd ~/live2.0/results/phase2b_additional
tar -czf stopped_simulations_backup_$(date +%Y%m%d).tar.gz miller_urey_extended/run_*/simulation.log
```

2. **Uruchom ponownie**:
```bash
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run
```

3. **Monitoruj** żeby upewnić się że działają:
```bash
# W osobnym terminalu
watch -n 60 'ps aux | grep python | grep -v grep | wc -l'
```

---

## ⚠️ Ważne Uwagi

- **Symulacje Phase 2 nie mają resume** - musisz uruchomić od nowa
- **Postęp jest stracony** - run_1 miał 37.8%, run_2 miał 19%
- **Sprawdź dlaczego się zatrzymały** przed restartem - może się powtórzyć
- **Rozważ użycie screen/tmux** żeby symulacje nie zatrzymały się przy rozłączeniu SSH

---

## 🔧 Użycie screen/tmux (Zalecane)

```bash
# Uruchom w screen
screen -S phase2b
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run

# Odłącz: Ctrl+A, potem D
# Podłącz ponownie: screen -r phase2b
```

