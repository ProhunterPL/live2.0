# 📋 Gdzie Są Logi Phase 2B

## 🔍 Problem

Service log (`phase2b_service.log`) jest pusty, ale procesy działają. Logi są zapisywane gdzie indziej!

---

## ✅ Gdzie Są Właściwe Logi

### 1. **Główny Log Runnera** (Najważniejszy):
```bash
tail -f ~/live2.0/results/phase2b_additional/logs/phase2b_runner.log
```

### 2. **Logi Poszczególnych Symulacji**:
```bash
# Sprawdź logi symulacji
tail -f ~/live2.0/results/phase2b_additional/miller_urey_extended/run_*/simulation.log
```

### 3. **Service Log** (może być pusty na początku):
```bash
tail -f ~/live2.0/aws_test/phase2b_service.log
```

### 4. **Systemd Journal**:
```bash
sudo journalctl -u phase2b -f
```

---

## 🎯 Sprawdź Gdzie Są Logi

```bash
# Sprawdź główny log runnera
ls -lh ~/live2.0/results/phase2b_additional/logs/

# Sprawdź ostatnie wpisy
tail -50 ~/live2.0/results/phase2b_additional/logs/phase2b_runner.log

# Sprawdź czy są nowe symulacje
ls -lh ~/live2.0/results/phase2b_additional/miller_urey_extended/run_*/simulation.log | tail -5
```

---

## 📊 Monitorowanie

### Sprawdź Postęp:
```bash
cd ~/live2.0
python3 aws_test/scripts/quick_diagnose.py
```

### Sprawdź Procesy:
```bash
ps aux | grep python | grep run_phase2
```

### Sprawdź Główny Log:
```bash
tail -f ~/live2.0/results/phase2b_additional/logs/phase2b_runner.log
```

---

## ⚠️ Ważne

- **Service log może być pusty** - logi są w `results/phase2b_additional/logs/`
- **Sprawdź główny log** - `phase2b_runner.log` zawiera wszystkie informacje
- **Logi symulacji** - każda symulacja ma swój `simulation.log`

---

## ✅ Szybkie Sprawdzenie

```bash
# Sprawdź czy logi się aktualizują
tail -20 ~/live2.0/results/phase2b_additional/logs/phase2b_runner.log

# Sprawdź postęp
cd ~/live2.0
python3 aws_test/scripts/quick_diagnose.py
```

