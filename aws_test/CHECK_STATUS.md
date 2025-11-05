# 🔍 Sprawdzanie Statusu Phase 2B na AWS

## ⚠️ Podejrzany Szybki Sukces

Skrypt `run_phase2b_master.py --mode run` zwrócił sukces **natychmiast**, co jest podejrzane.  
30 symulacji po 500K kroków każda powinno trwać **kilka godzin**, nie sekundy!

---

## 🔍 Co Sprawdzić Na AWS

### **1. Sprawdź czy symulacje rzeczywiście się uruchomiły:**

```bash
# Na AWS (SSH)
cd ~/live2.0/aws_test/results/phase2b_additional

# Sprawdź logi
cat logs/phase2b_runner.log | tail -50

# Sprawdź czy są procesy uruchomione
ps aux | grep python | grep run_phase2

# Sprawdź czy są jakieś wyniki w run directories
ls -la miller_urey_extended/run_1/
ls -la hydrothermal_extended/run_1/
ls -la formamide_extended/run_1/
```

### **2. Sprawdź czy są pliki results.json:**

```bash
# Na AWS
find . -name "results.json" -type f | head -10
find . -name "simulation.log" -type f | head -10
find . -name "snapshots" -type d | head -10
```

### **3. Sprawdź phase2b_results.json:**

```bash
# Na AWS
cat phase2b_results.json | python3 -m json.tool | head -100
```

Szukaj statusów `"status": "success"` vs `"status": "crashed"` lub `"status": "failed"`.

---

## 💡 Możliwe Scenariusze

### **Scenariusz A: Symulacje Failed, ale skrypt zgłosił sukces**

Jeśli `phase2b_results.json` pokazuje wszystkie runs jako `"crashed"` lub `"failed"`, to:
- Skrypt błędnie zgłosił sukces
- Symulacje się nie uruchomiły lub failed natychmiast
- Trzeba poprawić skrypt żeby lepiej wykrywał błędy

### **Scenariusz B: Symulacje są uruchomione w tle**

Jeśli `ps aux | grep python` pokazuje procesy, to:
- Symulacje mogą działać w tle
- Sprawdź czy są logi z postępem
- Monitoruj procesy i zużycie CPU

### **Scenariusz C: Symulacje się nie uruchomiły**

Jeśli nie ma procesów i nie ma wyników, to:
- Skrypt błędnie zgłosił sukces
- Symulacje się nie uruchomiły
- Sprawdź błędy w logach

---

## 🔧 Diagnostyka

### **Sprawdź logi master script:**

```bash
# Na AWS
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run 2>&1 | tee master_run.log
```

### **Sprawdź czy skrypt rzeczywiście uruchamia symulacje:**

```bash
# Na AWS - sprawdź kod
cd ~/live2.0/aws_test/scripts
cat run_phase2b_additional.py | grep -A 10 "run_single_simulation"
```

### **Uruchom ręcznie jedną symulację:**

```bash
# Na AWS - test jednej symulacji
cd ~/live2.0
python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \
    --output results/test_single_run \
    --steps 10000 \
    --seed 42
```

Jeśli to działa, to problem jest w skrypcie `run_phase2b_additional.py`.

---

## ✅ Co Zrobić Teraz

1. **Sprawdź status na AWS** (użyj komend powyżej)
2. **Sprawdź logi** - czy są błędy?
3. **Sprawdź procesy** - czy symulacje działają?
4. **Sprawdź wyniki** - czy są jakieś pliki results.json?

**Następnie**:
- Jeśli symulacje failed → popraw skrypt i uruchom ponownie
- Jeśli symulacje działają → monitoruj postęp
- Jeśli skrypt błędnie zgłosił sukces → popraw logikę wykrywania sukcesu

---

## 📊 Szybka Weryfikacja

Najszybsze sprawdzenie:

```bash
# Na AWS - wszystko w jednej komendzie
cd ~/live2.0/aws_test/results/phase2b_additional && \
echo "=== Processy ===" && ps aux | grep python | grep -v grep | wc -l && \
echo "=== Pliki results.json ===" && find . -name "results.json" | wc -l && \
echo "=== Logi ===" && tail -20 logs/phase2b_runner.log
```

To pokaże:
- Liczbę działających procesów Python
- Liczbę ukończonych symulacji (results.json)
- Ostatnie wpisy z logów

