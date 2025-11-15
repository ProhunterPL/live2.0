# 🔍 Debug: Symulacje Się Nie Uruchamiają

## Problem
- `ps aux | grep run_phase2_full.py | grep -v grep` pokazuje zero procesów
- Próby uruchomienia przez `screen` lub `nohup` nie działają
- Symulacje się nie startują

---

## 🔧 Krok 1: Diagnostyka Podstawowa

Uruchom na AWS:

```bash
cd ~/live2.0
bash aws_test/scripts/diagnose_startup_issue.sh
```

Ten skrypt sprawdzi:
- ✅ Czy Python działa
- ✅ Czy skrypty istnieją i są czytelne
- ✅ Czy config istnieje
- ✅ Czy importy działają
- ✅ Czy można uruchomić testową symulację

**Jeśli skrypt pokazuje błędy** → napraw je przed kontynuowaniem.

---

## 🚀 Krok 2: Test Ręcznego Uruchomienia (Foreground)

Uruchom symulację w foreground, żeby zobaczyć błędy:

```bash
cd ~/live2.0

# Test z run_3 (seed 103)
bash aws_test/scripts/start_single_simulation.sh 3 103
```

**Co sprawdzić:**
- Czy pojawiają się błędy importu?
- Czy pojawiają się błędy konfiguracji?
- Czy proces się uruchamia i od razu kończy?
- Czy są błędy Taichi/GPU?

**Jeśli działa** → przejdź do Kroku 3 (background).

**Jeśli nie działa** → zobacz sekcję "Typowe Problemy" poniżej.

---

## 🔄 Krok 3: Uruchomienie w Tle

Jeśli foreground działa, uruchom w tle:

```bash
cd ~/live2.0

# Uruchom run_3 w tle
bash aws_test/scripts/start_simulation_background.sh 3 103
```

**Sprawdź czy działa:**
```bash
# Sprawdź proces
ps aux | grep run_phase2_full.py | grep -v grep

# Sprawdź logi
tail -f ~/live2.0/results/phase2b_additional/miller_urey_extended/run_3/simulation.log
```

---

## 🐛 Typowe Problemy i Rozwiązania

### Problem 1: "Python3 not found"

**Rozwiązanie:**
```bash
# Sprawdź czy Python jest zainstalowany
which python3
python3 --version

# Jeśli nie ma, zainstaluj
sudo apt-get update
sudo apt-get install python3 python3-pip
```

---

### Problem 2: "ModuleNotFoundError: No module named 'taichi'"

**Rozwiązanie:**
```bash
cd ~/live2.0
pip3 install taichi numpy
```

Lub jeśli używasz venv:
```bash
source venv/bin/activate  # lub . venv/bin/activate
pip install taichi numpy
```

---

### Problem 3: "FileNotFoundError: Config file"

**Rozwiązanie:**
```bash
cd ~/live2.0

# Sprawdź czy config istnieje
ls -la aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml

# Jeśli nie ma, sprawdź dostępne configi
ls -la aws_test/configs/*.yaml

# Użyj istniejącego configu lub stwórz nowy
```

---

### Problem 4: "Permission denied"

**Rozwiązanie:**
```bash
cd ~/live2.0

# Sprawdź uprawnienia
ls -la scripts/run_phase2_full.py
ls -la aws_test/scripts/*.sh

# Nadaj uprawnienia wykonywania
chmod +x scripts/run_phase2_full.py
chmod +x aws_test/scripts/*.sh

# Sprawdź uprawnienia katalogu wyników
mkdir -p results/phase2b_additional/miller_urey_extended
chmod -R 755 results/
```

---

### Problem 5: "Process dies immediately"

**Sprawdź logi:**
```bash
# Sprawdź log symulacji
tail -100 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_3/simulation.log

# Sprawdź systemowe logi
dmesg | tail -20
journalctl -n 50
```

**Możliwe przyczyny:**
- Brak pamięci (OOM) → zmniejsz równoległość
- Błąd w kodzie → sprawdź traceback w logach
- Problem z Taichi → spróbuj `--force-cpu`

---

### Problem 6: "Screen/nohup nie działa"

**Rozwiązanie dla screen:**
```bash
# Sprawdź czy screen jest zainstalowany
which screen

# Jeśli nie ma
sudo apt-get install screen

# Użyj screen
screen -S phase2b
# W screen:
cd ~/live2.0
bash aws_test/scripts/start_simulation_background.sh 3 103
# Odłącz: Ctrl+A, potem D
```

**Rozwiązanie dla nohup:**
```bash
# nohup powinien działać zawsze, ale sprawdź:
which nohup

# Jeśli nie działa, użyj bezpośrednio:
cd ~/live2.0
python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \
    --output results/phase2b_additional/miller_urey_extended/run_3 \
    --seed 103 \
    --steps 500000 \
    --force-cpu \
    > results/phase2b_additional/miller_urey_extended/run_3/simulation.log 2>&1 &
```

---

### Problem 7: "ImportError: cannot import name 'SimulationConfig'"

**Rozwiązanie:**
```bash
cd ~/live2.0

# Sprawdź czy backend istnieje
ls -la backend/sim/config.py

# Sprawdź PYTHONPATH
echo $PYTHONPATH

# Ustaw PYTHONPATH jeśli potrzeba
export PYTHONPATH="$HOME/live2.0:$PYTHONPATH"

# Test importu
python3 -c "from backend.sim.config import SimulationConfig; print('OK')"
```

---

## 📋 Checklist Przed Uruchomieniem

Przed uruchomieniem symulacji sprawdź:

- [ ] Python3 działa: `python3 --version`
- [ ] Zależności zainstalowane: `pip3 list | grep taichi`
- [ ] Skrypty istnieją: `ls scripts/run_phase2_full.py`
- [ ] Config istnieje: `ls aws_test/configs/*SUPER_FAST.yaml`
- [ ] Katalog wyników istnieje: `ls results/phase2b_additional/`
- [ ] Uprawnienia OK: `ls -la scripts/run_phase2_full.py`
- [ ] Importy działają: `python3 -c "from backend.sim.config import SimulationConfig"`
- [ ] Testowa symulacja działa (100 kroków)

---

## 🎯 Szybkie Rozwiązanie (Jeśli Wszystko Inne Zawodzi)

```bash
cd ~/live2.0

# 1. Sprawdź środowisko
python3 --version
pip3 list | grep -E "taichi|numpy"

# 2. Zainstaluj zależności jeśli potrzeba
pip3 install --user taichi numpy

# 3. Test minimalny (100 kroków)
python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \
    --output results/test_minimal \
    --seed 999 \
    --steps 100 \
    --force-cpu

# 4. Jeśli działa, uruchom pełną symulację
python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \
    --output results/phase2b_additional/miller_urey_extended/run_3 \
    --seed 103 \
    --steps 500000 \
    --force-cpu \
    >> results/phase2b_additional/miller_urey_extended/run_3/simulation.log 2>&1 &

# 5. Sprawdź czy działa
sleep 5
ps aux | grep run_phase2_full.py | grep -v grep
tail -20 results/phase2b_additional/miller_urey_extended/run_3/simulation.log
```

---

## 📞 Jeśli Nadal Nie Działa

Zbierz informacje diagnostyczne:

```bash
cd ~/live2.0

# Zbierz informacje
{
    echo "=== System Info ==="
    uname -a
    python3 --version
    pip3 list
    
    echo ""
    echo "=== File Check ==="
    ls -la scripts/run_phase2_full.py
    ls -la aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml
    
    echo ""
    echo "=== Import Test ==="
    python3 -c "import sys; sys.path.insert(0, '/home/ubuntu/live2.0'); from backend.sim.config import SimulationConfig; print('OK')" 2>&1
    
    echo ""
    echo "=== Process Check ==="
    ps aux | grep python | grep -v grep
    
} > ~/diagnostic_info.txt 2>&1

# Wyświetl
cat ~/diagnostic_info.txt
```

Wyślij zawartość `~/diagnostic_info.txt` do debugowania.

