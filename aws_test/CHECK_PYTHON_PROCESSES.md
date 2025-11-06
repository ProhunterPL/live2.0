# 🔍 Sprawdzanie Procesów Python na AWS

Output z `top` pokazuje tylko procesy systemowe na początku listy. Procesy Python mogą być niżej w liście (top sortuje po CPU użyciu).

## ✅ Polecenia do Sprawdzenia Procesów Python

### 1. Sprawdź wszystkie procesy Python:
```bash
ps aux | grep python | grep -v grep
```

### 2. Sprawdź tylko procesy związane z symulacjami:
```bash
ps aux | grep -E "run_phase2|run_phase2b" | grep -v grep
```

### 3. Sprawdź procesy Python z użyciem CPU (top dla Python):
```bash
top -b -n 1 | grep python
```

### 4. Sprawdź ile procesów Python działa:
```bash
ps aux | grep python | grep -v grep | wc -l
```

### 5. Sprawdź szczegóły procesów Python (z czasem działania):
```bash
ps aux | grep python | grep -v grep | awk '{print $2, $9, $10, $11, $12, $13}'
```

### 6. Sprawdź procesy Python w top (interaktywny):
```bash
# W top, naciśnij:
# 'c' - pokaż pełną komendę
# '/' - wyszukaj "python"
# 'o' - filtruj po nazwie procesu (wpisz: COMMAND=python)
```

---

## 🔍 Sprawdź Czy Symulacje Działają (Alternatywa)

Jeśli nie widzisz procesów Python w `top`, sprawdź czy symulacje rzeczywiście działają:

### 1. Sprawdź logi - czy są aktualizowane?
```bash
cd ~/live2.0/results/phase2b_additional
find . -name "simulation.log" -exec ls -lh {} \; | head -5
```

### 2. Sprawdź ostatnie kroki w logach:
```bash
tail -5 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_*/simulation.log | grep "Step"
```

### 3. Sprawdź kiedy ostatnio były aktualizowane logi:
```bash
find ~/live2.0/results/phase2b_additional -name "simulation.log" -mmin -10
```
(Pokazuje logi zaktualizowane w ostatnich 10 minutach)

### 4. Użyj skryptu do sprawdzenia postępu:
```bash
cd ~/live2.0
python3 aws_test/scripts/check_phase2b_progress.py --results-dir results/phase2b_additional
```

---

## 💡 Możliwe Powody Braku Procesów w Top

1. **Procesy są niżej w liście** - `top` sortuje po CPU, jeśli Python nie używa dużo CPU w danej chwili, będzie niżej
2. **Procesy działają w tle** - mogą być uruchomione jako background jobs (`jobs`)
3. **Procesy działają w screen/tmux** - sprawdź `screen -ls` lub `tmux ls`
4. **Procesy się zakończyły** - sprawdź logi czy są błędy

---

## 🎯 Najlepsze Polecenie do Szybkiego Sprawdzenia

```bash
# Wszystko w jednej komendzie:
cd ~/live2.0/results/phase2b_additional && \
echo "=== Procesy Python ===" && \
ps aux | grep python | grep -v grep && \
echo "" && \
echo "=== Liczba procesów ===" && \
ps aux | grep python | grep -v grep | wc -l && \
echo "" && \
echo "=== Ostatnie kroki w logach ===" && \
tail -3 miller_urey_extended/run_*/simulation.log 2>/dev/null | grep "Step" | tail -5
```

