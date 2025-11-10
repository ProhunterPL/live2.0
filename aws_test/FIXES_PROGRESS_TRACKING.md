# 🔧 Naprawy Problemu z Brakiem Postępu w Symulacjach Phase2B

## 📋 Zidentyfikowane Problemy

1. **Błąd parsowania logów w `quick_diagnose.py`**
   - Skrypt szukał wzorca `"Step (\d+) completed"`, ale logi mają format `"Step X/Y"`
   - Skrypt nie mógł odczytać aktualnego postępu

2. **Brak flush logów w `run_phase2_full.py`**
   - Logi były buforowane i nie zapisywane od razu
   - Postęp nie był widoczny w czasie rzeczywistym

3. **Błąd składni w `run_phase2b_additional.py` (na AWS)**
   - Błąd: `parser.add_argument("--max_parallel=4,` - nieprawidłowa składnia
   - Prawidłowa składnia: `parser.add_argument("--max-parallel", type=int, default=4, ...)`

## ✅ Wprowadzone Naprawy

### 1. Naprawiono `aws_test/scripts/quick_diagnose.py`
- Dodano parsowanie wzorca `"Step X/Y"` (z obsługą przecinków)
- Dodano fallback dla starego formatu `"Step X completed"`

### 2. Naprawiono `scripts/run_phase2_full.py`
- Dodano `FlushingFileHandler` - automatycznie flush po każdym logu
- Dodano ręczny flush po każdym logu postępu (co 10,000 kroków)
- Postęp będzie teraz widoczny w czasie rzeczywistym

### 3. Weryfikacja `aws_test/scripts/run_phase2b_additional.py`
- Lokalna wersja jest poprawna (linia 318)
- Na AWS może być inna wersja z błędem składni

## 🚀 Instrukcje dla AWS

### Krok 1: Zsynchronizuj pliki z repozytorium

```bash
cd ~/live2.0
git pull origin main  # lub odpowiednia gałąź
```

### Krok 2: Sprawdź czy błąd składni został naprawiony

```bash
python3 -m py_compile aws_test/scripts/run_phase2b_additional.py
```

Jeśli nie ma błędów, plik jest poprawny.

### Krok 3: Jeśli błąd nadal występuje, napraw ręcznie

```bash
nano aws_test/scripts/run_phase2b_additional.py
```

Znajdź linię 318 i upewnij się, że wygląda tak:
```python
parser.add_argument("--max-parallel", type=int, default=4,
                   help="Maximum parallel simulations (default: 4 for 64 CPU cores)")
```

**NIE** tak:
```python
parser.add_argument("--max_parallel=4,  # ❌ BŁĄD
```

### Krok 4: Sprawdź postęp z poprawionym skryptem

```bash
python3 aws_test/scripts/quick_diagnose.py
```

Teraz powinien poprawnie odczytywać postęp z logów.

### Krok 5: Dla nowych symulacji - postęp będzie widoczny od razu

Nowe symulacje uruchomione po synchronizacji będą miały:
- ✅ Logi flushowane natychmiast
- ✅ Postęp widoczny w czasie rzeczywistym
- ✅ Poprawne odczytywanie przez `quick_diagnose.py`

### Krok 6: Dla istniejących symulacji

Istniejące symulacje (run_5, run_6, run_7, run_8) będą nadal działać, ale:
- Logi mogą być buforowane (stary kod)
- `quick_diagnose.py` powinien teraz poprawnie odczytywać postęp z istniejących logów

## 🔍 Diagnostyka

### Sprawdź czy symulacje rzeczywiście działają:

```bash
# Sprawdź procesy
ps aux | grep run_phase2_full.py | grep -v grep

# Sprawdź ostatnie logi (powinny być aktualizowane co 10,000 kroków)
tail -20 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_5/simulation.log

# Sprawdź postęp przez poprawiony skrypt
python3 aws_test/scripts/quick_diagnose.py
```

### Jeśli postęp nadal nie jest widoczny:

1. **Sprawdź czy logi są zapisywane:**
   ```bash
   ls -lh ~/live2.0/results/phase2b_additional/miller_urey_extended/run_*/simulation.log
   ```

2. **Sprawdź czy procesy używają CPU:**
   ```bash
   top -p $(pgrep -f run_phase2_full.py | tr '\n' ',' | sed 's/,$//')
   ```

3. **Sprawdź czy symulacje nie są zawieszone:**
   ```bash
   # Sprawdź ostatnie 50 linii logu dla błędów
   tail -50 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_5/simulation.log | grep -i error
   ```

## 📝 Uwagi

- **Buforowanie logów**: Python może buforować logi. Nowy kod wymusza flush po każdym logu.
- **Częstotliwość logowania**: Postęp jest logowany co 10,000 kroków. Dla 500,000 kroków = 50 wpisów.
- **Czas między logami**: Przy ~10-12 kroków/sekundę, logi powinny pojawiać się co ~15-20 minut.

## ✅ Oczekiwane Rezultaty

Po zastosowaniu napraw:
1. ✅ `quick_diagnose.py` poprawnie odczyta postęp z logów
2. ✅ Nowe symulacje będą miały logi flushowane natychmiast
3. ✅ Postęp będzie widoczny w czasie rzeczywistym
4. ✅ Błąd składni w `run_phase2b_additional.py` zostanie naprawiony

