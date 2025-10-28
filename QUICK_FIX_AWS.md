# Quick Fix - AWS Phase 2B

## Problem
Script używał `python` zamiast `python3` na AWS Ubuntu.

## Rozwiązanie
Naprawiono w pliku `aws_test/run_phase2b_master.py`.

## Co teraz zrobić na AWS:

### 1. Zaktualizuj plik na AWS

Na AWS, skopiuj poprawiony plik:

```bash
# Na AWS, w katalogu aws_test:
nano run_phase2b_master.py

# Znajdź wszystkie linie z "python" i zmień na "python3"
# Lub przepisz cały plik z poprawek
```

### LUB wykonaj te komendy na AWS:

```bash
cd ~/live2.0/aws_test

# Zmień wszystkie python na python3 w skrypcie
sed -i 's/"python"/"python3"/g' run_phase2b_master.py

# Sprawdź czy są jeszcze wywołania python
grep -n '"python"' run_phase2b_master.py
```

### 2. Sprawdź czy istnieją foldery

```bash
cd ~/live2.0/aws_test

# Sprawdź strukturę
ls -la scripts/

# Jeśli brakuje, stwórz
mkdir -p results/phase2b_additional
```

### 3. Uruchom ponownie

```bash
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode all
```

---

## Jeśli nadal nie działa

Sprawdź czy wszystkie zależności są zainstalowane:

```bash
# Sprawdź python3
which python3
python3 --version

# Sprawdź czy masz pip
python3 -m pip --version

# Zainstaluj zależności jeśli brakuje
cd ~/live2.0
pip3 install -r requirements.txt
```

---

## Alternatywa: Uruchom bezpośrednio skrypty

Jeśli master script nadal nie działa, uruchom skrypty bezpośrednio:

```bash
cd ~/live2.0/aws_test

# 1. Debug formamide
python3 scripts/debug_formamide.py --output-dir results/phase2b_additional/formamide_debug

# 2. Uruchom 30 symulacji
python3 scripts/run_phase2b_additional.py --output-dir results/phase2b_additional

# 3. Monitoruj postęp (w osobnej sesji)
python3 scripts/monitor_runs.py --results-dir results/phase2b_additional

# 4. Analizuj wyniki
python3 scripts/analyze_additional_results.py
```

---

**Status**: Naprawiono lokalnie, trzeba zaktualizować na AWS.

