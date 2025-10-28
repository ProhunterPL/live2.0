# Sprawdzenie Statusu na AWS

## Kroki diagnostyczne:

```bash
# 1. Sprawdź strukturę projektu
cd ~/live2.0
pwd
ls -la

# 2. Sprawdź aws_test
ls -la aws_test/
ls -la aws_test/scripts/
ls -la aws_test/configs/

# 3. Sprawdź czy są już jakieś wyniki
find . -name "phase2b_additional" -type d 2>/dev/null
find . -name "results" -type d 2>/dev/null

# 4. Sprawdź backend
ls -la backend/
ls -la backend/sim/

# 5. Sprawdź czy skrypty run_phase2_full istnieją
find . -name "run_phase2_full.py" -o -name "*.py" | grep phase2
```

## Uruchomienie poprawne:

```bash
cd ~/live2.0/aws_test

# Zobacz dokładnie co skrypty próbują uruchomić
cat scripts/debug_formamide.py | grep -A5 "scripts/run_phase2_full"

# Sprawdź czy ten plik istnieje
find ~/live2.0 -name "run_phase2_full.py"

# Jeśli nie istnieje, szukaj podobnych
find ~/live2.0 -name "*phase2*.py"
```

