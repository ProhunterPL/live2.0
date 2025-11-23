# Instrukcje dla AWS - Final Solution

## Problem
Skrypty Phase 2B oczekują `scripts/run_phase2_full.py`, ale ten plik nie jest w repo na AWS.

## Rozwiązanie - Uruchom Ręcznie przez scripts/

```bash
cd ~/live2.0

# Sprawdź czy masz skrypty
ls scripts/ | grep run_phase2

# Jeśli masz scripts/run_phase2_full.py, to uruchom:
python3 scripts/run_phase2_full.py \
  --config aws_test/configs/phase2_miller_urey_test.yaml \
  --output aws_test/results/test_1 \
  --steps 10000 \
  --seed 42

# Sprawdź czy symulacja działa
ls aws_test/results/test_1/

# Jeśli działa, to uruchom batch:
cd scripts
python3 run_phase2_batch.py --all --runs 2 --output ../aws_test/results/phase2_batch
```

## Alternatywnie - Sprawdź co masz w backend/

```bash
cd ~/live2.0

# Szukaj skryptów symulacyjnych
find . -name "*.py" | grep -i "sim\|phase2" | grep -v __pycache__

# Sprawdź backend
ls backend/sim/*.py
```

## Jeśli nie ma skryptów - Musimy użyć backend bezpośrednio

Jeśli nie masz `run_phase2_full.py`, to musisz stworzyć prosty wrapper:

```bash
cd ~/live2.0

# Stwórz prosty test script
cat > test_phase2.py << 'EOF'
import sys
sys.path.insert(0, '.')

from backend.sim.phase2_config import load_phase2_config
import taichi as ti

# Load config
config = load_phase2_config('aws_test/configs/phase2_miller_urey_test.yaml')
print(f"Loaded scenario: {config.scenario_name}")
print("Config OK!")
EOF

python3 test_phase2.py
```

## Podsumowanie - Co teraz zrobić

**Najpierw** - sprawdź strukturę na AWS:

```bash
cd ~/live2.0
ls -la
ls scripts/ 2>/dev/null
ls backend/sim/ 2>/dev/null
ls aws_test/configs/
```

**Potem** - daj mi znać co widzisz, a przygotuję właściwe komendy.

