# Uruchom Prawdziwe Symulacje

## Problem
Skrypty Phase 2B to template'y - nie wykonują rzeczywistych symulacji.

## Sprawdź dostępne konfiguracje:

```bash
cd ~/live2.0

# Sprawdź co jest w aws_test/configs/
ls -la aws_test/configs/

# Sprawdź co jest w głównym configs/
ls -la configs/ | grep phase2
```

## Rozwiązanie - Użyj run_phase2_batch.py

```bash
cd ~/live2.0

# Najpierw test - sprawdź czy configs istnieją
ls -la configs/phase2_miller_urey_test.yaml

# Uruchom batch (jeśli plik istnieje)
python3 scripts/run_phase2_batch.py \
  --scenario miller_urey \
  --runs 3 \
  --output aws_test/results/test_batch \
  --config configs/phase2_miller_urey_test.yaml

# Sprawdź wyniki
ls aws_test/results/test_batch/
```

## Alternatywnie - Użyj configs w aws_test

Jeśli `aws_test/configs/` mają pliki:

```bash
cd ~/live2.0

# Sprawdź co jest
ls aws_test/configs/

# Uruchom z tym plikiem
python3 scripts/run_phase2_full.py \
  --config aws_test/configs/phase2_miller_urey_extended.yaml \
  --output aws_test/test_1 \
  --steps 10000 \
  --seed 42
```

