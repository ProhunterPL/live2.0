# Naprawa Ścieżki - Szybkie Rozwiązanie

## Problem
Skrypty Phase 2B szukają: `/home/ubuntu/live2.0/aws_test/scripts/run_phase2_full.py`
Ale plik jest w: `/home/ubuntu/live2.0/scripts/run_phase2_full.py`

## Rozwiązanie - Na AWS:

```bash
cd ~/live2.0/aws_test

# Stwórz symlink do właściwego skryptu
ln -s ../scripts/run_phase2_full.py scripts/run_phase2_full.py

# Sprawdź
ls -la scripts/run_phase2_full.py

# Teraz uruchom ponownie
python3 run_phase2b_master.py --mode all
```

## Alternatywnie - Użyj działający skrypt batch

Jeśli to nie zadziała, użyj bezpośrednio:

```bash
cd ~/live2.0

# Uruchom 5 symulacji Miller-Urey (test)
python3 scripts/run_phase2_batch.py \
  --scenario miller_urey \
  --runs 5 \
  --output aws_test/results/phase2_real \
  --config aws_test/configs/phase2_miller_urey.yaml

# Jeśli to zadziała, uruchom wszystkie scenariusze
python3 scripts/run_phase2_batch.py --all --runs 10 \
  --output aws_test/results/phase2_real
```

