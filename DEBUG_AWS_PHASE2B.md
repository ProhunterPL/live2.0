# Debug AWS Phase 2B

## Problem
Symulacje się NIE uruchomiły, mimo że skrypt mówi "completed successfully".

## Sprawdź logi:

```bash
cd ~/live2.0/aws_test
cat results/phase2b_additional/logs/phase2b_runner.log
```

## Sprawdź co faktycznie jest w folderach run:

```bash
# Zobacz co jest w jednym z runów
ls -la results/phase2b_additional/miller_urey_extended/run_1/

# Zobacz logi formamide debug
ls -la results/phase2b_additional/formamide_debug/short_test/
```

## Najprawdopodobniej skrypty są template'ami

Musimy uruchomić **prawdziwą symulację** ręcznie:

```bash
cd ~/live2.0

# Test - czy podstawowy skrypt działa
python3 scripts/run_phase2_full.py \
  --config aws_test/configs/phase2_miller_urey_test.yaml \
  --output aws_test/test_manual \
  --steps 1000 \
  --seed 42
```

Jeśli to działa, to możemy:
1. Uruchomić batch przez `scripts/run_phase2_batch.py`
2. Lub użyć `scripts/phase2_master.py`

