# ⚡ Szybki Start - Phase 2B Local

## 🎯 TL;DR

```powershell
# 1. Test (5 min)
python scripts/run_phase2_full.py --config aws_test/configs/phase2_miller_urey_extended.yaml --output results/test_local_miller_urey --steps 10000 --seed 42

# 2. Pełna symulacja (1-2h)
python scripts/run_phase2_full.py --config aws_test/configs/phase2_miller_urey_extended.yaml --output results/phase2b_local/miller_urey/run_01 --steps 500000 --seed 100

# 3. Batch 10 symulacji (10-20h)
python run_phase2b_local.py --scenario miller_urey --runs 10

# 4. Wszystko (30 symulacji, 2-3 dni)
python run_phase2b_local.py --all --runs 10
```

## 📊 Commands

| Komenda | Czas | Opis |
|---------|------|------|
| `test_gpu_quick.py` | 1 min | Test GPU |
| 10K steps | 5 min | Test symulacji |
| 500K steps | 1-2h | 1 pełna symulacja |
| 10 runs | 10-20h | Batch scenariusz |
| 30 runs | 2-3 dni | Kompletne Phase 2B |

## ✅ Status

- ✅ GPU działa (RTX 5070 zweryfikowane)
- ✅ Skrypty gotowe
- ✅ Konfiguracje przygotowane
- ✅ Gotowe do uruchomienia

Zacznij od testu (5 min)!

