# Uruchomienie na RTX 5070 - Lokalnie

## Przygotowanie

### 1. Aktywuj Środowisko
```powershell
cd D:\live2.0
.\activate_live_env.ps1  # Jeśli masz virtualenv
# lub
pip install numpy taichi pyyaml
```

### 2. Test GPU
```powershell
python -c "import taichi as ti; ti.init(arch=ti.cuda); print('GPU OK!')"
```

### 3. Uruchom Phase 2 Batch

```powershell
cd D:\live2.0

# Test - 1 symulacja (10K kroków)
python scripts/run_phase2_full.py `
  --config configs/phase2_miller_urey_test.yaml `
  --output results/test_local `
  --steps 10000 `
  --seed 42
```

### 4. Jeśli to działa - Uruchom Pełny Batch

```powershell
# 10 symulacji Miller-Urey
python scripts/run_phase2_batch.py `
  --scenario miller_urey `
  --runs 10 `
  --output results/phase2_local/miller_urey `
  --config configs/phase2_miller_urey.yaml

# 10 symulacji Hydrothermal  
python scripts/run_phase2_batch.py `
  --scenario hydrothermal `
  --runs 10 `
  --output results/phase2_local/hydrothermal `
  --config configs/phase2_hydrothermal.yaml

# 10 symulacji Formamide
python scripts/run_phase2_batch.py `
  --scenario formamide `
  --runs 10 `
  --output results/phase2_local/formamide `
  --config configs/phase2_formamide.yaml
```

### 5. Edytuj run_phase2_full.py żeby używał GPU

Edytuj `scripts/run_phase2_full.py` linia 203:

```python
def initialize_taichi(self):
    logger.info("Initializing Taichi...")
    
    try:
        ti.init(arch=ti.cuda)
        logger.info("✅ Using GPU acceleration")
    except Exception as e:
        logger.warning(f"⚠️ GPU not available: {e}")
        import multiprocessing
        max_threads = multiprocessing.cpu_count()
        ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)
        logger.info(f"Using CPU with {max_threads} threads")
```

## Szacowany Czas na RTX 5070

- 1 symulacja (500K kroków): ~1-2 godziny (vs 20+ godzin na CPU)
- 30 symulacji: ~30-60 godzin (2-3 dni)
- Możesz uruchomić równolegle 2-3 symulacje

## Równoległe Uruchomienie

```powershell
# Uruchom 3 scenariusze równolegle (PowerShell)
Start-Job -ScriptBlock { python scripts/run_phase2_batch.py --scenario miller_urey --runs 10 --output results/phase2_local/miller_urey --config configs/phase2_miller_urey.yaml }
Start-Job -ScriptBlock { python scripts/run_phase2_batch.py --scenario hydrothermal --runs 10 --output results/phase2_local/hydrothermal --config configs/phase2_hydrothermal.yaml }  
Start-Job -ScriptBlock { python scripts/run_phase2_batch.py --scenario formamide --runs 10 --output results/phase2_local/formamide --config configs/phase2_formamide.yaml }

# Sprawdź status
Get-Job
```

