# Optymalizacja Symulacji

## Problem
Symulacja działa, ale detekcja klastrów (`detect_novel_substances`) trwa ~50 sekund co 500 kroków.

## Rozwiązanie

### 1. Wyłącz lub Zmniejsz Częstotliwość Detekcji

Edytuj backend/sim/core/stepper.py i znajdź:
```python
# Zwiększ interwał detekcji z 500 na 10000 kroków
if step % 500 == 0:  # zmień na:
if step % 10000 == 0:
```

### 2. Włącz GPU (jeśli instancja ma GPU)

Edytuj scripts/run_phase2_full.py linia 203:
```python
# Zmień z:
ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)

# Na:
try:
    ti.init(arch=ti.cuda)
    logger.info("Using GPU acceleration")
except:
    ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)
    logger.info("Using CPU (GPU not available)")
```

### 3. Na AWS - Sprawdź czy masz GPU

```bash
nvidia-smi
```

### 4. Jeśli NIE masz GPU - Użyj większej instancji

```bash
# Na AWS, zatrzymaj obecną instancję
# Uruchom nową: g4dn.xlarge (GPU NVIDIA T4)
# lub: p3.2xlarge (GPU NVIDIA V100)
```

### 5. Albo - Uruchom Bez Detekcji (najszybsze)

Stwórz minimalny runner który NIE włącza detekcji klastrów:

```python
# quick_run.py
from backend.sim.core.stepper import SimulationStepper
from backend.sim.run_simulation import create_simulation

# Uruchom bez detekcji
sim = create_simulation(config)
stepper = SimulationStepper(sim_config)

# WYŁĄCZ detekcję dla szybkości
stepper.detection_interval = 100000000  # Prawie nigdy

for i in range(steps):
    stepper.step()
```

