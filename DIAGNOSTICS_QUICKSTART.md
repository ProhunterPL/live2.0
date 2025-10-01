# Diagnostics Quick Start

## Enable Diagnostics

Diagnostics are **enabled by default**. CSV files will be written to `./diagnostics/` directory.

## Running with Diagnostics

```python
from backend.sim.config import SimulationConfig
from backend.sim.core.stepper import SimulationStepper

# Create simulation (diagnostics enabled by default)
config = SimulationConfig()
sim = SimulationStepper(config)

# Run simulation
sim.start()
for _ in range(1000):
    sim.step()
sim.stop()  # Automatically closes diagnostics

print("Diagnostics saved to ./diagnostics/")
```

## Configuration Options

```python
config = SimulationConfig(
    enable_diagnostics=True,       # Enable/disable
    diagnostics_dir="my_run",      # Output directory
    diagnostics_frequency=10       # Log every N steps
)
```

## Analyzing Results

### Quick Analysis

```bash
python analyze_diagnostics.py diagnostics/
```

This generates:
- `bond_dynamics.png` - Bond formation/breaking dynamics
- `cluster_dynamics.png` - Cluster size evolution
- `phase_analysis.png` - Phase transitions and activity

### Python Analysis

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv('diagnostics/metrics_TIMESTAMP.csv')

# Plot bonds over time
plt.plot(df['sim_time'], df['num_bonds_total'])
plt.xlabel('Time')
plt.ylabel('Bonds')
plt.show()
```

## Output Files

- `metrics_*.csv` - Main time series (bonds, clusters, events)
- `bond_types_*.csv` - Bond type distribution
- `cluster_dist_*.csv` - Cluster size histogram
- `bond_lifetimes_*.csv` - Bond lifetime distribution
- `session_summary.txt` - Summary statistics

## What to Look For

### Healthy Simulation
- ✅ Oscillating bond count (dynamic equilibrium)
- ✅ Sustained event rate (formation/breaking)
- ✅ Cluster growth and fragmentation cycles
- ✅ Non-zero variance in energy

### Problem Signs
- ⚠️ Monotonic energy growth (not conserved)
- ⚠️ Zero events after startup (frozen system)
- ⚠️ Exponential bond growth (runaway)
- ⚠️ All particles isolated (no binding)

## Disable for Performance

```python
config = SimulationConfig(enable_diagnostics=False)
```

Saves ~5-10% CPU time.

## Full Documentation

See [docs/DIAGNOSTICS.md](docs/DIAGNOSTICS.md) for complete details.

