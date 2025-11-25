# Diagnostics and Observables System

## Overview

The diagnostics system provides comprehensive time-series logging of simulation observables to CSV files. This allows detailed analysis of system dynamics, phase transitions, emergent behaviors, and energy conservation.

## Features

### Tracked Metrics

#### Bond Metrics
- **num_bonds_total**: Total number of active bonds
- **bonds_by_type**: Count of bonds by type (homogeneous, weakly heterogeneous, strongly heterogeneous)
- **avg_bond_length**: Average distance between bonded particles
- **avg_bond_tension**: Average deviation from equilibrium bond length
- **bond_lifetime_hist**: Distribution of bond lifetimes (logged periodically)

#### Cluster Metrics
- **num_clusters**: Total number of clusters
- **cluster_size_hist**: Distribution of cluster sizes (top-K clusters)
- **largest_cluster_size**: Size of the largest cluster
- **cluster_energy_mean**: Average energy per cluster
- **cluster_energy_var**: Variance in cluster energies
- **R_g_mean**: Average radius of gyration (spatial extent of clusters)

#### Event Tracking
- **events_formed**: Number of bonds formed this step
- **events_broken**: Number of bonds broken this step
- **events_merged**: Number of cluster merge events
- **events_split**: Number of cluster split events

## Configuration

### Enable/Disable Diagnostics

```python
from sim.config import SimulationConfig

config = SimulationConfig(
    enable_diagnostics=True,      # Enable logging
    diagnostics_dir="diagnostics", # Output directory
    diagnostics_frequency=10       # Log every N steps
)
```

### Configuration Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable_diagnostics` | bool | `True` | Enable/disable diagnostics logging |
| `diagnostics_dir` | str | `"diagnostics"` | Directory for CSV files |
| `diagnostics_frequency` | int | `10` | Log data every N steps |

## Output Files

The system generates timestamped CSV files in the diagnostics directory:

### 1. `metrics_TIMESTAMP.csv`

Main metrics file with columns:
- `step`: Simulation step number
- `sim_time`: Simulation time
- `wall_time`: Wall clock time (seconds since start)
- `num_particles`: Number of active particles
- `num_bonds_total`: Total bonds
- `avg_bond_length`: Average bond length
- `avg_bond_tension`: Average bond tension
- `num_clusters`: Number of clusters
- `largest_cluster_size`: Largest cluster size
- `cluster_energy_mean`: Mean cluster energy
- `cluster_energy_var`: Variance in cluster energy
- `R_g_mean`: Mean radius of gyration
- `events_formed`: Bonds formed
- `events_broken`: Bonds broken
- `events_merged`: Clusters merged
- `events_split`: Clusters split

### 2. `bond_types_TIMESTAMP.csv`

Bond type distribution:
- `step`: Simulation step
- `sim_time`: Simulation time
- `bond_type`: Type of bond (homogeneous, weakly_heterogeneous, strongly_heterogeneous)
- `count`: Number of bonds of this type

### 3. `cluster_dist_TIMESTAMP.csv`

Cluster size distribution (top-K clusters):
- `step`: Simulation step
- `sim_time`: Simulation time
- `cluster_size`: Size of cluster
- `count`: Number of clusters of this size

### 4. `bond_lifetimes_TIMESTAMP.csv`

Bond lifetime histogram (logged every 100 steps):
- `step`: Simulation step
- `sim_time`: Simulation time
- `lifetime`: Bond lifetime (binned)
- `count`: Number of bonds in this lifetime bin

### 5. `session_summary.txt`

Summary statistics for the entire session:
- Total steps
- Total simulation time
- Total bonds formed/broken
- Mean and median bond lifetime
- Active bonds at session end

## Usage Examples

### Basic Usage

```python
from sim.config import SimulationConfig
from sim.core.stepper import SimulationStepper

# Create configuration with diagnostics enabled
config = SimulationConfig(
    enable_diagnostics=True,
    diagnostics_dir="my_experiment",
    diagnostics_frequency=5  # Log every 5 steps
)

# Create and run simulation
sim = SimulationStepper(config)
sim.start()

# Run for 1000 steps
for _ in range(1000):
    sim.step()

# Stop (automatically closes diagnostics and writes summary)
sim.stop()
```

### Disable Diagnostics for Performance

```python
config = SimulationConfig(
    enable_diagnostics=False  # No logging overhead
)
```

### Custom Diagnostics Directory

```python
config = SimulationConfig(
    diagnostics_dir=f"experiments/run_{run_id}"
)
```

## Analysis Examples

### Python Analysis

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load main metrics
df = pd.read_csv('diagnostics/metrics_20250101_120000.csv')

# Plot bond count over time
plt.figure(figsize=(12, 6))
plt.plot(df['sim_time'], df['num_bonds_total'])
plt.xlabel('Simulation Time')
plt.ylabel('Number of Bonds')
plt.title('Bond Formation Dynamics')
plt.show()

# Plot cluster size evolution
plt.figure(figsize=(12, 6))
plt.plot(df['sim_time'], df['largest_cluster_size'], label='Largest')
plt.plot(df['sim_time'], df['num_clusters'], label='Total Clusters')
plt.xlabel('Simulation Time')
plt.ylabel('Count')
plt.legend()
plt.title('Cluster Dynamics')
plt.show()

# Analyze events
events = df[['events_formed', 'events_broken', 'events_merged', 'events_split']]
events.sum().plot(kind='bar')
plt.title('Total Events')
plt.show()
```

### Detecting Phases

```python
# Detect condensation phase
df['bond_rate'] = df['num_bonds_total'].diff() / df['sim_time'].diff()
condensation_phase = df[df['bond_rate'] > df['bond_rate'].quantile(0.9)]

# Detect oscillations
from scipy.signal import find_peaks
peaks, _ = find_peaks(df['num_bonds_total'], distance=20)
print(f"Oscillation period: ~{np.mean(np.diff(peaks))} steps")
```

## Performance Considerations

### Logging Frequency

- **Higher frequency (1-5 steps)**: More detailed data, higher overhead
- **Lower frequency (10-50 steps)**: Less overhead, sufficient for most analyses
- **Recommended**: 10 steps for standard runs, 1-5 for detailed studies

### Memory Usage

- CSV files are written incrementally and flushed every 10 steps
- Typical file sizes for 10,000 steps:
  - `metrics_*.csv`: ~500 KB
  - `bond_types_*.csv`: ~200 KB
  - `cluster_dist_*.csv`: ~300 KB
  - `bond_lifetimes_*.csv`: ~100 KB

### CPU Overhead

- Diagnostics logging adds ~5-10% overhead at frequency=10
- Overhead scales linearly with frequency
- Most overhead from numpy array conversions

## What to Look For

### System is "Alive" Indicators

1. **Oscillations** in bond count → dynamic equilibrium
2. **Phase transitions** in cluster size → condensation/fragmentation
3. **Sustained event rate** → ongoing reactions
4. **Stable energy with variance** → active but conserved

### Problem Indicators

1. **Monotonic energy growth** → energy not conserved
2. **Zero events after initial phase** → system frozen
3. **Exponential bond growth** → runaway aggregation
4. **All clusters size 1** → no binding occurring

## Advanced Features

### Custom Analysis Script

See `analyze_diagnostics.py` for a complete analysis pipeline including:
- Time series plots
- Phase detection
- Event rate analysis
- Energy conservation checks
- Cluster size distributions

### Integration with Jupyter

```python
import pandas as pd
from IPython.display import display

# Interactive analysis
df = pd.read_csv('diagnostics/metrics_*.csv')
display(df.describe())
df.plot(x='sim_time', y=['num_bonds_total', 'num_clusters'], 
        subplots=True, figsize=(12, 8))
```

## Troubleshooting

### Issue: No CSV files generated

**Solution**: Check `enable_diagnostics=True` in config

### Issue: Files are empty

**Solution**: Run simulation for at least one logging interval

### Issue: High memory usage

**Solution**: Reduce `diagnostics_frequency` or disable logging for production runs

### Issue: ValueError in diagnostics

**Solution**: Check that simulation has particles and bonds before logging starts

## Future Extensions

Planned features:
- Real-time plotting dashboard
- HDF5 format for large datasets
- Configurable metric subsets
- Automatic anomaly detection
- Integration with ML pipelines

