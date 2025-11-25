# Performance Tuning Guide 🚀

This guide helps you optimize LIVE 2.0 performance for your specific hardware.

## Quick Start

### 1. Check Your Current Backend

```bash
python scripts/check_current_backend.py
```

This shows which Taichi backend (CUDA/Vulkan/CPU) your system will use.

### 2. Run Benchmark

**Windows:**
```powershell
.\run_benchmark.ps1
```

**Linux/Mac:**
```bash
python tests/benchmark_gpu_vs_cpu.py
```

The benchmark tests:
- ✅ Simulation steps per second
- ✅ Visualization frame time
- ✅ Memory usage

Takes ~5-10 minutes to complete.

## Understanding Results

### Simulation Performance

**Steps/sec** = How many simulation timesteps can be computed per second

- **100-500 steps/sec**: Excellent (GPU)
- **50-100 steps/sec**: Good (fast CPU or older GPU)
- **10-50 steps/sec**: Acceptable (CPU with many cores)
- **<10 steps/sec**: Slow (single-threaded CPU)

### Visualization Performance

**ms/frame** = How long it takes to render one frame

- **<50ms**: Excellent (20+ FPS)
- **50-100ms**: Good (10-20 FPS)
- **100-200ms**: Acceptable (5-10 FPS)
- **>200ms**: Slow (stuttering)

## Hardware-Specific Recommendations

### NVIDIA GPU (CUDA)

**Best for:**
- ✅ Large particle counts (>1000)
- ✅ Real-time visualization
- ✅ Interactive simulations

**Configuration:**
```python
# backend/api/server.py uses this automatically
ti.init(arch=ti.cuda, device_memory_GB=4.0)
```

**Expected performance:**
- 100-500 steps/sec depending on GPU
- 10-50ms visualization time

### AMD/Intel GPU (Vulkan)

**Best for:**
- ✅ Systems without NVIDIA GPU
- ✅ Moderate particle counts

**Expected performance:**
- 50-200 steps/sec
- 50-100ms visualization time

### CPU (Many Cores)

**Best for:**
- ✅ Chemistry-heavy workloads
- ✅ Batch processing
- ✅ Cloud instances with many vCPUs

**Configuration:**
```python
import multiprocessing
ti.init(arch=ti.cpu, cpu_max_num_threads=multiprocessing.cpu_count())
```

**Expected performance:**
- 10-100 steps/sec (depends on core count)
- 100-500ms visualization time

**Note:** In our Phase 2B tests on AWS (96 vCPUs), CPU was **faster** than GPU because:
- Chemistry calculations don't parallelize well on GPU
- Many CPU cores can handle parallel independent tasks
- Lower GPU-CPU memory transfer overhead

## Optimization Tips

### For GPU Users

1. **Reduce memory transfers:**
   - Cache visualization data (already implemented)
   - Reduce frequency of `to_numpy()` calls

2. **Tune particle count:**
   - GPU excels at >1000 particles
   - Below that, CPU may be competitive

3. **Monitor GPU memory:**
   ```python
   ti.init(arch=ti.cuda, device_memory_GB=4.0)  # Adjust as needed
   ```

### For CPU Users

1. **Use all available cores:**
   ```python
   ti.init(arch=ti.cpu, cpu_max_num_threads=multiprocessing.cpu_count())
   ```

2. **Reduce visualization frequency:**
   - Set `energy_field_cache_interval` higher
   - Update particles every 10-20 steps instead of every step

3. **Disable heavy operations during batch runs:**
   ```python
   config.detect_novel_substances = False  # Until needed
   config.metrics_update_interval = 1000   # Less frequent
   config.novelty_check_interval = 5000    # Much less frequent
   ```

## Real-World Scenarios

### Interactive Local Development

**Best setup:**
- GPU (CUDA or Vulkan)
- 500-1000 particles
- Real-time visualization

**Expected:**
- Smooth 10-30 FPS
- <100ms per frame

### Batch Production Runs

**Best setup:**
- CPU with many cores (if available)
- 500-2000 particles
- Visualization disabled or minimal

**Expected:**
- 20-100 steps/sec
- Can run 24/7 for long experiments

### Cloud Computing

**Best setup depends on instance:**
- **GPU instance (g4dn, p3)**: Use CUDA
- **CPU-heavy instance (c6i, c7i)**: Use CPU with all cores
- **Balanced instance (m5, m6i)**: Test both!

## Benchmark Results Interpretation

### Example Output

```
📊 Simulation Performance (steps/sec):
1. CUDA GPU          - 312.5 steps/sec (3.2ms/step) [8.9x faster]
2. CPU (96 threads)  - 87.3 steps/sec (11.5ms/step) [2.5x faster]
3. CPU (48 threads)  - 52.1 steps/sec (19.2ms/step) [1.5x faster]
4. CPU (24 threads)  - 35.2 steps/sec (28.4ms/step) [1.0x faster]

📊 Visualization Performance (lower is better):
1. CUDA GPU          - 42.3ms avg [3.2x faster]
2. CPU (96 threads)  - 98.7ms avg [1.4x faster]
3. CPU (48 threads)  - 124.5ms avg [1.1x faster]
4. CPU (24 threads)  - 135.6ms avg [1.0x faster]
```

**Interpretation:**
- GPU is fastest for both simulation and visualization
- More CPU threads = better performance (until memory bandwidth limit)
- For this system: **Use CUDA GPU**

### Example Output (CPU Wins)

```
📊 Simulation Performance (steps/sec):
1. CPU (96 threads)  - 124.3 steps/sec (8.0ms/step) [2.1x faster]
2. CUDA GPU          - 89.7 steps/sec (11.1ms/step) [1.5x faster]
3. CPU (48 threads)  - 67.2 steps/sec (14.9ms/step) [1.1x faster]
4. CPU (24 threads)  - 59.8 steps/sec (16.7ms/step) [1.0x faster]
```

**Interpretation:**
- CPU with many threads beats GPU for simulation steps
- Likely due to chemistry-heavy workload (not compute-bound)
- For batch runs: **Use CPU (96 threads)**
- For interactive use: Test visualization performance too

## Troubleshooting

### GPU Not Detected

```bash
python scripts/check_current_backend.py
```

If CUDA fails:
1. Check NVIDIA driver: `nvidia-smi`
2. Install CUDA toolkit
3. Reinstall taichi: `pip install --upgrade taichi`

### Poor GPU Performance

Possible causes:
- GPU is busy with other tasks (close apps)
- Old GPU driver (update)
- Small particle count (GPU overhead dominates)
- Memory transfer bottleneck (reduce visualization frequency)

### Poor CPU Performance

Possible causes:
- Not using all threads (check initialization)
- Memory bottleneck (reduce particle count)
- Other processes competing (close apps)
- Thermal throttling (check CPU temperature)

## Advanced Configuration

### Custom Backend Selection

Edit `backend/api/server.py`:

```python
# Force specific backend
ti.init(arch=ti.cuda)  # Always GPU
# or
ti.init(arch=ti.cpu, cpu_max_num_threads=96)  # Always CPU
```

### Per-Scenario Optimization

Different scenarios may benefit from different backends:

**Visualization-heavy (real-time UI):**
→ Use GPU

**Chemistry-heavy (complex reactions):**
→ Test both, CPU may win with many cores

**Large particle count (>2000):**
→ GPU usually wins

**Small particle count (<500):**
→ CPU may be competitive

## Conclusion

**The best backend depends on:**
1. Your hardware (GPU model, CPU core count)
2. Your workload (interactive vs batch, particle count)
3. Your priorities (speed vs cost vs energy efficiency)

**Always benchmark your specific setup!**

Run: `.\run_benchmark.ps1` (Windows) or `python tests/benchmark_gpu_vs_cpu.py` (Linux/Mac)

---

For questions or issues, see:
- [README.md](../README.md) - Main documentation
- [QUICK_START.md](QUICK_START.md) - Getting started
- [GitHub Issues](https://github.com/yourusername/live2.0/issues) - Report problems

