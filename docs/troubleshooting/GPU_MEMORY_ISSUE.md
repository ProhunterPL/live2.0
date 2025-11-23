# GPU Memory Issue - RTX 5070 + Taichi

## Problem

RTX 5070 kończy pamięć GPU podczas benchmarków LIVE 2.0:

```
Taichi JIT:0: allocate_from_reserved_memory: block: [0,0,0], thread: [0,0,0] 
Assertion `Out of CUDA pre-allocated memory.
Consider using ti.init(device_memory_fraction=0.9) or ti.init(device_memory_GB=4) 
to allocate more GPU memory` failed.
```

## Analiza

### GPU ma 12GB VRAM
```
Memory-Usage | GPU-Util
983MiB / 12227MiB | 0%
```

**Problem nie jest w braku pamięci!** (11GB wolne)

### Prawdziwy problem

1. **Taichi pre-allocation bug**
   - Taichi rezerwuje pamięć z góry
   - Podczas długiego benchmarku alokacja się kończy
   - Bug w zarządzaniu pamięcią Taichi

2. **GPU jest WOLNE na chemii**
   ```
   Bonds/Clusters: 11659.1ms (11.6 sekund!) ❌
   ```
   GPU nie jest dobre w:
   - Złożonej logice (if/else)
   - Grafach (pointer chasing)
   - Nieregularnym dostępie do pamięci

## Potwierdzenie z Phase 2B

Z testów AWS:
- **CPU (96 threads)**: Szybszy dla chemii
- **GPU**: Wolniejszy dla operacji grafowych

**Wniosek:** CPU jest lepszy dla tego workloadu!

## Rozwiązanie

### Opcja 1: Użyj Pure CPU (REKOMENDOWANE)

```python
import taichi as ti
import multiprocessing

ti.init(arch=ti.cpu, cpu_max_num_threads=multiprocessing.cpu_count())

from backend.sim.core.stepper import SimulationStepper
stepper = SimulationStepper(config)
```

**Dlaczego:**
- ✅ Brak problemów z pamięcią
- ✅ Szybszy dla chemii (Phase 2B proof)
- ✅ Stabilny
- ✅ Masz wiele rdzeni CPU

### Opcja 2: GPU tylko dla wizualizacji

Użyj GPU tylko do renderowania, CPU do symulacji:

```python
# Symulacja na CPU (stabilna)
ti.init(arch=ti.cpu, cpu_max_num_threads=multiprocessing.cpu_count())

# Frontend używa WebGL dla renderowania
# (GPU używane tylko przez przeglądarkę)
```

### Opcja 3: Hybrid (gdy GPU zadziała)

Jeśli naprawisz GPU:
```python
from backend.sim.core.hybrid_stepper import HybridSimulationStepper

# GPU: fizyka
# CPU: chemia
stepper = HybridSimulationStepper(config)
```

Ale wymaga działającego GPU!

## Próby naprawy GPU

### ❌ Nie pomogło:

1. **Downgrade do Taichi 1.6.0**
   - Nadal ten sam problem
   - To nie jest bug wersji

2. **Więcej pamięci (`device_memory_GB=4.0`)**
   - Nadal crashuje
   - Problem nie w ilości pamięci

3. **Mniej cząstek (30-50)**
   - Pomaga trochę
   - Ale nadal crashuje w długich testach

### 💡 Co może pomóc:

1. **Taichi z CUDA 11 zamiast CUDA 12**
   ```bash
   pip install taichi-nightly
   ```

2. **Starszy driver NVIDIA**
   - RTX 5070 jest nowe (2025)
   - Może być niekompatybilne

3. **Różne dystrybucje Taichi**
   ```bash
   # Próbuj różnych wersji
   pip install taichi==1.5.0
   pip install taichi==1.4.0
   ```

## Benchmark bez GPU

Uruchom benchmark tylko CPU:

```powershell
.\run_cpu_benchmark.ps1
```

Lub:
```powershell
python tests\benchmark_hybrid.py --modes pure_cpu
```

## Statystyki z testów

### Pure GPU (przed crashem)
- Fizyka: OK (~2-5ms/step)
- **Bonds/Clusters: 11659ms** ❌ BARDZO WOLNE
- Crash po ~200 steps

### Pure CPU (przewidywane)
- Fizyka: ~10-20ms/step (wolniej niż GPU)
- **Bonds/Clusters: ~50-100ms** ✅ ZNACZNIE SZYBCIEJ
- Stabilne, bez crashy

### Hybrid (teoretycznie)
- Fizyka (GPU): ~5ms
- Chemia (CPU): ~100ms (async)
- **Total effective: ~5ms** ✅
- Ale wymaga działającego GPU!

## Rekomendacja

**Dla Twojego systemu: Użyj Pure CPU**

```powershell
# Benchmark CPU
.\run_cpu_benchmark.ps1

# W kodzie produkcyjnym
ti.init(arch=ti.cpu, cpu_max_num_threads=multiprocessing.cpu_count())
```

**Dlaczego:**
1. ✅ Stabilne (brak crashy)
2. ✅ Szybsze dla chemii (z Phase 2B)
3. ✅ RTX 5070 ma problemy z Taichi
4. ✅ Masz wiele rdzeni CPU

**GPU możesz używać tylko do:**
- Wizualizacji w przeglądarce (WebGL)
- Innych projektów (nie Taichi)

## Przyszłość

Gdy RTX 5070 będzie lepiej wspierane przez Taichi (6-12 miesięcy):
- Spróbuj ponownie GPU mode
- Lub użyj Hybrid mode
- Monitoruj Taichi release notes

## Linki

- [Taichi GPU Memory Issues](https://github.com/taichi-dev/taichi/issues)
- [RTX 50xx Support Status](https://github.com/taichi-dev/taichi/discussions)
- [CPU Performance Guide](PERFORMANCE_TUNING.md)

---

**TL;DR:** RTX 5070 + Taichi = problemy. Użyj CPU, jest stabilniejsze i szybsze dla chemii.

