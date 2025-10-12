# LIVE 2.0 - Analiza Wydajności i Optymalizacje

## 🎯 Cel: 2 FPS → 60 FPS (30x przyspieszenie)

## 📊 Zidentyfikowane Bottlenecki

### 1. **MetricsCollector - Ręczna pętla Python** ⚠️ KRYTYCZNY
**Lokalizacja:** `backend/sim/core/metrics.py:178-195`

**Problem:**
```python
for i in range(min(particle_count, 500)):  # Python loop!
    if active[i] == 1:
        # Obliczenia w Python zamiast Taichi kernel
```

**Koszt:** ~50-100ms na wywołanie (w Python!)

**Rozwiązanie:** Przenieść do Taichi kernel
- ✅ GPU parallelization
- ✅ 100-1000x szybciej

### 2. **Energy Field Transfer** ⚠️ WYSOKIE
**Lokalizacja:** `backend/sim/core/stepper.py:1201`

**Problem:**
```python
data['energy_field'] = self.energy_manager.energy_system.energy_field.to_numpy()
```

**Koszt:** ~20-50ms (kopiowanie GPU→CPU przy każdym framie)

**Rozwiązanie:**
- Cache energy_field (aktualizuj co N kroków)
- Zmniejsz rozdzielczość dla przesyłania (downsample)
- Kompresja przed wysłaniem

### 3. **Particles/Bonds Get Operations** ⚠️ ŚREDNIE
**Lokalizacja:** `backend/sim/core/stepper.py:1226, 1251`

**Problem:**
- `get_active_particles()` - kopiowanie GPU→CPU
- `get_bonds()` - kopiowanie dużych struktur

**Koszt:** ~10-30ms łącznie

**Rozwiązanie:**
- Już jest cache (co 5/20 kroków) ✅
- Dalsze zwiększenie interwału cache

### 4. **WebSocket msgpack Encoding** ⚠️ ŚREDNIE
**Lokalizacja:** `backend/api/server.py` (broadcast)

**Problem:**
- msgpack serialization dużych numpy arrays
- Wysyłanie pełnych danych co frame

**Koszt:** ~20-40ms

**Rozwiązanie:**
- Delta encoding (wysyłaj tylko zmiany)
- Kompresja (zlib/lz4)
- Downsampling danych

### 5. **Throttled Operations** ℹ️ NISKI
**Obecny stan:**
- Bonds: co 100 kroków ✅
- Clusters: co 200 kroków ✅
- Mutations: co 200 kroków ✅
- Metrics: co 200 kroków ✅
- Novelty: co 500 kroków ✅

**Możliwa optymalizacja:** Dalsze zwiększenie interwałów

## 🚀 Plan Optymalizacji (Priorytet)

### Faza 1: Quick Wins (2x-5x przyspieszenie)
1. **MetricsCollector → Taichi kernel** (50-100ms → 1-2ms)
2. **Cache energy_field** (50ms → 5ms średnio)
3. **Zwiększ cache intervals** (particles: 5→10, bonds: 20→50)

### Faza 2: Medium (2x-3x dodatkowe)
4. **Downsample energy_field** dla wizualizacji (256x256 → 128x128)
5. **WebSocket compression** (zlib level 1)
6. **Optymalizuj particles.get_active_particles()**

### Faza 3: Advanced (1.5x-2x dodatkowe)
7. **Delta encoding** dla WebSocket
8. **Adaptive quality** (dynamiczna rozdzielczość)
9. **GPU-accelerated compression**

## 📈 Szacowane Wyniki

| Faza | Obecny FPS | Cel FPS | Przyspieszenie |
|------|-----------|---------|----------------|
| Start | 2 | - | 1x |
| Faza 1 | 2 | 8-15 | 4-7x |
| Faza 2 | 8-15 | 20-35 | 2-3x |
| Faza 3 | 20-35 | 30-50 | 1.5-2x |
| **TOTAL** | **2** | **30-50** | **15-25x** |

## 🔍 Pomiary Wydajności (Obecne)

### Breakdown czasów (estymacja):
```
Total Frame Time: ~500ms (2 FPS)
├─ Simulation Step: ~200ms (40%)
│  ├─ Physics (compute_forces): 80ms
│  ├─ Binding update: 60ms
│  ├─ Energy operations: 40ms
│  └─ Other: 20ms
│
├─ Get Visualization Data: ~200ms (40%)
│  ├─ MetricsCollector (Python loop): 80ms ⚠️
│  ├─ Energy field to_numpy(): 50ms ⚠️
│  ├─ Particles extraction: 30ms
│  ├─ Bonds/Clusters extraction: 30ms
│  └─ Other: 10ms
│
└─ WebSocket Broadcast: ~100ms (20%)
   ├─ msgpack encoding: 60ms ⚠️
   └─ Network send: 40ms
```

## 🛠️ Implementacja

### 1. MetricsCollector Taichi Kernel

```python
@ti.kernel
def update_particle_metrics_kernel_optimized(
    active: ti.template(),
    attributes: ti.template(), 
    energy: ti.template(),
    velocities: ti.template(),
    particle_count: ti.i32,
    result: ti.template()  # [total_energy, total_mass, active_count]
):
    # Reset accumulators
    result[0] = 0.0  # total_energy
    result[1] = 0.0  # total_mass
    result[2] = 0    # active_count
    
    # Parallel reduction
    for i in range(particle_count):
        if active[i] == 1:
            ti.atomic_add(result[2], 1)
            mass = attributes[i][0]
            ti.atomic_add(result[1], mass)
            
            # Energy
            particle_energy = energy[i]
            vx, vy = velocities[i][0], velocities[i][1]
            kinetic = 0.5 * mass * (vx * vx + vy * vy)
            ti.atomic_add(result[0], particle_energy + kinetic)
```

### 2. Energy Field Cache

```python
class EnergyFieldCache:
    def __init__(self, update_interval=10):
        self.update_interval = update_interval
        self.cached_field = None
        self.last_update_step = -999
        
    def get(self, step_count, energy_system):
        if step_count - self.last_update_step >= self.update_interval:
            self.cached_field = energy_system.energy_field.to_numpy()
            self.last_update_step = step_count
        return self.cached_field
```

### 3. Zwiększone Throttling

```python
# W _perform_step():
if self.step_count % 200 == 0:  # było 100
    self.binding.update_bonds(...)
    
if self.step_count % 400 == 0:  # było 200
    self.binding.update_clusters(...)

# W get_visualization_data():
if self.step_count % 10 == 0:  # było 5
    positions, velocities, attributes, active_mask, energies = self.particles.get_active_particles()
    
if self.step_count % 50 == 0:  # było 20
    bonds = self.binding.get_bonds()
    clusters = self.binding.get_clusters()
```

## ⚡ Instant Actions (bez zmiany kodu)

1. **Zmniejsz rozdzielczość gridu**: 256x256 → 128x128
2. **Zmniejsz max_particles**: 10000 → 5000
3. **Zwiększ dt**: 0.01 → 0.02 (większe kroki czasowe)

## 📝 Monitoring

Dodaj do logów:
```python
logger.info(f"Performance breakdown:")
logger.info(f"  Step: {step_time*1000:.1f}ms")
logger.info(f"  Viz: {viz_time*1000:.1f}ms")
logger.info(f"  Broadcast: {broadcast_time*1000:.1f}ms")
logger.info(f"  FPS: {1.0/(step_time+viz_time+broadcast_time):.1f}")
```

## 🎯 Następne Kroki

1. Zacznij od Fazy 1 (największy zysk)
2. Mierz FPS po każdej zmianie
3. Iteruj do osiągnięcia 60 FPS
4. Zachowaj wszystkie funkcjonalności

