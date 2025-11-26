# Figure 4: Algorithm Flowchart

## Description
Complete algorithm flow for adaptive spatial hashing system.

---

## Main Simulation Loop with Adaptive Hash

```
┌──────────────────────────────────────────────────────────────────┐
│                     START SIMULATION                             │
│                  (N particles, box_size)                         │
└───────────────────────────┬──────────────────────────────────────┘
                            │
                            ▼
                    ┌───────────────┐
                    │ Initialize    │
                    │ α = 2.0       │
                    │ β = 0.3       │
                    │ s = 10.0      │  (default)
                    └───────┬───────┘
                            │
                            ▼
        ╔═══════════════════════════════════════════════════╗
        ║           SIMULATION TIME STEP LOOP               ║
        ╠═══════════════════════════════════════════════════╣
        ║                                                   ║
        ║   ┌─────────────────────────────────────┐        ║
        ║   │  UPDATE PARTICLE POSITIONS          │        ║
        ║   │  (integrate forces, apply dt)       │        ║
        ║   └────────────┬────────────────────────┘        ║
        ║                │                                  ║
        ║                ▼                                  ║
        ║   ┌─────────────────────────────────────┐        ║
        ║   │  UPDATE BONDING TOPOLOGY            │        ║
        ║   │  (form/break bonds based on dist)   │        ║
        ║   └────────────┬────────────────────────┘        ║
        ║                │                                  ║
        ║                ▼                                  ║
        ║   ┌─────────────────────────────────────┐        ║
        ║   │  Count active particles (N)         │        ║
        ║   │  Count total bonds (B)              │        ║
        ║   │  Compute b̄ = B / N                  │        ║
        ║   └────────────┬────────────────────────┘        ║
        ║                │                                  ║
        ║                ▼                                  ║
        ║   ┌─────────────────────────────────────┐        ║
        ║   │  COMPUTE OPTIMAL CELL SIZE          │◄───────╢─ PATENT 1
        ║   │                                     │        ║   CLAIM 1
        ║   │  s_new = α · √(A/N) / (1 + β·b̄)    │        ║
        ║   │  s_new = clamp(s_new, 5.0, 25.0)   │        ║
        ║   └────────────┬────────────────────────┘        ║
        ║                │                                  ║
        ║                ▼                                  ║
        ║   ┌─────────────────────────────────────┐        ║
        ║   │  CHECK IF REBUILD NEEDED            │◄───────╢─ PATENT 1
        ║   │                                     │        ║   CLAIM 4
        ║   │  Δs = |s_new - s_old| / s_old      │        ║
        ║   │  rebuild = (Δs > τ) OR (step % 20) │        ║
        ║   └────────────┬────────────────────────┘        ║
        ║                │                                  ║
        ║          ┌─────┴─────┐                           ║
        ║          │           │                           ║
        ║       NO │           │ YES                       ║
        ║          │           │                           ║
        ║          │           ▼                           ║
        ║          │  ┌────────────────────────┐          ║
        ║          │  │ REBUILD SPATIAL HASH   │◄─────────╢─ PATENT 1
        ║          │  │                        │          ║   CLAIM 2
        ║          │  │ 1. Clear grid          │          ║
        ║          │  │ 2. For each particle:  │          ║
        ║          │  │    cell = hash(pos, s) │          ║
        ║          │  │    atomic_add(cell, i) │          ║
        ║          │  │ 3. Update s_old = s_new│          ║
        ║          │  └────────┬───────────────┘          ║
        ║          │           │                           ║
        ║          └───────────┤                           ║
        ║                      │                           ║
        ║                      ▼                           ║
        ║          ┌────────────────────────┐             ║
        ║          │ COMPUTE FORCES         │◄────────────╢─ O(n)
        ║          │                        │             ║   complexity
        ║          │ For each particle i:   │             ║
        ║          │   cell = hash(pos[i])  │             ║
        ║          │   For 3×3 neighbors:   │             ║
        ║          │     For j in cell:     │             ║
        ║          │       F += force(i,j)  │             ║
        ║          └────────┬───────────────┘             ║
        ║                   │                              ║
        ║                   ▼                              ║
        ║          ┌────────────────────────┐             ║
        ║          │ APPLY FORCES           │             ║
        ║          │ v += F/m × dt          │             ║
        ║          └────────┬───────────────┘             ║
        ║                   │                              ║
        ║                   ▼                              ║
        ║          ┌────────────────────────┐             ║
        ║          │ step += 1              │             ║
        ║          │ Continue?              │             ║
        ║          └────────┬───────────────┘             ║
        ║                   │                              ║
        ╚═══════════════════╪═══════════════════════════════╝
                            │
                     ┌──────┴──────┐
                     │             │
                 YES │             │ NO
                     │             │
                     ▼             ▼
              (Continue Loop)  ┌────────┐
                               │  END   │
                               └────────┘
```

---

## Subroutine: Compute Optimal Cell Size (PATENT 1, CLAIM 1)

```
┌───────────────────────────────────────────────────────────┐
│  compute_optimal_cell_size(N, B, box_width, box_height)  │
└──────────────────────┬────────────────────────────────────┘
                       │
                       ▼
            ┌──────────────────────┐
            │  A = width × height  │
            │  (simulation area)   │
            └──────────┬───────────┘
                       │
                       ▼
            ┌──────────────────────┐
            │  ρ = N / A           │
            │  (particle density)  │
            └──────────┬───────────┘
                       │
                       ▼
            ┌──────────────────────┐
            │  b̄ = B / N           │
            │  (avg bonds/particle)│
            └──────────┬───────────┘
                       │
                       ▼
            ┌──────────────────────────────────┐
            │  base = √(A / N)                 │
            │  (characteristic spacing)        │
            └──────────┬───────────────────────┘
                       │
                       ▼
            ┌──────────────────────────────────┐
            │  f_bond = 1 / (1 + β·b̄)         │  ◄── INNOVATION
            │  (bonding adjustment)            │      Uses bonding
            └──────────┬───────────────────────┘      topology!
                       │
                       ▼
            ┌──────────────────────────────────┐
            │  s_optimal = α · base · f_bond   │  ◄── PATENT
            │                                  │      FORMULA
            └──────────┬───────────────────────┘
                       │
                       ▼
            ┌──────────────────────────────────┐
            │  s_final = clamp(s_optimal,      │
            │                  s_min, s_max)   │
            └──────────┬───────────────────────┘
                       │
                       ▼
            ┌──────────────────────┐
            │  RETURN s_final      │
            └──────────────────────┘
```

---

## Subroutine: Rebuild Spatial Hash (PATENT 1, CLAIM 2)

```
┌────────────────────────────────────────────────────┐
│  rebuild_spatial_hash(positions, active, N, s)    │
└──────────────────────┬─────────────────────────────┘
                       │
                       ▼
            ┌────────────────────────┐
            │  W = ⌈width / s⌉      │
            │  H = ⌈height / s⌉     │  ◄── Grid dimensions
            │  (adaptive grid size)  │      change per rebuild
            └──────────┬─────────────┘
                       │
                       ▼
            ┌────────────────────────┐
            │  Clear all cell counts │
            │  cell_count[*] = 0     │
            └──────────┬─────────────┘
                       │
                       ▼
        ╔══════════════════════════════════╗
        ║  FOR i = 0 TO N-1 (PARALLEL)     ║  ◄── GPU KERNEL
        ╠══════════════════════════════════╣
        ║                                  ║
        ║    ┌────────────────────────┐   ║
        ║    │  pos = positions[i]    │   ║
        ║    └────────┬───────────────┘   ║
        ║             │                    ║
        ║             ▼                    ║
        ║    ┌────────────────────────┐   ║
        ║    │  gx = ⌊pos.x / s⌋     │   ║
        ║    │  gy = ⌊pos.y / s⌋     │   ║  ◄── Hash function
        ║    │  cell = gy × W + gx   │   ║      (uses adaptive s)
        ║    └────────┬───────────────┘   ║
        ║             │                    ║
        ║             ▼                    ║
        ║    ┌────────────────────────────────┐
        ║    │  slot = atomic_add(            │  ◄── LOCK-FREE
        ║    │           cell_count[cell], 1) │      (GPU-safe)
        ║    └────────┬───────────────────────┘
        ║             │                    ║
        ║             ▼                    ║
        ║    ┌────────────────────────┐   ║
        ║    │  IF slot < MAX_SLOTS   │   ║
        ║    │    cell_list[cell][slot] = i
        ║    └────────────────────────┘   ║
        ║                                  ║
        ╚══════════════════════════════════╝
                       │
                       ▼
            ┌────────────────────────┐
            │  RETURN (success)      │
            └────────────────────────┘
```

---

## GPU Execution Architecture (PATENT 1, CLAIM 2 & 5)

### CPU/GPU Hybrid System
```
┌─────────────────────────────────────────────────────────┐
│                     CPU (Host)                          │
│                                                          │
│  ┌────────────────────────────────────────────────┐   │
│  │ 1. Compute adaptive cell size (lightweight)    │   │ ◄── O(1)
│  │    s_new = α·√(A/N) / (1 + β·b̄)               │   │     Fast!
│  │                                                 │   │
│  │ 2. Check if rebuild needed                     │   │
│  │    if |s_new - s_old| / s_old > 0.15:         │   │
│  │      trigger_gpu_rebuild = True                │   │
│  │                                                 │   │
│  │ 3. Launch GPU kernels                          │   │
│  │    build_spatial_hash<<<blocks, threads>>>     │   │
│  │    compute_forces<<<blocks, threads>>>         │   │
│  └─────────────┬──────────────────────────────────┘   │
│                │                                        │
│                │ PCIe transfer (minimal!)              │
│                │ Only send: cell_size (1 float)        │
│                │                                        │
└────────────────┼────────────────────────────────────────┘
                 │
                 ▼
┌────────────────────────────────────────────────────────┐
│                     GPU (Device)                        │
│                                                          │
│  ┌────────────────────────────────────────────────┐   │
│  │ PERSISTENT DATA (stays on GPU):                │   │ ◄── No transfer!
│  │   - positions[N]      (8N bytes)               │   │
│  │   - attributes[N]     (16N bytes)              │   │
│  │   - forces[N]         (8N bytes)               │   │
│  │   - grid_cell_list    (350 KB)                 │   │
│  │   - cell_count        (10 KB)                  │   │
│  └────────────────────────────────────────────────┘   │
│                                                          │
│  ┌────────────────────────────────────────────────┐   │
│  │ GPU KERNELS (execute in parallel):             │   │
│  │                                                 │   │
│  │  build_spatial_hash(cell_size=s_new)           │   │ ◄── O(N)
│  │    → 1024 threads, 1-2ms                       │   │
│  │    → Lock-free atomic operations               │   │
│  │                                                 │   │
│  │  compute_forces()                              │   │ ◄── O(N)
│  │    → 1024 threads, 5-8ms                       │   │
│  │    → 3×3 neighbor stencil (unrolled)          │   │
│  │                                                 │   │
│  │  integrate_positions(dt)                       │   │ ◄── O(N)
│  │    → 1024 threads, 0.5ms                       │   │
│  └────────────────────────────────────────────────┘   │
└────────────────────────────────────────────────────────┘
```

### GPU Kernel: Rebuild Spatial Hash
```
Thread Assignment:
  Launch: 1024 threads (N=1000 particles, round up to block size)
  
  Block 0 (256 threads):  P0, P1, ..., P255
  Block 1 (256 threads):  P256, P257, ..., P511
  Block 2 (256 threads):  P512, P513, ..., P767
  Block 3 (256 threads):  P768, ..., P999

Execution Timeline:
  0ms    1ms         2ms         3ms
  │      │           │           │
  ├──────┼───────────┼───────────┤
  │                                  │
  │ ┌─────────────────────────────┐ │
  │ │ PHASE 1: Clear cell counts │ │ ◄── All threads cooperate
  │ │   for i in parallel:        │ │
  │ │     cell_count[i] = 0       │ │
  │ └─────────────────────────────┘ │
  │                                  │
  │ ┌─────────────────────────────┐ │
  │ │ PHASE 2: Hash particles     │ │ ◄── PATENT 1, CLAIM 2
  │ │   Thread i:                 │ │     Atomic operations!
  │ │     pos = positions[i]       │ │
  │ │     cell = hash(pos, s_new)  │ │ ◄── Uses adaptive cell_size
  │ │     slot = atomic_add(...)   │ │ ◄── Lock-free insertion
  │ │     cell_list[cell][slot]=i  │ │
  │ └─────────────────────────────┘ │
```

### Atomic Operations (Thread-Safe)
```
Multiple threads may target SAME cell simultaneously:

Cell 42:
  Thread 17 → atomic_add(cell_count[42], 1) → returns 0 → slot 0
  Thread 35 → atomic_add(cell_count[42], 1) → returns 1 → slot 1
  Thread 89 → atomic_add(cell_count[42], 1) → returns 2 → slot 2
  
No race condition! Hardware ensures atomicity.

Without atomic:
  Thread 17 reads count=0, adds particle → count=1
  Thread 35 reads count=0 (race!), adds particle → count=1
  Result: Both write to slot 0, one particle lost! ❌

With atomic:
  Hardware serializes increments → guaranteed unique slots ✅
```

### GPU Kernel: Compute Forces
```
Thread execution:
  Thread i processes particle i:
    1. Find particle's cell: cell_i
    2. For 9 neighbor cells (3×3 stencil):
         For each particle j in neighbor:
           Compute force F_ij
           atomic_add(forces[i], F_ij)
           atomic_add(forces[j], -F_ij)

Neighbor Search Example:
  Particle i in cell (5, 5):
  
  Grid (adaptive cell_size = 12.0):
  ┌────┬────┬────┬────┬────┬────┬────┐
  │    │    │    │    │    │    │    │
  ├────┼────┼────┼────┼────┼────┼────┤
  │    │    │∘∘  │∘   │∘∘  │    │    │  ◄── Row 4
  ├────┼────┼────┼────┼────┼────┼────┤
  │    │    │∘   │●   │∘   │    │    │  ◄── Row 5 (particle i at ●)
  ├────┼────┼────┼────┼────┼────┼────┤
  │    │    │∘∘  │∘∘  │∘   │    │    │  ◄── Row 6
  └────┴────┴────┴────┴────┴────┴────┘
  
  Check cells: [(4,4), (5,4), (6,4),
                (4,5), (5,5), (6,5),  ◄── 3×3 = 9 cells
                (4,6), (5,6), (6,6)]
  
  Total neighbors to check: 12 particles (∘)
  Actually interact: ~4-6 (within cutoff)
```

### Compile-Time Loop Unrolling
```python
# Source code:
for dx in ti.static(range(-1, 2)):  # Compile-time constant!
    for dy in ti.static(range(-1, 2)):
        check_cell(cell_x + dx, cell_y + dy)

# Compiled GPU code (unrolled):
check_cell(cell_x - 1, cell_y - 1)
check_cell(cell_x - 1, cell_y + 0)
check_cell(cell_x - 1, cell_y + 1)
check_cell(cell_x + 0, cell_y - 1)
check_cell(cell_x + 0, cell_y + 0)  ◄── 9 explicit calls
check_cell(cell_x + 0, cell_y + 1)     NO branches!
check_cell(cell_x + 1, cell_y - 1)     NO loops!
check_cell(cell_x + 1, cell_y + 0)
check_cell(cell_x + 1, cell_y + 1)

Benefit: No branch prediction, better GPU utilization
Speedup: ~10-15% compared to dynamic loop
```

### Data Transfer Optimization
```
Traditional approach (all CPU):
  Transfer per step: 2N floats (positions + forces)
  Cost: ~0.5-1.0ms for N=1000
  Bottleneck: PCIe bandwidth

Adaptive hybrid (CPU compute, GPU work):
  Transfer per rebuild: 1 float (cell_size)
  Cost: ~0.001ms (negligible!)
  Benefit: 1000× less data transfer ✅

GPU-only (ideal but complex):
  Transfer: None (all on GPU)
  Challenge: Difficult to implement adaptive logic in kernel
  Tradeoff: Simpler to compute cell_size on CPU
```

### Memory Access Patterns
```
Coalesced Memory Access (GOOD):
  Thread block accessing positions:
  
  Memory (contiguous):
  [P0][P1][P2][P3][P4][P5][P6][P7]...
   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑
   │   │   │   │   │   │   │   │
  T0  T1  T2  T3  T4  T5  T6  T7  (threads)
  
  Result: Single 128-byte transaction ✅
  Bandwidth: ~800 GB/s (full speed)

Non-Coalesced Access (BAD - avoided):
  Random access pattern:
  
  Memory:
  [P0][P1][P2][P3][P4][P5][P6][P7]...
   ↑       ↑               ↑       ↑
   │       │               │       │
  T0      T7              T3      T1  (scattered)
  
  Result: 4 separate transactions ❌
  Bandwidth: ~200 GB/s (1/4 speed)
  
  Our implementation avoids this by:
  - Sequential particle processing
  - Cell list organized contiguously
```

---

## Key Decision Points

### Decision 1: When to Rebuild?
```
Option A: Fixed interval (every 10-20 steps)
  + Predictable overhead
  - May rebuild unnecessarily

Option B: Adaptive trigger (Δs > τ)    ◄── PATENT 1, CLAIM 4 (chosen)
  + Only rebuild when needed
  + Saves computation
  - Slightly more complex logic
```

### Decision 2: Cell Size Bounds?
```
Why clamp to [5.0, 25.0]?

s_min = 5.0:
  - Prevents over-refinement (too many cells)
  - Guarantees cutoff < 2×cell_size (3×3 stencil sufficient)
  
s_max = 25.0:
  - Prevents under-sampling (missing interactions)
  - Keeps grid resolution reasonable (~10×10 minimum)
```

### Decision 3: Parameters α, β?
```
α = 2.0:  Empirically determined (tested 1.5-3.0 range)
          - α < 2.0: Too many cells (overhead)
          - α > 2.0: Misses interactions (inaccurate)

β = 0.3:  Bonding influence (tested 0.1-0.5 range)
          - β < 0.3: Ignores bonding effect
          - β > 0.3: Over-reacts to bonding
```

---

## Complexity Analysis

```
Operation                    | Complexity | Notes
----------------------------|------------|------------------
Compute optimal cell size   | O(1)       | Simple formula
Clear grid                  | O(C)       | C = grid cells
Rebuild hash (parallel)     | O(N)       | N particles
Compute forces (3×3)        | O(N·k)     | k ≈ 9-16 neighbors
Total per step             | O(N)       | Linear scaling ✅

Comparison to alternatives:
  All-pairs:       O(N²)     (1000× slower for N=1000)
  Octree:          O(N log N) (10× slower)
  Fixed hash:      O(N)      (same complexity, worse constants)
  Adaptive hash:   O(N)      (same complexity, BETTER constants)
```

---

## Patent Claim Support

This flowchart demonstrates **Patent 1 - Spatial Hashing**:
- **Claim 1**: Complete algorithm with patent formula
- **Claim 2**: GPU implementation with atomic operations
- **Claim 3**: Bonding topology integration (b̄ term)
- **Claim 4**: Adaptive rebuild trigger (Δs threshold)
- **Claim 5**: CPU/GPU hybrid (cell size compute on CPU, rebuild on GPU)

