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
        ║   │  COMPUTE OPTIMAL CELL SIZE          │◄───────╢─ PATENT
        ║   │                                     │        ║   CLAIM 1
        ║   │  s_new = α · √(A/N) / (1 + β·b̄)    │        ║
        ║   │  s_new = clamp(s_new, 5.0, 25.0)   │        ║
        ║   └────────────┬────────────────────────┘        ║
        ║                │                                  ║
        ║                ▼                                  ║
        ║   ┌─────────────────────────────────────┐        ║
        ║   │  CHECK IF REBUILD NEEDED            │◄───────╢─ PATENT
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
        ║          │  │ REBUILD SPATIAL HASH   │◄─────────╢─ PATENT
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

## Subroutine: Compute Optimal Cell Size (PATENT CLAIM 1)

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

## Subroutine: Rebuild Spatial Hash (PATENT CLAIM 2)

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

## Key Decision Points

### Decision 1: When to Rebuild?
```
Option A: Fixed interval (every 10-20 steps)
  + Predictable overhead
  - May rebuild unnecessarily

Option B: Adaptive trigger (Δs > τ)    ◄── PATENT CLAIM 4 (chosen)
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

This flowchart demonstrates:
- **Claim 1**: Complete algorithm with patent formula
- **Claim 2**: GPU implementation with atomic operations
- **Claim 3**: Bonding topology integration (b̄ term)
- **Claim 4**: Adaptive rebuild trigger (Δs threshold)
- **Claim 5**: CPU/GPU hybrid (cell size compute on CPU, rebuild on GPU)

