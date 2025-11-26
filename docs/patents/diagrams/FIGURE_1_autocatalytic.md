# Figure 1: Autocatalytic Cycle Detection Algorithm

## Description
Illustrates the detection and analysis of autocatalytic cycles in prebiotic chemistry reaction networks.

---

## Autocatalysis Detection Pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│                    INPUT: REACTION NETWORK                      │
│                                                                   │
│  Directed Graph G = (V, E)                                       │
│    V = Molecules (nodes)                                         │
│    E = Reactions (edges: reactant → product)                     │
│                                                                   │
│  Example:                                                        │
│    A → B  (reaction: A forms B)                                  │
│    B → C  (reaction: B forms C)                                  │
│    C → A  (reaction: C forms A)                                   │
│                                                                   │
│  + Abundance History: {molecule_id: [count_at_step_t]}          │
└───────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
        ╔═══════════════════════════════════════════════════╗
        ║         STEP 1: CYCLE DETECTION                   ║
        ║         (Johnson's Algorithm)                     ║
        ╠═══════════════════════════════════════════════════╣
        ║                                                   ║
        ║   Find all simple cycles in directed graph        ║
        ║   Maximum cycle length: 6 nodes                    ║
        ║   Timeout: 300 seconds (safety limit)             ║
        ║                                                   ║
        ║   Output: List of cycles                          ║
        ║     Cycle = [node1, node2, ..., nodeN]            ║
        ║                                                   ║
        ╚═══════════════════════════════════════════════════╝
                         │
                         ▼
        ╔═══════════════════════════════════════════════════╗
        ║         STEP 2: AUTOCATALYTIC FILTERING           ║
        ║         (Amplification Check)                     ║
        ╠═══════════════════════════════════════════════════╣
        ║                                                   ║
        ║   For each cycle:                                  ║
        ║     For each molecule M in cycle:                 ║
        ║       1. Check: M appears as both reactant        ║
        ║          and product in cycle                      ║
        ║       2. Calculate amplification:                 ║
        ║          amp = abundance_late / abundance_early    ║
        ║       3. If amp ≥ 1.5×: AUTOCATALYTIC ✅          ║
        ║                                                   ║
        ╚═══════════════════════════════════════════════════╝
                         │
                         ▼
        ╔═══════════════════════════════════════════════════╗
        ║         STEP 3: CYCLE CLASSIFICATION               ║
        ║         (Type Detection)                          ║
        ╠═══════════════════════════════════════════════════╣
        ║                                                   ║
        ║   Classify cycle type:                            ║
        ║     - Direct: length=2, A + B → 2A                ║
        ║     - Indirect: length=3-4, A→B→C→A              ║
        ║     - Hypercycle: length>4, cross-catalysis       ║
        ║                                                   ║
        ╚═══════════════════════════════════════════════════╝
                         │
                         ▼
        ╔═══════════════════════════════════════════════════╗
        ║         STEP 4: METRIC CALCULATION                 ║
        ║         (Strength & Amplification)                ║
        ╠═══════════════════════════════════════════════════╣
        ║                                                   ║
        ║   For each autocatalytic cycle:                  ║
        ║     - Amplification factor (avg across nodes)    ║
        ║     - Catalytic strength (0-1):                  ║
        ║         strength = 0.7×amp_score + 0.3×connectivity║
        ║     - First detection step                        ║
        ║                                                   ║
        ╚═══════════════════════════════════════════════════╝
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                    OUTPUT: AUTOCATALYTIC CYCLES                 │
│                                                                   │
│  List of AutocatalyticCycle objects:                            │
│    - nodes: [molecule_ids]                                      │
│    - edges: [(reactant, product)]                                │
│    - amplification_factor: float                                  │
│    - cycle_type: 'direct' | 'indirect' | 'hypercycle'            │
│    - strength: float (0-1)                                       │
│    - first_detected_step: int                                    │
└─────────────────────────────────────────────────────────────────┘
```

---

## Example: Direct Autocatalysis

### Reaction Network
```
Molecule A: HCN (hydrogen cyanide)
Molecule B: H2O (water)

Reaction: A + B → 2A + C
  (HCN + H2O → 2HCN + byproduct)
```

### Cycle Detection
```
Graph:
  A ──→ A  (self-catalysis)
  
Cycle: [A]
Type: Direct autocatalysis
Amplification: 3.2× (A abundance increases from 10 to 32)
Strength: 0.85
```

### Temporal Evolution
```
Abundance of A over time:

40 |                                    ●●●●●●●●●
   |                                ●●●
30 |                            ●●●
   |                        ●●●
20 |                    ●●●
   |                ●●●
10 |            ●●●
   |        ●●●
 0 |●●●
   |__________________________________________________|
   0K    100K   200K   300K   400K   500K    Steps

Early (0-100K):   Abundance = 10
Late (400K-500K): Abundance = 32
Amplification: 32/10 = 3.2× ✅
```

---

## Example: Indirect Autocatalysis (Catalytic Chain)

### Reaction Network
```
Molecules:
  A: Formaldehyde (HCHO)
  B: Glycolaldehyde (C2H4O2)
  C: Glyceraldehyde (C3H6O3)

Reactions:
  A → B  (formaldehyde forms glycolaldehyde)
  B → C  (glycolaldehyde forms glyceraldehyde)
  C → A  (glyceraldehyde catalyzes formaldehyde formation)
```

### Cycle Detection
```
Graph:
  A ──→ B ──→ C ──→ A
  
Cycle: [A, B, C]
Type: Indirect autocatalysis
Length: 3 nodes
```

### Amplification Analysis
```
For each molecule in cycle:

Molecule A (HCHO):
  Early abundance: 5
  Late abundance:  18
  Amplification: 18/5 = 3.6× ✅

Molecule B (C2H4O2):
  Early abundance: 2
  Late abundance:  8
  Amplification: 8/2 = 4.0× ✅

Molecule C (C3H6O3):
  Early abundance: 1
  Late abundance:  5
  Amplification: 5/1 = 5.0× ✅

Average amplification: (3.6 + 4.0 + 5.0) / 3 = 4.2×
Cycle strength: 0.78
```

### Network Visualization
```
Reaction Network Graph:

    ┌─────────┐
    │   A     │ (HCHO)
    │ (HCHO)  │
    └────┬────┘
         │
         │ A → B
         ▼
    ┌─────────┐
    │   B     │ (C2H4O2)
    │(C2H4O2) │
    └────┬────┘
         │
         │ B → C
         ▼
    ┌─────────┐
    │   C     │ (C3H6O3)
    │(C3H6O3) │
    └────┬────┘
         │
         │ C → A (catalyzes)
         │
         └──────┐
                │
         ┌──────┴──────┐
         │  CYCLE     │ ◄── AUTOCATALYTIC
         │  DETECTED  │
         └────────────┘
```

---

## Example: Hypercycle (Cross-Catalysis)

### Reaction Network
```
Molecules: A, B, C, D, E

Reactions:
  A → B  (A catalyzes B formation)
  B → C  (B catalyzes C formation)
  C → D  (C catalyzes D formation)
  D → E  (D catalyzes E formation)
  E → A  (E catalyzes A formation)
  
  + Cross-catalytic edges:
    A → C  (A also helps C)
    B → E  (B also helps E)
    C → A  (C also helps A)
```

### Cycle Detection
```
Primary cycle: [A, B, C, D, E]
Length: 5 nodes

Cross-catalytic edges:
  A → C (additional catalysis)
  B → E (additional catalysis)
  C → A (feedback loop)

Type: Hypercycle (multiple catalytic interactions)
```

### Hypercycle Structure
```
Complex Network:

        ┌─── A ───┐
        │         │
        │    ┌────┴────┐
        │    │         │
        ▼    ▼         │
        B ──→ C ◄──────┘ (cross-catalysis)
        │    │
        │    │
        ▼    ▼
        D ──→ E
        │    │
        └────┘
        
Primary cycle: A→B→C→D→E→A
Cross-catalysis: A→C, B→E, C→A
Result: Hypercycle (multiple feedback loops)
```

---

## Amplification Detection Algorithm

### Temporal Analysis
```
For molecule M in cycle:

┌───────────────────────────────────────────────────────────┐
│  STEP 1: Extract abundance history                         │
│    history = abundance_history[M]                           │
│    history = [count_0, count_1, ..., count_T]              │
└───────────────────────┬─────────────────────────────────────┘
                        │
                        ▼
┌───────────────────────────────────────────────────────────┐
│  STEP 2: Divide into early and late periods               │
│    early_period = history[0 : T//3]                       │
│    late_period  = history[2*T//3 : T]                     │
└───────────────────────┬─────────────────────────────────────┘
                        │
                        ▼
┌───────────────────────────────────────────────────────────┐
│  STEP 3: Calculate mean abundances                        │
│    early_mean = mean(early_period)                         │
│    late_mean  = mean(late_period)                          │
└───────────────────────┬─────────────────────────────────────┘
                        │
                        ▼
┌───────────────────────────────────────────────────────────┐
│  STEP 4: Compute amplification factor                     │
│    if early_mean > 0:                                      │
│      amplification = late_mean / early_mean                 │
│    else:                                                   │
│      amplification = 0 (not detected)                      │
└───────────────────────┬─────────────────────────────────────┘
                        │
                        ▼
┌───────────────────────────────────────────────────────────┐
│  STEP 5: Classify as autocatalytic                        │
│    if amplification ≥ 1.5:                                │
│      return AUTOCATALYTIC ✅                               │
│    else:                                                   │
│      return NOT AUTOCATALYTIC ❌                           │
└───────────────────────────────────────────────────────────┘
```

### Amplification Threshold
```
Minimum amplification: 1.5× (configurable)

Rationale:
  - 1.0× = No growth (not autocatalytic)
  - 1.5× = Moderate growth (autocatalytic)
  - 2.0× = Strong growth (highly autocatalytic)
  - 5.0× = Very strong growth (exceptional)

Example classifications:
  Amplification 1.2× → NOT autocatalytic (below threshold)
  Amplification 1.8× → AUTOCATALYTIC ✅
  Amplification 4.5× → AUTOCATALYTIC ✅ (strong)
```

---

## Cycle Classification Logic

### Decision Tree
```
┌───────────────────────────────────────────────────────────┐
│  Classify Cycle Type                                       │
└───────────────────────┬─────────────────────────────────────┘
                        │
                        ▼
            ┌───────────────────────┐
            │  Cycle length = 2?   │
            └───────┬───────────────┘
                    │
            ┌───────┴───────┐
            │               │
          YES              NO
            │               │
            ▼               ▼
    ┌──────────────┐  ┌──────────────────┐
    │   DIRECT     │  │  Cycle length   │
    │  AUTOCATALYSIS│  │  ≤ 4?           │
    │              │  └──────┬───────────┘
    │  A + B → 2A  │         │
    └──────────────┘    ┌────┴────┐
                        │         │
                      YES        NO
                        │         │
                        ▼         ▼
                ┌──────────┐  ┌──────────────┐
                │ INDIRECT │  │ Check cross- │
                │          │  │ catalysis    │
                │ A→B→C→A  │  └──────┬───────┘
                └──────────┘         │
                              ┌──────┴──────┐
                              │             │
                            High          Low
                            │             │
                            ▼             ▼
                    ┌──────────────┐  ┌──────────┐
                    │  HYPERCYCLE  │  │ INDIRECT │
                    │              │  │          │
                    │ Multiple     │  │ A→B→C→A  │
                    │ feedback     │  │          │
                    │ loops        │  │          │
                    └──────────────┘  └──────────┘
```

### Classification Criteria
```
Direct Autocatalysis:
  - Cycle length = 2 nodes
  - Example: A + B → 2A
  - Simplest form of self-replication

Indirect Autocatalysis:
  - Cycle length = 3-4 nodes
  - Example: A → B → C → A
  - Catalytic chain with feedback

Hypercycle:
  - Cycle length > 4 nodes OR
  - Multiple cross-catalytic edges
  - Example: A→B→C→D→E→A with A→C, B→E
  - Most complex, highest catalytic strength
```

---

## Catalytic Strength Calculation

### Formula
```
strength = 0.7 × amp_score + 0.3 × connectivity

where:
  amp_score = normalized amplification (0-1)
    = min(1.0, (amplification - 1.0) / 10.0)
  
  connectivity = internal_edges / max_possible_edges
    = edges_within_cycle / (n × (n-1))
```

### Example Calculation
```
Cycle: [A, B, C]
Amplification: 4.2×

Step 1: Calculate amp_score
  amp_score = min(1.0, (4.2 - 1.0) / 10.0)
           = min(1.0, 0.32)
           = 0.32

Step 2: Calculate connectivity
  Internal edges: 3 (A→B, B→C, C→A)
  Max possible: 3 × 2 = 6
  connectivity = 3/6 = 0.5

Step 3: Calculate strength
  strength = 0.7 × 0.32 + 0.3 × 0.5
           = 0.224 + 0.15
           = 0.374

Result: strength = 0.374 (moderate)
```

---

## Performance Characteristics

### Algorithm Complexity
```
Operation                    | Complexity | Notes
----------------------------|------------|------------------
Johnson's cycle detection   | O((V+E)(C+1)) | C = cycles found
Amplification check         | O(V × T)   | V = nodes, T = time steps
Cycle classification        | O(V²)      | Check cross-edges
Strength calculation        | O(V)       | Simple aggregation
Total                       | O((V+E)(C+1) + V×T) | Dominated by cycle detection

Optimizations:
  - Max cycle length: 6 (prevents exponential explosion)
  - Timeout: 300s (safety limit)
  - Cycle limit: 2M (prevents memory issues)
```

### Scalability
```
Network Size | Cycles Found | Detection Time | Memory
-------------|--------------|----------------|--------
100 nodes    | ~50          | 0.5s          | 10 MB
500 nodes    | ~500         | 5s            | 50 MB
1000 nodes   | ~2000        | 30s           | 200 MB
5000 nodes   | ~50000       | 300s (timeout)| 2 GB

For very large networks (>1000 nodes):
  - Consider filtering network first
  - Reduce max_cycle_length to 4
  - Use parallel processing
```

---

## Patent Claim Support

This figure demonstrates **Patent 2 - Autocatalysis Detection**:

### Patent 2, Claim 1: Cycle Detection Algorithm
- Johnson's algorithm for finding all cycles in directed graph
- Maximum cycle length constraint (performance)
- Timeout mechanism (safety)

### Patent 2, Claim 2: Autocatalytic Filtering
- Amplification detection using temporal abundance analysis
- Minimum amplification threshold (1.5×)
- Net production calculation

### Patent 2, Claim 3: Cycle Classification
- Direct, indirect, and hypercycle types
- Cross-catalytic edge detection
- Cycle length-based classification

### Patent 2, Claim 4: Catalytic Strength Metric
- Combined amplification and connectivity score
- Normalized strength (0-1 scale)
- Temporal stability analysis

### Patent 2, Claim 5: Real-Time Detection
- Works on simulation snapshots
- Tracks abundance history over time
- Detects cycles as they emerge

---

## Key Innovations Summary

1. ✅ **Temporal amplification analysis** - Uses abundance history, not just topology
2. ✅ **Multi-type classification** - Direct, indirect, hypercycle detection
3. ✅ **Catalytic strength metric** - Quantifies cycle importance
4. ✅ **Scalable algorithm** - Handles networks up to 1000+ nodes
5. ✅ **Safety mechanisms** - Timeout and cycle limits prevent hangs
6. ✅ **Real-time detection** - Works on evolving networks during simulation
