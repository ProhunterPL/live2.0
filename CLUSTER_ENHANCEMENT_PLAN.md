# Cluster Detection & Management Enhancement Plan

## Current Implementation vs. Suggested Model

### Currently Have ✅
- Basic cluster detection (simple flood fill)
- `cluster_id` and `cluster_sizes` fields
- `find_cluster_root()` and `union_clusters()` functions (NOT USED!)
- Basic cluster stats (count, max size, average)
- Diagnostics logging (size histogram, largest cluster)

### Critical Issues ⚠️
1. **Union-Find NOT used** - has DSU code but uses slow flood fill
2. **No cluster metrics** - missing R_g, density, graph metrics
3. **No lifecycle tracking** - birth, split, merge, dissolve events not tracked
4. **No density filtering** - accepts any connected component
5. **No spatial merging** - close components not merged
6. **No event system** - no callbacks for cluster events

## Proposed Enhancements

### Phase 1: Proper Union-Find (High Priority) 🔥

**Problem:** Current `update_clusters_kernel()` uses O(N²) flood fill  
**Solution:** Use existing DSU code with GPU-friendly implementation

```python
@ti.kernel
def detect_clusters_dsu_kernel(active: ti.template(), particle_count: ti.i32):
    # Initialize: each particle is its own cluster
    for i in range(particle_count):
        if active[i] == 1:
            cluster_id_field[i] = i
        else:
            cluster_id_field[i] = -1
    
    # Union phase: process all bonds
    for i, j in ti.ndrange(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE):
        if bond_active_field[i, j] == 1:
            # Union the clusters
            root_i = find_root(i)
            root_j = find_root(j)
            if root_i != root_j:
                # Attach smaller to larger (union by rank)
                cluster_id_field[root_j] = root_i
    
    # Path compression pass
    for i in range(particle_count):
        if active[i] == 1:
            root = find_root(i)
            cluster_id_field[i] = root
```

**Benefits:**
- O(N*α(N)) instead of O(N²) - much faster
- GPU-friendly (2-3 parallel passes)
- Already have helper functions!

### Phase 2: Cluster Metrics (High Priority) 🔥

Add comprehensive cluster metrics for analysis and diagnostics:

```python
# Per-cluster metrics (tracked in separate fields or computed on-demand)
cluster_size[c]: i32              # |C| - already have
cluster_bond_count[c]: i32        # number of bonds in cluster
cluster_avg_degree[c]: f32        # 2*|Bonds(C)| / |C|
cluster_graph_density[c]: f32     # 2*|Bonds(C)| / (|C|*(|C|-1))
cluster_R_g[c]: f32               # radius of gyration (spatial extent)
cluster_centroid[c]: vec2         # center of mass
cluster_total_energy[c]: f32      # sum of particle energies
cluster_avg_energy[c]: f32        # average energy per particle
cluster_compactness[c]: f32       # avg distance to centroid
cluster_age[c]: i32               # steps since birth
cluster_stability[c]: f32         # rolling avg of bond changes
```

**Radius of Gyration (R_g):**
```python
@ti.func
def compute_R_g(cluster_members: ti.template(), count: ti.i32, 
                positions: ti.template()) -> ti.f32:
    # Compute center of mass
    com = ti.Vector([0.0, 0.0])
    for k in range(count):
        i = cluster_members[k]
        com += positions[i]
    com /= count
    
    # Compute R_g²
    R_g_sq = 0.0
    for k in range(count):
        i = cluster_members[k]
        r_vec = positions[i] - com
        R_g_sq += r_vec.dot(r_vec)
    R_g_sq /= count
    
    return ti.sqrt(R_g_sq)
```

### Phase 3: Density Filtering (Medium Priority) ⚡

Filter out low-density artifacts (random edges):

```python
# Configuration
min_density: float = 0.1  # minimum graph density
min_size: int = 4         # minimum cluster size

@ti.kernel
def filter_clusters_by_density(particle_count: ti.i32):
    # For each cluster, compute density and mark invalid ones
    for c in range(MAX_PARTICLES_COMPILE):
        if cluster_sizes_field[c] > 0:
            size = cluster_sizes_field[c]
            bonds = count_cluster_bonds(c)
            
            # Graph density = 2*E / (N*(N-1))
            if size >= 2:
                max_bonds = size * (size - 1) / 2
                density = ti.cast(bonds, ti.f32) / max_bonds
                
                if density < min_density or size < min_size:
                    # Mark cluster as invalid
                    cluster_valid_field[c] = 0
                else:
                    cluster_valid_field[c] = 1
```

### Phase 4: Spatial Merging (Medium Priority) ⚡

Merge close components (spatial bridges):

```python
# Configuration
d_merge: float = 5.0      # max centroid distance
r_prebond: float = 2.0    # max particle pair distance

def merge_close_clusters(clusters: List, positions: np.ndarray):
    # Compute centroids
    centroids = {}
    for c_id, members in clusters.items():
        centroid = np.mean([positions[i] for i in members], axis=0)
        centroids[c_id] = centroid
    
    # Find pairs to merge
    merges = []
    for c1 in clusters:
        for c2 in clusters:
            if c1 < c2:
                dist = np.linalg.norm(centroids[c1] - centroids[c2])
                if dist < d_merge:
                    # Check if any particle pair is close
                    if has_close_pair(clusters[c1], clusters[c2], positions, r_prebond):
                        merges.append((c1, c2))
    
    # Apply merges
    for c1, c2 in merges:
        union_clusters_py(c1, c2)
```

### Phase 5: Lifecycle Tracking (Medium Priority) ⚡

Track cluster birth, growth, split, merge, dissolution:

```python
class ClusterLifecycle:
    def __init__(self):
        self.clusters: Dict[int, ClusterInfo] = {}
        self.history: List[ClusterEvent] = []
    
    def update(self, current_clusters: Dict[int, set], step: int):
        # Detect new clusters (birth)
        new_ids = set(current_clusters.keys()) - set(self.clusters.keys())
        for c_id in new_ids:
            self.on_cluster_born(c_id, current_clusters[c_id], step)
        
        # Detect disappeared clusters (dissolved or split)
        old_ids = set(self.clusters.keys()) - set(current_clusters.keys())
        for c_id in old_ids:
            self.on_cluster_dissolved(c_id, step)
        
        # Detect splits (old cluster → multiple new)
        self.detect_splits(current_clusters, step)
        
        # Detect merges (multiple old → one new)
        self.detect_merges(current_clusters, step)
        
        # Update existing clusters
        for c_id, members in current_clusters.items():
            if c_id in self.clusters:
                self.clusters[c_id].update(members, step)
    
    def on_cluster_born(self, c_id: int, members: set, step: int):
        info = ClusterInfo(c_id, members, step)
        self.clusters[c_id] = info
        self.history.append(ClusterEvent('birth', c_id, step, len(members)))
    
    def on_cluster_dissolved(self, c_id: int, step: int):
        if c_id in self.clusters:
            lifetime = step - self.clusters[c_id].birth_step
            self.history.append(ClusterEvent('dissolve', c_id, step, lifetime=lifetime))
            del self.clusters[c_id]
    
    def detect_splits(self, current_clusters: Dict, step: int):
        # Complex: track particle membership changes
        # If cluster C splits into C1, C2: most particles from C go to C1 or C2
        pass
    
    def detect_merges(self, current_clusters: Dict, step: int):
        # Complex: track when two old clusters become one new
        pass

class ClusterInfo:
    def __init__(self, id: int, members: set, birth_step: int):
        self.id = id
        self.members = members
        self.birth_step = birth_step
        self.size_history = [len(members)]
        self.bond_change_history = []
    
    def update(self, new_members: set, step: int):
        self.size_history.append(len(new_members))
        # Track stability (bond changes)
        added = new_members - self.members
        removed = self.members - new_members
        self.bond_change_history.append(len(added) + len(removed))
        self.members = new_members

class ClusterEvent:
    def __init__(self, event_type: str, cluster_id: int, step: int, **kwargs):
        self.type = event_type  # 'birth', 'split', 'merge', 'dissolve'
        self.cluster_id = cluster_id
        self.step = step
        self.metadata = kwargs
```

### Phase 6: Event System (Low Priority) 📊

Callbacks for cluster events:

```python
# Event callbacks
def on_cluster_born(id: int, size: int, metrics: dict, step: int):
    log_event('cluster_birth', id, size, metrics, step)

def on_cluster_split(parent_id: int, child_ids: list, step: int):
    log_event('cluster_split', parent_id, len(child_ids), step)

def on_cluster_merge(id_a: int, id_b: int, id_new: int, step: int):
    log_event('cluster_merge', [id_a, id_b], id_new, step)

def on_cluster_dissolved(id: int, lifetime: int, step: int):
    log_event('cluster_dissolve', id, lifetime, step)
```

## Proposed Architecture

### Data Structures

```python
# Per-particle (existing)
cluster_id[i]: i32                    # cluster ID (root after DSU)
active[i]: i32                        # particle active flag

# Per-cluster (new)
cluster_sizes[c]: i32                 # existing
cluster_bond_count[c]: i32            # NEW
cluster_R_g[c]: f32                   # NEW - radius of gyration
cluster_centroid_x[c]: f32            # NEW
cluster_centroid_y[c]: f32            # NEW
cluster_total_energy[c]: f32          # NEW
cluster_avg_degree[c]: f32            # NEW
cluster_density[c]: f32               # NEW - graph density
cluster_age[c]: i32                   # NEW - steps since birth
cluster_valid[c]: i32                 # NEW - passes density filter
cluster_stability[c]: f32             # NEW - rolling avg stability

# Lifecycle tracking (Python-side)
cluster_lifecycle: ClusterLifecycle   # tracks events
cluster_history: List[ClusterEvent]   # event log
```

### Configuration

```python
class ClusterConfig:
    # Detection
    min_cluster_size: int = 4
    min_graph_density: float = 0.1
    use_dsu: bool = True              # use Union-Find vs flood fill
    
    # Spatial merging
    enable_spatial_merge: bool = False
    merge_distance: float = 5.0       # centroid distance
    merge_particle_radius: float = 2.0
    
    # Lifecycle
    track_lifecycle: bool = True
    dissolve_grace_period: int = 10   # steps before truly dissolved
    
    # Metrics
    compute_R_g: bool = True
    compute_graph_metrics: bool = True
    compute_stability: bool = True
```

## Integration with Diagnostics

Enhance diagnostics to include new cluster metrics:

```python
# In diagnostics.py
def _compute_cluster_metrics(self, positions, clusters, energies):
    metrics = {
        'num_clusters': len(clusters),
        'size_histogram': defaultdict(int),
        'largest_cluster_size': 0,
        'R_g_distribution': [],        # NEW
        'density_distribution': [],    # NEW
        'avg_degree_distribution': [], # NEW
        'energy_per_cluster': [],      # NEW
        'compactness_distribution': [] # NEW
    }
    
    for cluster_id, members in clusters.items():
        size = len(members)
        metrics['size_histogram'][size] += 1
        
        if size > metrics['largest_cluster_size']:
            metrics['largest_cluster_size'] = size
        
        # Compute R_g
        R_g = self._compute_cluster_R_g(members, positions)
        metrics['R_g_distribution'].append(R_g)
        
        # Compute graph density
        bonds_in_cluster = self._count_cluster_bonds(members)
        max_bonds = size * (size - 1) / 2
        density = bonds_in_cluster / max_bonds if max_bonds > 0 else 0
        metrics['density_distribution'].append(density)
        
        # Average degree
        avg_degree = 2 * bonds_in_cluster / size if size > 0 else 0
        metrics['avg_degree_distribution'].append(avg_degree)
    
    return metrics
```

## Implementation Steps

### Step 1: Fix Union-Find (1-2 hours) 🔥
1. Replace flood fill with proper DSU in `update_clusters_kernel`
2. Use existing `find_cluster_root` function
3. Add path compression pass
4. Benchmark performance improvement

### Step 2: Add Core Metrics (2-3 hours) 🔥
1. Add cluster metric fields
2. Implement R_g computation
3. Implement graph density
4. Implement average degree
5. Update diagnostics to log new metrics

### Step 3: Density Filtering (1 hour) ⚡
1. Add min_density configuration
2. Implement filtering kernel
3. Mark invalid clusters
4. Update get_clusters to respect filter

### Step 4: Lifecycle Tracking (3-4 hours) ⚡
1. Create ClusterLifecycle class
2. Track birth/dissolve events
3. Detect splits (complex)
4. Detect merges (complex)
5. Add to diagnostics CSV

### Step 5: Spatial Merging (2-3 hours) ⚡
1. Implement centroid computation
2. Implement close-pair detection
3. Implement merge logic
4. Add configuration options

### Step 6: Event System (1-2 hours) 📊
1. Define event callback interface
2. Integrate with lifecycle
3. Add event logging to diagnostics
4. Create visualization helpers

## Main Loop Integration

```python
def sim_step():
    # 1. Physics
    integrate_forces()
    
    # 2. Bonds (from Bond Enhancement Plan)
    try_form_bonds()
    apply_bond_forces_and_energy()
    break_weak_bonds()
    
    # 3. Clusters (NEW ENHANCED)
    detect_clusters()              # DSU-based, fast
    compute_cluster_metrics()      # R_g, density, etc.
    filter_by_density()            # remove artifacts
    merge_close_clusters()         # optional spatial bridges
    update_cluster_lifecycle()     # track events
    
    # 4. Visualization & Logging
    update_visuals()
    log_observables()
```

## Sensible Defaults

```python
min_cluster_size = 4              # minimum 4 particles
min_graph_density = 0.15          # 15% of max possible bonds
dissolve_grace_period = 25        # 25 steps tolerance
merge_distance = 5.0              # 5 particle radii
merge_particle_radius = 2.0       # 2 particle radii
```

## Benefits

1. **Performance** 
   - O(N*α(N)) vs O(N²) - massive speedup
   - GPU-friendly parallel DSU
   - Better scaling to large systems

2. **Physical Accuracy**
   - R_g tracks spatial extent (important for polymers)
   - Graph density distinguishes tight vs loose clusters
   - Stability metric shows dynamic equilibrium

3. **Scientific Insight**
   - Lifecycle tracking reveals phase transitions
   - Split/merge events show reorganization
   - Density filtering removes noise

4. **Analysis Power**
   - Rich metrics for diagnostics
   - Event history for trajectory analysis
   - Better visualization (color by density, R_g, etc.)

## Visualization Enhancements

```python
# Cluster rendering
for cluster in clusters:
    # Convex hull or "blob"
    hull = compute_convex_hull(cluster.members)
    draw_hull(hull, color=cluster_color(cluster.density))
    
    # Label
    label = f"#{cluster.id}: n={cluster.size}, ρ={cluster.density:.2f}, R_g={cluster.R_g:.1f}"
    draw_text(cluster.centroid, label)
    
    # Bonds (thickness ~ tension)
    for bond in cluster.bonds:
        thickness = compute_bond_tension(bond)
        color = bond_type_color(bond.type)
        alpha = 1.0 - bond.age / max_age  # fade with age
        draw_line(bond.a, bond.b, thickness, color, alpha)
```

## Next Steps

1. ✅ Review this plan
2. Implement Step 1 (DSU) - immediate performance gain
3. Implement Step 2 (metrics) - rich analysis
4. Add configuration options
5. Update diagnostics system
6. Create visualization demo
7. Document new features

## Related Files

- `backend/sim/core/binding.py` - cluster detection
- `backend/sim/core/diagnostics.py` - metrics logging
- `backend/sim/config.py` - configuration
- `docs/DIAGNOSTICS.md` - analysis documentation
- `BOND_ENHANCEMENT_PLAN.md` - related bond improvements

