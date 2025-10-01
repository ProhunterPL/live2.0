# Bond System Enhancement Plan

## Current Implementation vs. Suggested Model

### Currently Have ✅
- `bond_age` - bond lifetime tracking
- `bond_energy` - bond energy
- `bond_matrix` - bond strength
- `bond_active` - active bond flag
- Distance-based formation
- Distance & strength-based breaking
- Basic cluster detection
- Diagnostics logging of bonds

### Missing Features from Suggested Model 🎯

#### 1. Bond Types & Parameters
```python
# Currently: single type, no differentiation
# Proposed: multiple bond types with distinct physics

bond_type: int         # 0=van der Waals, 1=covalent, 2=H-bond, 3=metallic
k_spring: float        # spring constant per bond
rest_len: float        # rest length per bond  
damping: float         # damping coefficient
strength: float        # breaking threshold
flags: int             # bitmask: directed, locked, catalytic
```

#### 2. Valence System
```python
# Currently: unlimited bonds per particle
# Proposed: valence constraints

valence_max: int       # max bonds per particle
bonds: set[int]        # efficient bond lookup per particle
```

#### 3. Advanced Formation Rules
```python
# Currently: simple distance + basic compatibility
# Proposed:

✓ Distance: dist(a,b) < r_capture(typ_a, typ_b)
✓ Energy: energy_a + energy_b > E_min
✓ Valence: len(bonds[a]) < valence_max[a]
✓ Type compatibility: allowed[(typ_a, typ_b)] -> bond_type
✓ Probabilistic: p = sigmoid((E - E_min)/sigma) * f(dist)
```

#### 4. Advanced Breaking Rules
```python
# Currently: distance threshold only
# Proposed multiple conditions:

✓ Overload: |F_spring| > strength OR |dist-rest_len|/rest_len > eps
✓ Aging: age > max_age with probability
✓ Energy deficit: energy_a + energy_b < E_keep
✓ Valence conflict: stronger bond wants to form
```

#### 5. Physical Forces
```python
# Currently: no direct spring forces
# Proposed:

F_spring = -k*(d - rest_len)*dir - c*(v_rel·dir)*dir
Energy_diffusion = α*(E_mean_neighbors - E_i)
Angular_constraints (for rigid structures)
```

#### 6. Event System
```python
# Currently: implicit in diagnostics
# Proposed: explicit callbacks

on_bond_created(a, b, type, step, reason)
on_bond_broken(a, b, reason, step)
```

## Implementation Priority

### Phase 1: Core Enhancements (High Priority) 🔥
1. **Bond types** - Add `bond_type_field[i,j]`
2. **Rest length** - Add `bond_rest_len_field[i,j]`
3. **Spring constant** - Add `bond_k_spring_field[i,j]`
4. **Strength threshold** - Use for breaking condition
5. **Spring forces** - Implement F = -k*(d - rest_len)*dir

### Phase 2: Valence System (Medium Priority) ⚡
1. **Valence max** - Add per-particle valence limit
2. **Bond count** - Track bonds per particle efficiently
3. **Valence-based formation** - Respect max bonds
4. **Bond competition** - Handle valence conflicts

### Phase 3: Advanced Physics (Medium Priority) ⚙️
1. **Damping** - Add velocity-dependent damping
2. **Energy diffusion** - Spread energy along bonds
3. **Probabilistic formation** - Sigmoid probability
4. **Aging breakage** - Probability-based aging

### Phase 4: Events & Analysis (Low Priority) 📊
1. **Event callbacks** - Structured event logging
2. **Bond type analysis** - Track type distributions
3. **Force analysis** - Track overload events
4. **Valence saturation metrics** - Track valence usage

## Proposed Architecture

### Data Structures (Taichi fields)

```python
# Per-bond properties (N×N matrices)
bond_active[i,j]: i32         # existing
bond_type[i,j]: i32           # NEW: 0=vdW, 1=covalent, 2=H, 3=metal
bond_strength[i,j]: f32       # existing (bond_matrix)
bond_rest_len[i,j]: f32       # NEW
bond_k_spring[i,j]: f32       # NEW
bond_damping[i,j]: f32        # NEW
bond_age[i,j]: f32            # existing
bond_energy[i,j]: f32         # existing
bond_flags[i,j]: i32          # NEW: bitmask

# Per-particle bond management
particle_valence_max[i]: i32  # NEW
particle_bond_count[i]: i32   # NEW
```

### Configuration

```python
class BondingConfig:
    # Bond type parameters
    bond_types: Dict[int, BondTypeParams]
    
    # Formation rules
    r_capture: Dict[tuple, float]  # (type_a, type_b) -> capture radius
    E_min_formation: float
    allowed_types: Dict[tuple, int]  # (type_a, type_b) -> bond_type
    
    # Breaking rules
    overload_threshold: float  # strain threshold
    max_age: Dict[int, float]  # by bond type
    E_keep: float  # min energy to maintain bond
    
    # Physics
    enable_spring_forces: bool
    enable_damping: bool
    enable_energy_diffusion: bool

class BondTypeParams:
    id: int
    name: str
    k_spring: float
    rest_len: float
    damping: float
    strength: float
    max_age: float
    color: tuple  # for visualization
```

## Benefits

1. **Physical Realism**
   - Realistic spring forces with damping
   - Type-specific bond behavior
   - Energy conservation via diffusion

2. **Chemical Accuracy**
   - Valence limits (like real atoms)
   - Different bond types (weak/strong)
   - Formation/breaking rules based on chemistry

3. **Emergent Complexity**
   - Complex molecules from simple rules
   - Bond competition → interesting dynamics
   - Aging → dynamic equilibrium

4. **Analysis & Debugging**
   - Track bond type distributions
   - Identify overload events
   - Monitor valence saturation
   - Spring force analysis

## Backward Compatibility

- All enhancements are **additive**
- Default configuration matches current behavior
- Can enable features incrementally
- Diagnostics system already compatible

## Next Steps

1. ✅ Review this plan
2. Implement Phase 1 (bond types + spring forces)
3. Add configuration options
4. Update diagnostics to track new metrics
5. Add visualization for bond types
6. Write tests for new features
7. Document new capabilities

## Related Files

- `backend/sim/core/binding.py` - core bond logic
- `backend/sim/core/diagnostics.py` - metrics logging
- `backend/sim/config.py` - configuration
- `backend/sim/core/potentials.py` - force calculations
- `docs/DIAGNOSTICS.md` - analysis documentation

