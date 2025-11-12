"""
Binding system for Live 2.0 simulation
Handles particle binding, bond formation/breaking, and cluster detection
"""

import taichi as ti
import numpy as np
from typing import List, Tuple, Dict, Set, Optional
from ..config import SimulationConfig

# Compile-time constants
MAX_PARTICLES_COMPILE = 10000
PARTICLE_RADIUS_COMPILE = 0.5

# Global Taichi fields
bond_matrix_field = None
bond_active_field = None
cluster_id_field = None
cluster_sizes_field = None
next_cluster_id_field = None
bond_energy_field = None
bond_age_field = None
# New bond type fields
bond_type_field = None
bond_k_spring_field = None
bond_rest_len_field = None
bond_damping_field = None
bond_strength_field = None
# Valence system fields
particle_valence_max_field = None
particle_bond_count_field = None
# Cluster metrics fields
cluster_bond_count_field = None
cluster_R_g_field = None
cluster_density_field = None

def init_binding_fields():
    """Initialize global Taichi fields for binding"""
    global bond_matrix_field, bond_active_field, cluster_id_field
    global cluster_sizes_field, next_cluster_id_field, bond_energy_field, bond_age_field
    global bond_type_field, bond_k_spring_field, bond_rest_len_field
    global bond_damping_field, bond_strength_field
    global particle_valence_max_field, particle_bond_count_field
    global cluster_bond_count_field, cluster_R_g_field, cluster_density_field
    
    bond_matrix_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    bond_active_field = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    cluster_id_field = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    cluster_sizes_field = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    next_cluster_id_field = ti.field(dtype=ti.i32, shape=())
    bond_energy_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    bond_age_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    # Bond type fields
    bond_type_field = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    bond_k_spring_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    bond_rest_len_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    bond_damping_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    bond_strength_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    # Valence system fields
    particle_valence_max_field = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    particle_bond_count_field = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    # Cluster metrics fields
    cluster_bond_count_field = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    cluster_R_g_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    cluster_density_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))

# Module-level kernels
@ti.kernel
def reset_binding_kernel():
    """Reset binding system - module-level kernel"""
    for i, j in ti.ndrange(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE):
        bond_matrix_field[i, j] = 0.0
        bond_active_field[i, j] = 0
        bond_energy_field[i, j] = 0.0
        bond_age_field[i, j] = 0.0
        bond_type_field[i, j] = 0
        bond_k_spring_field[i, j] = 0.0
        bond_rest_len_field[i, j] = 0.0
        bond_damping_field[i, j] = 0.0
        bond_strength_field[i, j] = 0.0
    
    for i in range(MAX_PARTICLES_COMPILE):
        cluster_id_field[i] = -1
        cluster_sizes_field[i] = 0
        particle_valence_max_field[i] = 20  # MASSIVE valence limit for huge clusters
        particle_bond_count_field[i] = 0
    
    next_cluster_id_field[None] = 0

@ti.kernel
def update_bond_counts_kernel(particle_count: ti.i32):
    """Update bond counts for valence system"""
    # Reset counts
    for i in range(MAX_PARTICLES_COMPILE):
        particle_bond_count_field[i] = 0
    
    # Count bonds per particle
    max_check = ti.min(particle_count, MAX_PARTICLES_COMPILE)
    for i in range(max_check):
        for j in range(i + 1, max_check):
            if bond_active_field[i, j] == 1:
                ti.atomic_add(particle_bond_count_field[i], 1)
                ti.atomic_add(particle_bond_count_field[j], 1)

@ti.kernel
def compute_cluster_metrics_kernel(positions: ti.template(), particle_count: ti.i32):
    """Compute cluster metrics: bond count, R_g, density"""
    # Reset cluster metrics
    for c in range(MAX_PARTICLES_COMPILE):
        cluster_bond_count_field[c] = 0
        cluster_R_g_field[c] = 0.0
        cluster_density_field[c] = 0.0
    
    # Count bonds per cluster
    max_check = ti.min(particle_count, MAX_PARTICLES_COMPILE)
    for i in range(max_check):
        for j in range(i + 1, max_check):
            if bond_active_field[i, j] == 1:
                cluster_i = cluster_id_field[i]
                cluster_j = cluster_id_field[j]
                if cluster_i == cluster_j and cluster_i >= 0:
                    ti.atomic_add(cluster_bond_count_field[cluster_i], 1)
    
    # Compute R_g and density for each cluster
    for c in range(MAX_PARTICLES_COMPILE):
        size = cluster_sizes_field[c]
        if size > 1:
            # Compute center of mass
            com_x = 0.0
            com_y = 0.0
            count = 0
            
            for i in range(max_check):
                if cluster_id_field[i] == c:
                    com_x += positions[i][0]
                    com_y += positions[i][1]
                    count += 1
            
            if count > 0:
                com_x /= count
                com_y /= count
                
                # Compute R_g²
                R_g_sq = 0.0
                for i in range(max_check):
                    if cluster_id_field[i] == c:
                        dx = positions[i][0] - com_x
                        dy = positions[i][1] - com_y
                        R_g_sq += dx * dx + dy * dy
                
                R_g_sq /= count
                cluster_R_g_field[c] = ti.sqrt(R_g_sq)
                
                # Compute graph density
                bonds = cluster_bond_count_field[c]
                max_bonds = size * (size - 1) / 2
                if max_bonds > 0:
                    cluster_density_field[c] = ti.cast(bonds, ti.f32) / max_bonds

@ti.kernel
def update_bonds_kernel(positions: ti.template(), attributes: ti.template(),
                       active: ti.template(), particle_count: ti.i32, dt: ti.f32):
    """Update bond formation and breaking - module-level kernel"""
    # Update bond ages
    for i, j in ti.ndrange(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE):
        if bond_active_field[i, j] == 1:
            bond_age_field[i, j] += dt
    
    # Check for new bonds
    for i in range(particle_count):
        if active[i] == 1:
            for j in range(i + 1, particle_count):
                if active[j] == 1 and bond_active_field[i, j] == 0:
                    # Check if particles should form a bond
                    bond_type = should_form_bond_func(i, j, positions, attributes)
                    if bond_type >= 0:  # -1 means no bond, >= 0 is bond type
                        # Calculate current distance for rest length
                        r_vec = positions[i] - positions[j]
                        dist = r_vec.norm()
                        
                        # Initialize variables (required for Taichi kernels)
                        k_spring = 5.0  # default value
                        rest_len = dist
                        damping = 0.15
                        strength = 10.0
                        
                        # Set parameters based on bond type
                        # Type 0: vdW (weak), Type 1: covalent (strong), Type 2: H-bond, Type 3: metallic
                        # FIXED: Use realistic equilibrium distances, not current stretched distance
                        if bond_type == 0:  # van der Waals - equilibrium ~3.0-3.5 Å = 1.5-1.75 units
                            k_spring = 5.0  # Increased from 2.0 for stronger restoring force
                            rest_len = PARTICLE_RADIUS_COMPILE * 3.0  # Realistic vdW equilibrium distance
                            damping = 0.15  # Increased from 0.1
                            strength = 5.0
                        elif bond_type == 1:  # covalent - equilibrium ~1.5-2.0 Å = 0.75-1.0 units
                            k_spring = 20.0  # Increased from 10.0 for much stronger covalent bonds
                            rest_len = PARTICLE_RADIUS_COMPILE * 2.0  # Realistic covalent bond length
                            damping = 0.25  # Increased from 0.2
                            strength = 20.0
                        elif bond_type == 2:  # H-bond - equilibrium ~2.5-3.0 Å = 1.25-1.5 units
                            k_spring = 8.0  # Increased from 5.0
                            rest_len = PARTICLE_RADIUS_COMPILE * 2.5  # Realistic H-bond distance
                            damping = 0.2  # Increased from 0.15
                            strength = 10.0
                        else:  # metallic or other
                            k_spring = 10.0  # Increased from 7.0
                            rest_len = PARTICLE_RADIUS_COMPILE * 2.5  # Realistic metallic bond distance
                            damping = 0.3  # Increased from 0.25
                            strength = 15.0
                        
                        # Check valence constraints BEFORE forming bond
                        valence_i = particle_bond_count_field[i]
                        valence_j = particle_bond_count_field[j]
                        max_valence_i = particle_valence_max_field[i]
                        max_valence_j = particle_valence_max_field[j]
                        
                        # Only form bond if both particles have valence available
                        if valence_i < max_valence_i and valence_j < max_valence_j:
                            # Form bond with type-specific parameters
                            bond_active_field[i, j] = 1
                            bond_active_field[j, i] = 1
                            bond_matrix_field[i, j] = 1.0
                            bond_matrix_field[j, i] = 1.0
                            bond_type_field[i, j] = bond_type
                            bond_type_field[j, i] = bond_type
                            bond_k_spring_field[i, j] = k_spring
                            bond_k_spring_field[j, i] = k_spring
                            bond_rest_len_field[i, j] = rest_len
                            bond_rest_len_field[j, i] = rest_len
                            bond_damping_field[i, j] = damping
                            bond_damping_field[j, i] = damping
                            bond_strength_field[i, j] = strength
                            bond_strength_field[j, i] = strength
                            bond_energy_field[i, j] = -1.0  # Negative binding energy
                            bond_energy_field[j, i] = -1.0
                            bond_age_field[i, j] = 0.0
                            bond_age_field[j, i] = 0.0
                            
                            # Update bond counts for valence system
                            ti.atomic_add(particle_bond_count_field[i], 1)
                            ti.atomic_add(particle_bond_count_field[j], 1)
    
    # Check for bond breaking
    for i, j in ti.ndrange(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE):
        if bond_active_field[i, j] == 1:
            if should_break_bond_func(i, j, positions, attributes):
                # Break bond inline
                bond_active_field[i, j] = 0
                bond_active_field[j, i] = 0
                bond_matrix_field[i, j] = 0.0
                bond_matrix_field[j, i] = 0.0
                
                # Update bond counts for valence system
                ti.atomic_sub(particle_bond_count_field[i], 1)
                ti.atomic_sub(particle_bond_count_field[j], 1)

@ti.func
def calculate_binding_probability(i: ti.i32, j: ti.i32, positions: ti.template(),
                                attributes: ti.template()) -> ti.f32:
    """Calculate binding probability based on particle compatibility"""
    pos_i = positions[i]
    pos_j = positions[j]
    
    # Calculate distance
    r_vec = pos_i - pos_j
    r = r_vec.norm()
    
    # Get particle properties
    mass_i = attributes[i][0]
    mass_j = attributes[j][0]
    charge_i = attributes[i][1]
    charge_j = attributes[j][1]
    
    # Calculate compatibility metrics
    mass_ratio = ti.min(mass_i, mass_j) / ti.max(mass_i, mass_j)
    charge_product = charge_i * charge_j
    charge_sum = ti.abs(charge_i + charge_j)
    
    # Base probability from distance (closer = higher probability)
    distance_factor = ti.exp(-r / 3.0)  # LESS steep decay - more permissive
    
    # Compatibility factor based on particle properties - MUCH MORE PERMISSIVE
    compatibility_factor = 0.0
    if mass_ratio > 0.5:  # REDUCED from 0.7 - more permissive
        compatibility_factor = 1.0  # INCREASED from 0.8 - maximum compatibility
    elif charge_product < -0.3:  # REDUCED from -0.5 - easier ionic bonds
        compatibility_factor = 0.8  # INCREASED from 0.6
    elif charge_sum > 0.3:  # REDUCED from 0.5 - easier metallic bonds
        compatibility_factor = 0.6  # INCREASED from 0.4
    elif mass_ratio > 0.1:  # REDUCED from 0.3 - very permissive vdW
        compatibility_factor = 0.4  # INCREASED from 0.2
    
    # Final probability - ensure minimum for close particles
    probability = distance_factor * compatibility_factor
    if r <= PARTICLE_RADIUS_COMPILE * 4.0:  # Very close particles get minimum probability
        probability = ti.max(probability, 0.15)  # Minimum 15% for very close particles
    return ti.min(probability, 1.0)  # Cap at 1.0

@ti.func
def should_form_bond_func(i: ti.i32, j: ti.i32, positions: ti.template(),
                         attributes: ti.template()) -> ti.i32:
    """Check if two particles should form a bond, returns bond type (-1 = no bond)"""
    pos_i = positions[i]
    pos_j = positions[j]
    
    # Calculate distance
    r_vec = pos_i - pos_j
    r = r_vec.norm()
    
    # Default: no bond
    bond_type = -1
    
    # FIXED: Bonds should only form when particles are close to equilibrium distance
    # This prevents the cycle of forming-stretching-breaking-reforming
    # vdW rest_len = 1.5, covalent = 1.0, H-bond = 1.25
    # Allow formation within 20% of rest length to ensure stability
    max_formation_dist_vdW = PARTICLE_RADIUS_COMPILE * 3.0 * 1.2  # 1.5 * 1.2 = 1.8 units
    max_formation_dist_covalent = PARTICLE_RADIUS_COMPILE * 2.0 * 1.2  # 1.0 * 1.2 = 1.2 units
    max_formation_dist_hbond = PARTICLE_RADIUS_COMPILE * 2.5 * 1.2  # 1.25 * 1.2 = 1.5 units
    
    # Calculate binding probability first
    binding_probability = calculate_binding_probability(i, j, positions, attributes)
    
    # Check distance thresholds based on potential bond type (determined by compatibility)
    # Get particle properties to determine bond type
    mass_i = attributes[i][0]
    mass_j = attributes[j][0]
    charge_i = attributes[i][1]
    charge_j = attributes[j][1]
    mass_ratio = ti.min(mass_i, mass_j) / ti.max(mass_i, mass_j)
    charge_product = charge_i * charge_j
    charge_sum = ti.abs(charge_i + charge_j)
    
    # Determine appropriate distance threshold based on compatibility
    max_formation_dist = max_formation_dist_vdW  # default
    if mass_ratio > 0.2:  # covalent candidate
        max_formation_dist = max_formation_dist_covalent
    elif charge_product < -0.1:  # H-bond candidate
        max_formation_dist = max_formation_dist_hbond
    
    if r <= max_formation_dist:
        # FIXED: Higher threshold prevents unstable bonds from forming too easily
        # Realistic bonds require stronger compatibility (0.15 instead of 0.005)
        if binding_probability > 0.15:  # INCREASED from 0.005 - more selective bond formation
            # Determine bond type based on compatibility (already calculated above)
            # SCIENTIFICALLY REALISTIC bonding conditions - balanced for stability
            if mass_ratio > 0.2:  # covalent bonds
                bond_type = 1
            elif charge_product < -0.1:  # H-bond
                bond_type = 2
            elif charge_sum > 0.1:  # metallic
                bond_type = 3
            elif mass_ratio > 0.05:  # vdW
                bond_type = 0
            else:
                bond_type = -1  # No bond if compatibility too low
    
    return bond_type

@ti.func
def should_break_bond_func(i: ti.i32, j: ti.i32, positions: ti.template(),
                          attributes: ti.template()) -> ti.i32:
    """Check if a bond should be broken (advanced: overload, aging, distance)"""
    pos_i = positions[i]
    pos_j = positions[j]
    
    # Calculate distance
    r_vec = pos_i - pos_j
    r = r_vec.norm()
    
    result = 0
    
    # Get bond parameters
    rest_len = bond_rest_len_field[i, j]
    strength = bond_strength_field[i, j]
    age = bond_age_field[i, j]
    bond_type = bond_type_field[i, j]
    
    # Condition 1: Overload - too much stretch/compression
    # FIXED: Realistic strain threshold - bonds break at 30-50% stretch, not 300%!
    strain = ti.abs(r - rest_len) / ti.max(rest_len, 0.1)
    max_strain = 0.5  # 50% strain threshold (realistic for most bonds)
    if bond_type == 0:  # vdW - more flexible
        max_strain = 0.8  # 80% strain
    elif bond_type == 1:  # covalent - less flexible
        max_strain = 0.3  # 30% strain
    elif bond_type == 2:  # H-bond - medium flexibility
        max_strain = 0.6  # 60% strain
    
    if strain > max_strain:
        result = 1
    
    # Condition 2: Distance too large (safety check)
    # FIXED: Break bonds if they exceed realistic maximum distances
    max_distance = PARTICLE_RADIUS_COMPILE * 4.0  # Reduced from 5.0
    if bond_type == 1:  # covalent bonds should be shorter
        max_distance = PARTICLE_RADIUS_COMPILE * 3.0
    if r > max_distance:
        result = 1
    
    # Condition 3: Aging - probabilistic breaking for old bonds
    # MUCH MORE STABLE: Bonds should last much longer for cluster growth
    max_age = 10000.0  # default (MUCH increased from 2000.0)
    if bond_type == 0:  # vdW - short lived but still stable
        max_age = 5000.0  # increased from 1000.0
    elif bond_type == 1:  # covalent - very long lived
        max_age = 20000.0  # increased from 5000.0
    elif bond_type == 2:  # H-bond - medium but stable
        max_age = 8000.0  # increased from 2000.0
    else:  # metallic
        max_age = 12000.0  # increased from 3000.0
    
    if age > max_age:
        # Probabilistic: higher age = higher chance to break
        # MUCH MORE CONSERVATIVE: break only if age > 2.0*max_age (was 1.5)
        if age > max_age * 2.0:
            result = 1
    
    return result

@ti.func
def find_root_dsu(i: ti.i32) -> ti.i32:
    """Find root with path compression (DSU)"""
    root = i
    # Find root
    while cluster_id_field[root] != root and cluster_id_field[root] >= 0:
        root = cluster_id_field[root]
    
    # Path compression
    current = i
    while current != root and cluster_id_field[current] >= 0:
        next_node = cluster_id_field[current]
        cluster_id_field[current] = root
        current = next_node
    
    return root

@ti.kernel
def compute_bond_forces_kernel(positions: ti.template(), velocities: ti.template(),
                               forces: ti.template(), active: ti.template(), 
                               particle_count: ti.i32):
    """Compute spring forces from bonds and add to force field"""
    # For each active bond, compute spring force
    for i in range(particle_count):
        if active[i] == 1:
            for j in range(i + 1, particle_count):
                if active[j] == 1 and bond_active_field[i, j] == 1:
                    # Get bond parameters
                    k_spring = bond_k_spring_field[i, j]
                    rest_len = bond_rest_len_field[i, j]
                    damping = bond_damping_field[i, j]
                    
                    # Calculate distance and direction
                    r_vec = positions[j] - positions[i]
                    r = r_vec.norm()
                    
                    if r > 1e-6:  # Avoid division by zero
                        dir = r_vec / r
                        
                        # Spring force: F = -k * (d - rest_len) * dir
                        spring_force_mag = k_spring * (r - rest_len)
                        spring_force = spring_force_mag * dir
                        
                        # Damping force: F = -c * (v_rel · dir) * dir
                        v_rel = velocities[j] - velocities[i]
                        v_rel_along_bond = v_rel.dot(dir)
                        damping_force = damping * v_rel_along_bond * dir
                        
                        # Total force
                        total_force = spring_force + damping_force
                        
                        # Apply to both particles (Newton's 3rd law)
                        ti.atomic_add(forces[i], total_force)
                        ti.atomic_sub(forces[j], total_force)

@ti.kernel
def update_clusters_kernel(active: ti.template(), particle_count: ti.i32):
    """Update cluster assignments using Union-Find (DSU) - O(N*α(N)) - OPTIMIZED"""
    # Phase 1: Initialize - each particle is its own cluster
    for i in range(MAX_PARTICLES_COMPILE):
        if i < particle_count and active[i] == 1:
            cluster_id_field[i] = i  # self-parent
        else:
            cluster_id_field[i] = -1  # inactive
    
    # Phase 2: Union - process all bonds (OPTIMIZED with bounds checking)
    max_check = ti.min(particle_count, MAX_PARTICLES_COMPILE)
    for i in range(max_check):
        if active[i] == 1:
            for j in range(i + 1, max_check):
                if active[j] == 1 and bond_active_field[i, j] == 1:
                    # Union the two clusters
                    root_i = find_root_dsu(i)
                    root_j = find_root_dsu(j)
                    if root_i != root_j:
                        # Attach smaller to larger (union by rank)
                        if root_i < root_j:
                            cluster_id_field[root_j] = root_i
                        else:
                            cluster_id_field[root_i] = root_j
    
    # Phase 3: Path compression pass - flatten all paths
    for i in range(max_check):
        if active[i] == 1:
            root = find_root_dsu(i)
            cluster_id_field[i] = root
    
    # Phase 4: Count cluster sizes
    for i in range(MAX_PARTICLES_COMPILE):
        cluster_sizes_field[i] = 0
    
    for i in range(max_check):
        if active[i] == 1:
            root = cluster_id_field[i]
            if root >= 0 and root < MAX_PARTICLES_COMPILE:
                ti.atomic_add(cluster_sizes_field[root], 1)

@ti.data_oriented
class BindingSystem:
    """Manages particle binding and bond formation"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.max_particles = config.max_particles
        
        # Initialize global fields if not already done
        if bond_matrix_field is None:
            init_binding_fields()
        
        # Use global fields
        self.bond_matrix = bond_matrix_field
        self.bond_active = bond_active_field
        self.cluster_id = cluster_id_field
        self.cluster_sizes = cluster_sizes_field
        self.next_cluster_id = next_cluster_id_field
        self.bond_energy = bond_energy_field
        self.bond_age = bond_age_field
        # Bond type fields
        self.bond_type = bond_type_field
        self.bond_k_spring = bond_k_spring_field
        self.bond_rest_len = bond_rest_len_field
        self.bond_damping = bond_damping_field
        self.bond_strength = bond_strength_field
        
        # Bond type parameters - SCIENTIFICALLY CALIBRATED from literature
        # Literature: C-C bond k=2255 kJ/(mol·Å²), D_e=348 kJ/mol (Luo 2007)
        # Using 1/4 of literature values for numerical stability (GROMACS/NAMD best practice)
        # 0 = van der Waals (weak)
        # 1 = covalent (strong)
        # 2 = hydrogen bond (medium)
        # 3 = metallic (medium-strong)
        self.bond_type_params = {
            0: {'k_spring': 2.0, 'rest_len': 1.0, 'damping': 0.1, 'strength': 5.0},     # vdW - unchanged
            1: {'k_spring': 500.0, 'rest_len': 0.8, 'damping': 0.2, 'strength': 100.0}, # covalent - was: 10, 20 -> 50x, 5x stronger (Luo 2007: 2255, 348)
            2: {'k_spring': 50.0, 'rest_len': 1.2, 'damping': 0.15, 'strength': 30.0},  # H-bond - was: 5, 10 -> 10x, 3x stronger
            3: {'k_spring': 100.0, 'rest_len': 0.9, 'damping': 0.25, 'strength': 50.0}  # metallic - was: 7, 15 -> 14x, 3.3x stronger
        }
        
        # Initialize
        self.reset()
    
    def reset(self):
        """Reset binding system"""
        reset_binding_kernel()
    
    def update_bonds(self, positions, attributes, active, particle_count: int, dt: float):
        """Update bond formation and breaking based on particle interactions"""
        update_bonds_kernel(positions, attributes, active, particle_count, dt)
        # Update bond counts for valence system
        update_bond_counts_kernel(particle_count)
    
    def apply_bond_forces(self, positions, velocities, forces, active, particle_count: int):
        """Compute and apply spring forces from bonds"""
        compute_bond_forces_kernel(positions, velocities, forces, active, particle_count)
    
    @ti.func
    def should_form_bond(self, i: ti.i32, j: ti.i32, positions: ti.template(),
                        attributes: ti.template()) -> ti.i32:
        """Check if two particles should form a bond"""
        pos_i = positions[i]
        pos_j = positions[j]
        
        # Calculate distance
        r_vec = pos_i - pos_j
        r = r_vec.norm()
        
        # NAPRAWIONE: UPROSZCZONE - tylko odległość!
        result = 0
        binding_distance = self.config.particle_radius * 3.0  # Tylko 3x promień
        if r <= binding_distance:
            result = 1  # Tworzymy wiązanie gdy cząstki są blisko!
        
        return result
    
    @ti.func
    def should_break_bond(self, i: ti.i32, j: ti.i32, positions: ti.template(),
                         attributes: ti.template()) -> ti.i32:
        """Check if a bond should break"""
        pos_i = positions[i]
        pos_j = positions[j]
        
        # Calculate distance
        r_vec = pos_i - pos_j
        r = r_vec.norm()
        
        # Break if too far apart
        result = 0
        if r > self.config.particle_radius * 4.0:
            result = 1
        else:
            # Break if bond is too weak
            bond_strength = self.bond_matrix[i, j]
            if bond_strength < self.config.unbinding_threshold:
                result = 1
        
        # Break if particles have incompatible properties
        if result == 0:
            mass_i = attributes[i][0]
            mass_j = attributes[j][0]
            charge_i = attributes[i][1]
            charge_j = attributes[j][1]
            
            mass_ratio = ti.min(mass_i, mass_j) / ti.max(mass_i, mass_j)
            if mass_ratio < 0.3:  # Masses became too different
                result = 1
            
            if result == 0:
                charge_product = charge_i * charge_j
                if charge_product > 0.2:  # Charges became incompatible
                    result = 1
        
        return result
    
    @ti.func
    def compute_binding_energy(self, i: ti.i32, j: ti.i32, r: ti.f32,
                             attributes: ti.template()) -> ti.f32:
        """Compute binding energy between two particles"""
        mass_i = attributes[i][0]
        mass_j = attributes[j][0]
        charge_i = attributes[i][1]
        charge_j = attributes[j][1]
        
        # Lennard-Jones contribution
        lj_energy = self.lennard_jones_potential(r, 1.0, 1.0)
        
        # Coulomb contribution
        coulomb_energy = self.coulomb_potential(r, charge_i, charge_j, 1.0)
        
        # Binding contribution
        binding_energy = -1.0 * ti.exp(-r / 1.0)  # Attractive at short range
        
        return lj_energy + coulomb_energy + binding_energy
    
    @ti.func
    def lennard_jones_potential(self, r: ti.f32, epsilon: ti.f32, sigma: ti.f32) -> ti.f32:
        """Lennard-Jones potential"""
        if r < 0.1:
            r = 0.1
        
        sr6 = (sigma / r) ** 6
        sr12 = sr6 * sr6
        return 4.0 * epsilon * (sr12 - sr6)
    
    @ti.func
    def coulomb_potential(self, r: ti.f32, q1: ti.f32, q2: ti.f32, k: ti.f32) -> ti.f32:
        """Coulomb potential"""
        if r < 0.1:
            r = 0.1
        
        return k * q1 * q2 / r
    
    @ti.kernel
    def form_bond(self, i: ti.i32, j: ti.i32):
        """Form a bond between two particles"""
        if i >= 0 and i < self.max_particles and j >= 0 and j < self.max_particles:
            self.bond_active[i, j] = 1
            self.bond_active[j, i] = 1
            self.bond_matrix[i, j] = 1.0
            self.bond_matrix[j, i] = 1.0
            self.bond_energy[i, j] = -1.0  # Negative binding energy
            self.bond_energy[j, i] = -1.0
            self.bond_age[i, j] = 0.0
            self.bond_age[j, i] = 0.0
    
    @ti.kernel
    def break_bond(self, i: ti.i32, j: ti.i32):
        """Break a bond between two particles"""
        if i >= 0 and i < self.max_particles and j >= 0 and j < self.max_particles:
            self.bond_active[i, j] = 0
            self.bond_active[j, i] = 0
            self.bond_matrix[i, j] = 0.0
            self.bond_matrix[j, i] = 0.0
            self.bond_energy[i, j] = 0.0
            self.bond_energy[j, i] = 0.0
            self.bond_age[i, j] = 0.0
            self.bond_age[j, i] = 0.0
    
    def update_clusters(self, positions, active, particle_count: int):
        """Update cluster assignments using union-find algorithm"""
        update_clusters_kernel(active, particle_count)
        
        # SPATIAL MERGING: Connect nearby clusters for larger structures
        # DEBUG: Log spatial merging activity
        if hasattr(self, '_merge_count'):
            prev_merge_count = self._merge_count
        else:
            prev_merge_count = 0
            self._merge_count = 0
        
        # Call spatial merging and get merge count
        merge_count = self._merge_nearby_clusters(positions, active, particle_count)
        
        # DEBUG: Log if any merges occurred
        if merge_count > 0:
            print(f"SPATIAL MERGING: {merge_count} bonds created at step")
        
        # Compute cluster metrics
        compute_cluster_metrics_kernel(positions, particle_count)
    
    @ti.func
    def find_cluster_root(self, i: ti.i32) -> ti.i32:
        """Find root of cluster using iterative approach (no recursion)"""
        # Iterative path compression to avoid recursion depth issues
        root = i
        while self.cluster_id[root] != root:
            root = self.cluster_id[root]
        
        # Path compression
        current = i
        while current != root:
            next_node = self.cluster_id[current]
            self.cluster_id[current] = root
            current = next_node
        
        return root
    
    @ti.func
    def union_clusters(self, i: ti.i32, j: ti.i32):
        """Union two clusters"""
        root_i = self.find_cluster_root(i)
        root_j = self.find_cluster_root(j)
        
        if root_i != root_j:
            self.cluster_id[root_j] = root_i
    
    @ti.kernel
    def _merge_nearby_clusters(self, positions: ti.template(), active: ti.template(), particle_count: ti.i32) -> ti.i32:
        """Merge nearby clusters by creating bonds between close particles from different clusters"""
        # DEBUG: Count merges
        merge_count = 0
        
        # Find particles from different clusters that are close enough to bond
        for i in range(particle_count):
            if active[i] == 1:
                cluster_i = cluster_id_field[i]
                
                for j in range(i + 1, particle_count):
                    if active[j] == 1:
                        cluster_j = cluster_id_field[j]
                        
                        # Only consider particles from different clusters
                        if cluster_i != cluster_j:
                            pos_i = positions[i]
                            pos_j = positions[j]
                            
                            # Calculate distance
                            dx = pos_j[0] - pos_i[0]
                            dy = pos_j[1] - pos_i[1]
                            dist = ti.sqrt(dx * dx + dy * dy)
                            
                            # If particles are close enough, create a bond to merge clusters
                            merge_distance = PARTICLE_RADIUS_COMPILE * 12.0  # BALANCED: Increased from 8.0 for easier merging
                            if dist <= merge_distance and bond_active_field[i, j] == 0:
                                # Check valence constraints
                                valence_i = particle_bond_count_field[i]
                                valence_j = particle_bond_count_field[j]
                                max_valence_i = particle_valence_max_field[i]
                                max_valence_j = particle_valence_max_field[j]
                                
                                if valence_i < max_valence_i and valence_j < max_valence_j:
                                    # Create bond to merge clusters
                                    bond_active_field[i, j] = 1
                                    bond_active_field[j, i] = 1
                                    bond_matrix_field[i, j] = 1.0
                                    bond_matrix_field[j, i] = 1.0
                                    bond_type_field[i, j] = 0  # vdW bond for merging
                                    bond_type_field[j, i] = 0
                                    bond_k_spring_field[i, j] = 3.0  # BALANCED: Stronger spring for merging (was 1.0)
                                    bond_k_spring_field[j, i] = 3.0
                                    bond_rest_len_field[i, j] = dist  # Use current distance
                                    bond_rest_len_field[j, i] = dist
                                    bond_damping_field[i, j] = 0.15  # BALANCED: Higher damping (was 0.05)
                                    bond_damping_field[j, i] = 0.15
                                    bond_strength_field[i, j] = 8.0  # BALANCED: Higher strength (was 2.0)
                                    bond_strength_field[j, i] = 8.0
                                    bond_energy_field[i, j] = -2.0  # BALANCED: Stronger binding energy (was -0.5)
                                    bond_energy_field[j, i] = -2.0
                                    bond_age_field[i, j] = 0.0
                                    bond_age_field[j, i] = 0.0
                                    
                                    # Update valence counts
                                    ti.atomic_add(particle_bond_count_field[i], 1)
                                    ti.atomic_add(particle_bond_count_field[j], 1)
                                    
                                    # DEBUG: Count merges
                                    merge_count += 1
        
        return merge_count
    
    def get_bonds(self) -> List[Tuple[int, int, float]]:
        """Get list of active bonds - OPTIMIZED with numpy vectorization"""
        # PERFORMANCE FIX: Use reasonable limit for visualization (max 500 particles to match get_clusters)
        max_check = min(500, self.max_particles)
        
        # OPTIMIZATION: Use numpy for fast vectorized operations instead of Python loops
        # Copy only upper triangle to numpy (much faster than Python loops)
        bond_active_np = self.bond_active.to_numpy()[:max_check, :max_check]
        bond_matrix_np = self.bond_matrix.to_numpy()[:max_check, :max_check]
        
        # Get upper triangle indices where bonds are active
        # Use numpy to find all active bonds at once
        i_indices, j_indices = np.where(np.triu(bond_active_np, k=1) == 1)
        
        # Extract strengths for active bonds
        strengths = bond_matrix_np[i_indices, j_indices]
        
        # Convert to list of tuples (much faster than appending in loop)
        bonds = [(int(i), int(j), float(strength)) for i, j, strength in zip(i_indices, j_indices, strengths)]
        
        return bonds
    
    def get_largest_cluster(self, min_size: int = 2) -> Optional[List[int]]:
        """Get only the largest cluster - OPTIMIZED for visualization"""
        # PERFORMANCE FIX: Use reasonable limit for full cluster visualization
        max_check = min(1000, self.max_particles)
        
        # Find largest cluster ID first
        largest_cluster_id = -1
        largest_size = 0
        
        for i in range(max_check):
            cid = int(self.cluster_id[i])
            if cid >= 0:
                cluster_size = int(self.cluster_sizes[cid])
                if cluster_size >= min_size and cluster_size > largest_size:
                    largest_size = cluster_size
                    largest_cluster_id = cid
        
        if largest_cluster_id == -1:
            return None
        
        # Extract only the largest cluster
        largest_cluster = []
        for i in range(max_check):
            cid = int(self.cluster_id[i])
            if cid == largest_cluster_id:
                largest_cluster.append(i)
        
        return largest_cluster
    
    def get_largest_cluster_fast(self, min_size: int = 2) -> Optional[List[int]]:
        """Get only the largest cluster - FAST version without numpy (safer)"""
        # PERFORMANCE FIX: Use reasonable limit for full cluster visualization
        max_check = min(1000, self.max_particles)
        
        # Find largest cluster ID first
        largest_cluster_id = -1
        largest_size = 0
        
        for i in range(max_check):
            cid = int(self.cluster_id[i])
            if cid >= 0:
                cluster_size = int(self.cluster_sizes[cid])
                if cluster_size >= min_size and cluster_size > largest_size:
                    largest_size = cluster_size
                    largest_cluster_id = cid
        
        if largest_cluster_id == -1:
            return None
        
        # Extract only the largest cluster
        largest_cluster = []
        for i in range(max_check):
            cid = int(self.cluster_id[i])
            if cid == largest_cluster_id:
                largest_cluster.append(i)
        
        return largest_cluster
    
    def get_clusters(self, min_size: int = 2) -> List[List[int]]:
        """Get list of clusters with minimum size - OPTIMIZED with numpy vectorization"""
        # PERFORMANCE FIX: Use reasonable limit for visualization (max 500 particles for our current config)
        from collections import defaultdict
        
        max_check = min(500, self.max_particles)
        
        # OPTIMIZATION: Use numpy for fast vectorized operations instead of Python loops
        cluster_id_np = self.cluster_id.to_numpy()[:max_check]
        cluster_sizes_np = self.cluster_sizes.to_numpy()
        
        # Find valid clusters (size >= min_size) using numpy
        valid_mask = cluster_id_np >= 0
        particle_indices = np.arange(max_check)[valid_mask]
        cluster_ids = cluster_id_np[valid_mask]
        
        # Filter by cluster size using numpy (much faster than checking in loop)
        cluster_sizes_for_particles = cluster_sizes_np[cluster_ids]
        size_valid_mask = cluster_sizes_for_particles >= min_size
        
        # Only process particles in valid clusters
        valid_particle_indices = particle_indices[size_valid_mask]
        valid_cluster_ids = cluster_ids[size_valid_mask]
        
        # Group particles by cluster ID (still need defaultdict for grouping)
        clusters_dict = defaultdict(list)
        for particle_idx, cluster_id in zip(valid_particle_indices, valid_cluster_ids):
            clusters_dict[int(cluster_id)].append(int(particle_idx))
        
        return list(clusters_dict.values())
    
    def get_cluster_stats(self) -> Dict:
        """Get cluster statistics"""
        clusters = self.get_clusters()
        cluster_sizes = [len(cluster) for cluster in clusters]
        
        if not cluster_sizes:
            return {
                'num_clusters': 0,
                'max_cluster_size': 0,
                'average_cluster_size': 0,
                'total_clustered_particles': 0,
                'clusters': {}
            }
        
        # Convert list of lists to dict for diagnostics
        clusters_dict = {i: set(cluster) for i, cluster in enumerate(clusters)}
        
        return {
            'num_clusters': len(clusters),
            'max_cluster_size': max(cluster_sizes),
            'average_cluster_size': sum(cluster_sizes) / len(cluster_sizes),
            'total_clustered_particles': sum(cluster_sizes),
            'clusters': clusters_dict
        }
    
    def compute_cluster_metrics(self, positions_np: np.ndarray) -> Dict:
        """Compute detailed cluster metrics including R_g and graph density"""
        clusters = self.get_clusters()
        
        if not clusters:
            return {
                'R_g_values': [],
                'graph_densities': [],
                'avg_degrees': [],
                'cluster_sizes': []
            }
        
        bond_matrix = self.bond_matrix.to_numpy()
        bond_active = self.bond_active.to_numpy()
        
        R_g_values = []
        graph_densities = []
        avg_degrees = []
        cluster_sizes = []
        
        for cluster in clusters:
            size = len(cluster)
            cluster_sizes.append(size)
            
            # Compute R_g (radius of gyration)
            if size >= 2 and len(positions_np) > max(cluster):
                cluster_positions = positions_np[cluster]
                # Center of mass
                com = np.mean(cluster_positions, axis=0)
                # R_g² = mean((r_i - COM)²)
                distances_sq = np.sum((cluster_positions - com)**2, axis=1)
                R_g = np.sqrt(np.mean(distances_sq))
                R_g_values.append(R_g)
            else:
                R_g_values.append(0.0)
            
            # Compute graph density
            if size >= 2:
                # Count bonds within cluster
                bonds_in_cluster = 0
                for i in cluster:
                    for j in cluster:
                        if i < j and i < bond_active.shape[0] and j < bond_active.shape[1]:
                            if bond_active[i, j] == 1:
                                bonds_in_cluster += 1
                
                # Max possible bonds = N*(N-1)/2
                max_bonds = size * (size - 1) / 2
                density = bonds_in_cluster / max_bonds if max_bonds > 0 else 0
                graph_densities.append(density)
                
                # Average degree = 2 * edges / nodes
                avg_degree = 2 * bonds_in_cluster / size if size > 0 else 0
                avg_degrees.append(avg_degree)
            else:
                graph_densities.append(0.0)
                avg_degrees.append(0.0)
        
        return {
            'R_g_values': R_g_values,
            'graph_densities': graph_densities,
            'avg_degrees': avg_degrees,
            'cluster_sizes': cluster_sizes,
            'mean_R_g': np.mean(R_g_values) if R_g_values else 0.0,
            'mean_density': np.mean(graph_densities) if graph_densities else 0.0,
            'mean_avg_degree': np.mean(avg_degrees) if avg_degrees else 0.0
        }
    
    def get_bond_stats(self) -> Dict:
        """Get bond statistics"""
        bonds = self.get_bonds()
        
        if not bonds:
            return {
                'num_bonds': 0,
                'average_bond_strength': 0,
                'total_bond_energy': 0
            }
        
        bond_strengths = [bond[2] for bond in bonds]
        bond_energies = []
        
        bond_matrix = self.bond_matrix.to_numpy()
        bond_energy = self.bond_energy.to_numpy()
        
        for i, j, strength in bonds:
            bond_energies.append(bond_energy[i, j])
        
        return {
            'num_bonds': len(bonds),
            'average_bond_strength': sum(bond_strengths) / len(bond_strengths),
            'total_bond_energy': sum(bond_energies)
        }
