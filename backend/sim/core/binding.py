"""
Binding system for Live 2.0 simulation
Handles particle binding, bond formation/breaking, and cluster detection
"""

import taichi as ti
import numpy as np
from typing import List, Tuple, Dict, Set
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
        particle_valence_max_field[i] = 4  # Default valence limit
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
                        if bond_type == 0:  # van der Waals
                            k_spring = 2.0
                            rest_len = dist  # Use current distance
                            damping = 0.1
                            strength = 5.0
                        elif bond_type == 1:  # covalent
                            k_spring = 10.0
                            rest_len = dist * 0.9  # Slightly shorter
                            damping = 0.2
                            strength = 20.0
                        elif bond_type == 2:  # H-bond
                            k_spring = 5.0
                            rest_len = dist * 1.1  # Slightly longer
                            damping = 0.15
                            strength = 10.0
                        else:  # metallic or other
                            k_spring = 7.0
                            rest_len = dist
                            damping = 0.25
                            strength = 15.0
                        
                        # Check valence constraints BEFORE forming bond
                        valence_i = particle_bond_count_field[i]
                        valence_j = particle_bond_count_field[j]
                        max_valence_i = particle_valence_max_field[i]
                        max_valence_j = particle_valence_max_field[j]
                        
                        # Only form bond if both particles have valence available
                        if valence_i < max_valence_i and valence_j < max_valence_j:
                            # Probabilistic bond formation - MUCH LOWER THRESHOLD FOR MORE BONDS
                            theta_bind = 0.05  # Even lower for much easier bonding
                            dE = -dist  # Simplified energy change
                            
                            # Sigmoid probability with easier formation
                            p = 1.0 / (1.0 + ti.exp(-(-dE - theta_bind) * 2.0))
                            
                            # Energy-dependent noise - MORE AGGRESSIVE BONDING
                            # Get local energy (simplified)
                            E = 0.5  # Placeholder - should get from energy field
                            p *= (1.0 + 1.0 * E)  # Increased energy factor
                            
                            # Additional distance factor - closer particles bond easier
                            distance_factor = 1.0 / (1.0 + dist * 0.5)
                            p *= distance_factor
                            
                            if ti.random(ti.f32) < p:
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
                                bond_energy_field[i, j] = 0.0
                                bond_energy_field[j, i] = 0.0
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
    
    if r <= PARTICLE_RADIUS_COMPILE * 4.0:  # INCREASED binding range (was 2.5, now 4.0)
        # Get particle properties
        mass_i = attributes[i][0]
        mass_j = attributes[j][0]
        charge_i = attributes[i][1]
        charge_j = attributes[j][1]
        
        # Calculate compatibility metrics
        mass_ratio = ti.min(mass_i, mass_j) / ti.max(mass_i, mass_j)
        charge_product = charge_i * charge_j
        charge_sum = ti.abs(charge_i + charge_j)
        
        # Determine bond type based on properties - EXTREMELY RELAXED for more bonds
        if mass_ratio > 0.2 and ti.abs(charge_product) > -1.0:  # Very relaxed
            # Similar mass, moderate charge interaction → covalent-like (strong)
            bond_type = 1  # covalent
        elif charge_product < 0.0:  # Any opposite charge
            # Opposite charges → hydrogen bond-like
            bond_type = 2  # H-bond
        elif mass_ratio > 0.1 and charge_sum > 0.0:  # Very relaxed metal
            # Medium mass ratio, moderate total charge → metallic-like
            bond_type = 3  # metallic
        else:
            # ALWAYS form vdW bond if within range
            bond_type = 0  # vdW
    
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
    strain = ti.abs(r - rest_len) / ti.max(rest_len, 0.1)
    if strain > 2.0:  # 200% strain threshold (much more stable bonds)
        result = 1
    
    # Condition 2: Distance too large (safety check)
    if r > PARTICLE_RADIUS_COMPILE * 5.0:
        result = 1
    
    # Condition 3: Aging - probabilistic breaking for old bonds
    # Different bond types have different max ages (MUCH INCREASED for stability)
    max_age = 2000.0  # default (much increased from 500.0)
    if bond_type == 0:  # vdW - short lived
        max_age = 1000.0  # increased from 200.0
    elif bond_type == 1:  # covalent - long lived
        max_age = 5000.0  # increased from 1000.0
    elif bond_type == 2:  # H-bond - medium
        max_age = 2000.0  # increased from 400.0
    else:  # metallic
        max_age = 3000.0  # increased from 600.0
    
    if age > max_age:
        # Probabilistic: higher age = higher chance to break
        # Simple deterministic approximation: break if age > 1.5*max_age
        if age > max_age * 1.5:
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
        
        # Bond type parameters (defaults)
        # 0 = van der Waals (weak)
        # 1 = covalent (strong)
        # 2 = hydrogen bond (medium)
        # 3 = metallic (medium-strong)
        self.bond_type_params = {
            0: {'k_spring': 2.0, 'rest_len': 1.0, 'damping': 0.1, 'strength': 5.0},   # vdW
            1: {'k_spring': 10.0, 'rest_len': 0.8, 'damping': 0.2, 'strength': 20.0}, # covalent
            2: {'k_spring': 5.0, 'rest_len': 1.2, 'damping': 0.15, 'strength': 10.0}, # H-bond
            3: {'k_spring': 7.0, 'rest_len': 0.9, 'damping': 0.25, 'strength': 15.0}  # metallic
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
    
    def get_bonds(self) -> List[Tuple[int, int, float]]:
        """Get list of active bonds - OPTIMIZED with limited matrix size"""
        # OPTIMIZATION: Limit to reasonable number of particles (1000 max) before to_numpy()
        max_check = min(1000, self.max_particles)
        
        # Convert only limited portion to numpy
        bond_matrix = self.bond_matrix.to_numpy()[:max_check, :max_check]
        bond_active = self.bond_active.to_numpy()[:max_check, :max_check]
        
        # Use NumPy to find active bonds (fast)
        import numpy as np
        i_indices, j_indices = np.where(np.triu(bond_active, k=1) == 1)
        
        # Build bond list using vectorized operations
        bonds = [(int(i), int(j), float(bond_matrix[i, j])) 
                 for i, j in zip(i_indices, j_indices)]
        
        return bonds
    
    def get_clusters(self, min_size: int = 2) -> List[List[int]]:
        """Get list of clusters with minimum size - OPTIMIZED with limited array size"""
        # OPTIMIZATION: Limit to reasonable number of particles (1000 max) before to_numpy()
        max_check = min(1000, self.max_particles)
        
        # Convert only limited portion to numpy
        cluster_id = self.cluster_id.to_numpy()[:max_check]
        cluster_sizes = self.cluster_sizes.to_numpy()[:max_check]
        
        # Use NumPy boolean indexing for faster filtering
        import numpy as np
        
        # Find valid clusters (size >= min_size)
        valid_cluster_mask = (cluster_id >= 0) & (cluster_sizes[cluster_id] >= min_size)
        valid_particles = np.where(valid_cluster_mask)[0]
        valid_cluster_ids = cluster_id[valid_particles]
        
        # Group particles by cluster using defaultdict
        from collections import defaultdict
        clusters_dict = defaultdict(list)
        for particle_idx, cid in zip(valid_particles, valid_cluster_ids):
            clusters_dict[int(cid)].append(int(particle_idx))
        
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
