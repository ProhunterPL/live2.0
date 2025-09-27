"""
Binding system for Live 2.0 simulation
Handles particle binding, bond formation/breaking, and cluster detection
"""

import taichi as ti
import numpy as np
from typing import List, Tuple, Dict, Set
from .config import SimulationConfig

@ti.data_oriented
class BindingSystem:
    """Manages particle binding and bond formation"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.max_particles = config.max_particles
        
        # Bond tracking
        self.bond_matrix = ti.field(dtype=ti.f32, 
                                  shape=(self.max_particles, self.max_particles))
        self.bond_active = ti.field(dtype=ti.i32,
                                  shape=(self.max_particles, self.max_particles))
        
        # Cluster tracking
        self.cluster_id = ti.field(dtype=ti.i32, shape=(self.max_particles,))
        self.cluster_sizes = ti.field(dtype=ti.i32, shape=(self.max_particles,))
        self.next_cluster_id = ti.field(dtype=ti.i32, shape=())
        
        # Bond properties
        self.bond_energy = ti.field(dtype=ti.f32, 
                                   shape=(self.max_particles, self.max_particles))
        self.bond_age = ti.field(dtype=ti.f32,
                               shape=(self.max_particles, self.max_particles))
        
        # Initialize
        self.reset()
    
    @ti.kernel
    def reset(self):
        """Reset binding system"""
        for i, j in ti.ndrange(self.max_particles, self.max_particles):
            self.bond_matrix[i, j] = 0.0
            self.bond_active[i, j] = 0
            self.bond_energy[i, j] = 0.0
            self.bond_age[i, j] = 0.0
        
        for i in range(self.max_particles):
            self.cluster_id[i] = -1
            self.cluster_sizes[i] = 0
        
        self.next_cluster_id[None] = 0
    
    @ti.kernel
    def update_bonds(self, positions: ti.template(), attributes: ti.template(),
                    active: ti.template(), particle_count: ti.i32, dt: ti.f32):
        """Update bond formation and breaking based on particle interactions"""
        # Update bond ages
        for i, j in ti.ndrange(self.max_particles, self.max_particles):
            if self.bond_active[i, j] == 1:
                self.bond_age[i, j] += dt
        
        # Check for new bonds
        for i in range(particle_count):
            if active[i] == 1:
                for j in range(i + 1, particle_count):
                    if active[j] == 1 and self.bond_active[i, j] == 0:
                        # Check if particles should form a bond
                        if self.should_form_bond(i, j, positions, attributes):
                            self.form_bond(i, j)
        
        # Check for bond breaking
        for i, j in ti.ndrange(self.max_particles, self.max_particles):
            if self.bond_active[i, j] == 1:
                if self.should_break_bond(i, j, positions, attributes):
                    self.break_bond(i, j)
    
    @ti.func
    def should_form_bond(self, i: ti.i32, j: ti.i32, positions: ti.template(),
                        attributes: ti.template()) -> ti.i32:
        """Check if two particles should form a bond"""
        pos_i = positions[i]
        pos_j = positions[j]
        
        # Calculate distance
        r_vec = pos_i - pos_j
        r = r_vec.norm()
        
        # Check distance threshold
        if r > self.config.particle_radius * 2.5:  # Within binding range
            return 0
        
        # Check binding compatibility
        mass_i = attributes[i][0]
        mass_j = attributes[j][0]
        charge_i = attributes[i][1]
        charge_j = attributes[j][1]
        
        # Mass compatibility
        mass_ratio = ti.min(mass_i, mass_j) / ti.max(mass_i, mass_j)
        if mass_ratio < 0.5:  # Too different masses
            return 0
        
        # Charge compatibility (opposite charges attract)
        charge_product = charge_i * charge_j
        if charge_product > 0.1:  # Same charges repel
            return 0
        
        # Energy threshold
        binding_energy = self.compute_binding_energy(i, j, r, attributes)
        if binding_energy < -0.5:  # Sufficiently negative binding energy
            return 1
        
        return 0
    
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
        if r > self.config.particle_radius * 4.0:
            return 1
        
        # Break if bond is too weak
        bond_strength = self.bond_matrix[i, j]
        if bond_strength < self.config.unbinding_threshold:
            return 1
        
        # Break if particles have incompatible properties
        mass_i = attributes[i][0]
        mass_j = attributes[j][0]
        charge_i = attributes[i][1]
        charge_j = attributes[j][1]
        
        mass_ratio = ti.min(mass_i, mass_j) / ti.max(mass_i, mass_j)
        if mass_ratio < 0.3:  # Masses became too different
            return 1
        
        charge_product = charge_i * charge_j
        if charge_product > 0.2:  # Charges became incompatible
            return 1
        
        return 0
    
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
    
    @ti.kernel
    def update_clusters(self, active: ti.template(), particle_count: ti.i32):
        """Update cluster assignments using union-find algorithm"""
        # Initialize cluster IDs
        for i in range(particle_count):
            if active[i] == 1:
                self.cluster_id[i] = i
            else:
                self.cluster_id[i] = -1
        
        # Union connected particles
        for i in range(particle_count):
            if active[i] == 1:
                for j in range(i + 1, particle_count):
                    if active[j] == 1 and self.bond_active[i, j] == 1:
                        self.union_clusters(i, j)
        
        # Calculate cluster sizes
        for i in range(self.max_particles):
            self.cluster_sizes[i] = 0
        
        for i in range(particle_count):
            if active[i] == 1:
                root = self.find_cluster_root(i)
                self.cluster_sizes[root] += 1
    
    @ti.func
    def find_cluster_root(self, i: ti.i32) -> ti.i32:
        """Find root of cluster using path compression"""
        if self.cluster_id[i] != i:
            self.cluster_id[i] = self.find_cluster_root(self.cluster_id[i])
        return self.cluster_id[i]
    
    @ti.func
    def union_clusters(self, i: ti.i32, j: ti.i32):
        """Union two clusters"""
        root_i = self.find_cluster_root(i)
        root_j = self.find_cluster_root(j)
        
        if root_i != root_j:
            self.cluster_id[root_j] = root_i
    
    def get_bonds(self) -> List[Tuple[int, int, float]]:
        """Get list of active bonds"""
        bond_matrix = self.bond_matrix.to_numpy()
        bond_active = self.bond_active.to_numpy()
        
        bonds = []
        for i in range(bond_matrix.shape[0]):
            for j in range(i + 1, bond_matrix.shape[1]):
                if bond_active[i, j] == 1:
                    bonds.append((i, j, bond_matrix[i, j]))
        
        return bonds
    
    def get_clusters(self, min_size: int = 2) -> List[List[int]]:
        """Get list of clusters with minimum size"""
        cluster_id = self.cluster_id.to_numpy()
        cluster_sizes = self.cluster_sizes.to_numpy()
        
        # Group particles by cluster
        clusters = {}
        for i, cid in enumerate(cluster_id):
            if cid >= 0 and cluster_sizes[cid] >= min_size:
                if cid not in clusters:
                    clusters[cid] = []
                clusters[cid].append(i)
        
        return list(clusters.values())
    
    def get_cluster_stats(self) -> Dict:
        """Get cluster statistics"""
        clusters = self.get_clusters()
        cluster_sizes = [len(cluster) for cluster in clusters]
        
        if not cluster_sizes:
            return {
                'num_clusters': 0,
                'max_cluster_size': 0,
                'average_cluster_size': 0,
                'total_clustered_particles': 0
            }
        
        return {
            'num_clusters': len(clusters),
            'max_cluster_size': max(cluster_sizes),
            'average_cluster_size': sum(cluster_sizes) / len(cluster_sizes),
            'total_clustered_particles': sum(cluster_sizes)
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
