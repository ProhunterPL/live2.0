"""
Potential functions for Live 2.0 simulation
Handles particle-particle interactions and binding potentials
"""

import taichi as ti
import numpy as np
from typing import Tuple, Dict, Optional
from ..config import SimulationConfig

@ti.data_oriented
class PotentialSystem:
    """Manages potential functions and particle interactions"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        
        # Potential parameters
        self.potential_strength = ti.field(dtype=ti.f32, shape=())
        self.potential_range = ti.field(dtype=ti.f32, shape=())
        self.binding_threshold = ti.field(dtype=ti.f32, shape=())
        self.unbinding_threshold = ti.field(dtype=ti.f32, shape=())
        
        # Initialize parameters
        self.potential_strength[None] = 1.0
        self.potential_range[None] = 2.0
        self.binding_threshold[None] = config.binding_threshold
        self.unbinding_threshold[None] = config.unbinding_threshold
        
        # Force field for storing computed forces
        self.forces = ti.Vector.field(2, dtype=ti.f32, shape=(config.max_particles,))
        
        # Binding matrix for tracking particle bonds
        self.binding_matrix = ti.field(dtype=ti.f32, 
                                     shape=(config.max_particles, config.max_particles))
        
        # Initialize
        self.reset()
    
    @ti.kernel
    def reset(self):
        """Reset potential system"""
        for i in range(self.forces.shape[0]):
            self.forces[i] = ti.Vector([0.0, 0.0])
        
        for i, j in ti.ndrange(self.binding_matrix.shape[0], self.binding_matrix.shape[1]):
            self.binding_matrix[i, j] = 0.0
    
    @ti.func
    def lennard_jones_potential(self, r: ti.f32, epsilon: ti.f32 = 1.0, sigma: ti.f32 = 1.0) -> ti.f32:
        """Lennard-Jones potential: V(r) = 4ε[(σ/r)¹² - (σ/r)⁶]"""
        if r < 0.1:  # Avoid division by zero
            r = 0.1
        
        sr6 = (sigma / r) ** 6
        sr12 = sr6 * sr6
        return 4.0 * epsilon * (sr12 - sr6)
    
    @ti.func
    def lennard_jones_force(self, r: ti.f32, epsilon: ti.f32 = 1.0, sigma: ti.f32 = 1.0) -> ti.f32:
        """Force magnitude from Lennard-Jones potential"""
        if r < 0.1:
            r = 0.1
        
        sr6 = (sigma / r) ** 6
        sr12 = sr6 * sr6
        return 24.0 * epsilon * (2.0 * sr12 - sr6) / r
    
    @ti.func
    def coulomb_potential(self, r: ti.f32, q1: ti.f32, q2: ti.f32, k: ti.f32 = 1.0) -> ti.f32:
        """Coulomb potential: V(r) = k*q1*q2/r"""
        if r < 0.1:
            r = 0.1
        
        return k * q1 * q2 / r
    
    @ti.func
    def coulomb_force(self, r: ti.f32, q1: ti.f32, q2: ti.f32, k: ti.f32 = 1.0) -> ti.f32:
        """Force magnitude from Coulomb potential"""
        if r < 0.1:
            r = 0.1
        
        return k * q1 * q2 / (r * r)
    
    @ti.func
    def harmonic_potential(self, r: ti.f32, k: ti.f32, r0: ti.f32) -> ti.f32:
        """Harmonic potential: V(r) = 0.5*k*(r-r0)²"""
        dr = r - r0
        return 0.5 * k * dr * dr
    
    @ti.func
    def harmonic_force(self, r: ti.f32, k: ti.f32, r0: ti.f32) -> ti.f32:
        """Force magnitude from harmonic potential"""
        return k * (r0 - r)
    
    @ti.func
    def binding_potential(self, r: ti.f32, binding_strength: ti.f32, 
                         equilibrium_distance: ti.f32 = 1.0) -> ti.f32:
        """Binding potential for particle bonding"""
        if r < 0.1:
            r = 0.1
        
        # Attractive potential that becomes repulsive at very short distances
        dr = r - equilibrium_distance
        return binding_strength * (dr * dr - 0.5 * dr * dr * dr / equilibrium_distance)
    
    @ti.func
    def binding_force(self, r: ti.f32, binding_strength: ti.f32,
                     equilibrium_distance: ti.f32 = 1.0) -> ti.f32:
        """Force magnitude from binding potential"""
        if r < 0.1:
            r = 0.1
        
        dr = r - equilibrium_distance
        return binding_strength * (2.0 * dr - 1.5 * dr * dr / equilibrium_distance)
    
    @ti.kernel
    def compute_forces(self, positions: ti.template(), attributes: ti.template(),
                      active: ti.template(), particle_count: ti.i32):
        """Compute forces between all particle pairs"""
        # Clear forces
        for i in range(self.forces.shape[0]):
            self.forces[i] = ti.Vector([0.0, 0.0])
        
        # Compute pairwise forces
        for i in range(particle_count):
            if active[i] == 1:
                for j in range(i + 1, particle_count):
                    if active[j] == 1:
                        pos_i = positions[i]
                        pos_j = positions[j]
                        
                        # Calculate distance vector
                        r_vec = pos_i - pos_j
                        r = r_vec.norm()
                        
                        if r > 0.1:  # Avoid division by zero
                            # Normalize distance vector
                            r_hat = r_vec / r
                            
                            # Get particle attributes
                            mass_i = attributes[i][0]
                            mass_j = attributes[j][0]
                            charge_i = attributes[i][1]  # x-component of charge vector
                            charge_j = attributes[j][1]  # x-component of charge vector
                            
                            # Compute forces
                            # Lennard-Jones force
                            lj_force = self.lennard_jones_force(r, 1.0, 1.0)
                            
                            # Coulomb force
                            coulomb_force = self.coulomb_force(r, charge_i, charge_j, 1.0)
                            
                            # Total force magnitude
                            total_force = lj_force + coulomb_force
                            
                            # Apply force to both particles (Newton's third law)
                            force_vec = total_force * r_hat
                            
                            self.forces[i] += force_vec
                            self.forces[j] -= force_vec
    
    @ti.kernel
    def update_binding_matrix(self, positions: ti.template(), attributes: ti.template(),
                            active: ti.template(), particle_count: ti.i32):
        """Update binding matrix based on particle distances and properties"""
        # Clear binding matrix
        for i, j in ti.ndrange(self.binding_matrix.shape[0], self.binding_matrix.shape[1]):
            self.binding_matrix[i, j] = 0.0
        
        # Update binding strengths
        for i in range(particle_count):
            if active[i] == 1:
                for j in range(i + 1, particle_count):
                    if active[j] == 1:
                        pos_i = positions[i]
                        pos_j = positions[j]
                        
                        r_vec = pos_i - pos_j
                        r = r_vec.norm()
                        
                        if r < self.potential_range[None]:
                            # Compute binding strength based on distance and particle properties
                            binding_strength = self.compute_binding_strength(i, j, r, attributes)
                            
                            # Update binding matrix (symmetric)
                            self.binding_matrix[i, j] = binding_strength
                            self.binding_matrix[j, i] = binding_strength
    
    @ti.func
    def compute_binding_strength(self, i: ti.i32, j: ti.i32, r: ti.f32, 
                               attributes: ti.template()) -> ti.f32:
        """Compute binding strength between two particles"""
        # Get particle properties
        mass_i = attributes[i][0]
        mass_j = attributes[j][0]
        charge_i = attributes[i][1]
        charge_j = attributes[j][1]
        
        # Binding strength depends on:
        # 1. Distance (closer = stronger)
        # 2. Mass compatibility
        # 3. Charge compatibility
        
        distance_factor = ti.exp(-r / self.potential_range[None])
        mass_factor = ti.min(mass_i, mass_j) / ti.max(mass_i, mass_j)
        charge_factor = 1.0 - ti.abs(charge_i - charge_j) / ti.max(ti.abs(charge_i), ti.abs(charge_j), 1.0)
        
        binding_strength = distance_factor * mass_factor * charge_factor * self.potential_strength[None]
        
        return binding_strength
    
    @ti.kernel
    def apply_binding_forces(self, positions: ti.template(), attributes: ti.template(),
                           active: ti.template(), particle_count: ti.i32):
        """Apply binding forces based on current binding matrix"""
        for i in range(particle_count):
            if active[i] == 1:
                for j in range(i + 1, particle_count):
                    if active[j] == 1:
                        binding_strength = self.binding_matrix[i, j]
                        
                        if binding_strength > self.binding_threshold[None]:
                            pos_i = positions[i]
                            pos_j = positions[j]
                            
                            r_vec = pos_i - pos_j
                            r = r_vec.norm()
                            
                            if r > 0.1:
                                r_hat = r_vec / r
                                
                                # Binding force (attractive)
                                binding_force = self.binding_force(r, binding_strength, 1.0)
                                force_vec = binding_force * r_hat
                                
                                self.forces[i] += force_vec
                                self.forces[j] -= force_vec
    
    def get_forces(self) -> np.ndarray:
        """Get computed forces as numpy array"""
        return self.forces.to_numpy()
    
    def get_binding_matrix(self) -> np.ndarray:
        """Get binding matrix as numpy array"""
        return self.binding_matrix.to_numpy()
    
    def get_bonds(self, threshold: float = None) -> list:
        """Get list of particle bonds above threshold"""
        if threshold is None:
            threshold = self.binding_threshold[None]
        
        binding_matrix = self.binding_matrix.to_numpy()
        bonds = []
        
        for i in range(binding_matrix.shape[0]):
            for j in range(i + 1, binding_matrix.shape[1]):
                if binding_matrix[i, j] > threshold:
                    bonds.append((i, j, binding_matrix[i, j]))
        
        return bonds
    
    def set_potential_parameters(self, strength: float = None, range_val: float = None,
                                binding_threshold: float = None, unbinding_threshold: float = None):
        """Set potential parameters"""
        if strength is not None:
            self.potential_strength[None] = strength
        if range_val is not None:
            self.potential_range[None] = range_val
        if binding_threshold is not None:
            self.binding_threshold[None] = binding_threshold
        if unbinding_threshold is not None:
            self.unbinding_threshold[None] = unbinding_threshold
    
    def get_stats(self) -> Dict:
        """Get potential system statistics"""
        binding_matrix = self.binding_matrix.to_numpy()
        bonds = self.get_bonds()
        
        return {
            'potential_strength': self.potential_strength[None],
            'potential_range': self.potential_range[None],
            'binding_threshold': self.binding_threshold[None],
            'unbinding_threshold': self.unbinding_threshold[None],
            'total_bonds': len(bonds),
            'max_binding_strength': float(np.max(binding_matrix)),
            'average_binding_strength': float(np.mean(binding_matrix[binding_matrix > 0]))
        }
