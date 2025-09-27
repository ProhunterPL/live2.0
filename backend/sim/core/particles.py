"""
Particle management for Live 2.0 simulation
Handles particle properties, attributes, and basic operations
"""

import taichi as ti
import numpy as np
from typing import Dict, List, Tuple, Optional
from .config import SimulationConfig

@ti.data_oriented
class ParticleSystem:
    """Manages particles and their properties in the simulation"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.max_particles = config.max_particles
        
        # Particle data structures
        self.positions = ti.Vector.field(2, dtype=ti.f32, shape=(self.max_particles,))
        self.velocities = ti.Vector.field(2, dtype=ti.f32, shape=(self.max_particles,))
        self.attributes = ti.Vector.field(4, dtype=ti.f32, shape=(self.max_particles,))  # [mass, charge_x, charge_y, charge_z]
        self.active = ti.field(dtype=ti.i32, shape=(self.max_particles,))
        self.particle_count = ti.field(dtype=ti.i32, shape=())
        
        # Particle types and properties
        self.type_ids = ti.field(dtype=ti.i32, shape=(self.max_particles,))
        self.binding_sites = ti.field(dtype=ti.i32, shape=(self.max_particles,))
        self.binding_strength = ti.field(dtype=ti.f32, shape=(self.max_particles,))
        
        # Dynamic properties
        self.energy = ti.field(dtype=ti.f32, shape=(self.max_particles,))
        self.age = ti.field(dtype=ti.f32, shape=(self.max_particles,))
        self.last_mutation = ti.field(dtype=ti.f32, shape=(self.max_particles,))
        
        # Type registry for tracking particle types
        self.type_registry = {}
        self.next_type_id = 0
        
        # Initialize
        self.reset()
    
    @ti.kernel
    def reset(self):
        """Reset particle system to initial state"""
        for i in range(self.max_particles):
            self.active[i] = 0
            self.energy[i] = 0.0
            self.age[i] = 0.0
            self.last_mutation[i] = 0.0
        
        self.particle_count[None] = 0
    
    def register_particle_type(self, name: str, mass: float = 1.0, 
                            charge: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                            binding_sites: int = 2, binding_strength: float = 1.0) -> int:
        """Register a new particle type. Returns type ID."""
        type_id = self.next_type_id
        self.next_type_id += 1
        
        self.type_registry[name] = {
            'id': type_id,
            'mass': mass,
            'charge': charge,
            'binding_sites': binding_sites,
            'binding_strength': binding_strength
        }
        
        return type_id
    
    @ti.kernel
    def add_particle(self, pos: ti.template(), vel: ti.template(), 
                    attributes: ti.template(), type_id: ti.i32,
                    binding_sites: ti.i32, binding_strength: ti.f32) -> ti.i32:
        """Add a particle to the system. Returns particle index or -1 if failed"""
        if self.particle_count[None] >= self.max_particles:
            return -1
        
        idx = self.particle_count[None]
        self.positions[idx] = pos
        self.velocities[idx] = vel
        self.attributes[idx] = attributes
        self.type_ids[idx] = type_id
        self.binding_sites[idx] = binding_sites
        self.binding_strength[idx] = binding_strength
        self.active[idx] = 1
        self.energy[idx] = 0.0
        self.age[idx] = 0.0
        self.last_mutation[idx] = 0.0
        
        self.particle_count[None] += 1
        return idx
    
    @ti.kernel
    def remove_particle(self, idx: ti.i32):
        """Remove particle by index"""
        if idx < 0 or idx >= self.max_particles:
            return
        
        self.active[idx] = 0
    
    @ti.kernel
    def update_positions(self, dt: ti.f32):
        """Update particle positions using velocity"""
        for i in range(self.max_particles):
            if self.active[i] == 1:
                self.positions[i] += self.velocities[i] * dt
                self.age[i] += dt
    
    @ti.kernel
    def apply_forces(self, forces: ti.template(), dt: ti.f32):
        """Apply forces to particles and update velocities"""
        for i in range(self.max_particles):
            if self.active[i] == 1:
                mass = self.attributes[i][0]
                if mass > 0:
                    acceleration = forces[i] / mass
                    self.velocities[i] += acceleration * dt
    
    @ti.kernel
    def add_energy(self, energy_amount: ti.template()):
        """Add energy to particles"""
        for i in range(self.max_particles):
            if self.active[i] == 1:
                self.energy[i] += energy_amount[i]
    
    @ti.kernel
    def decay_energy(self, decay_rate: ti.f32):
        """Apply energy decay to all particles"""
        for i in range(self.max_particles):
            if self.active[i] == 1:
                self.energy[i] *= decay_rate
    
    @ti.kernel
    def mutate_particle(self, idx: ti.i32, mutation_strength: ti.f32, 
                       current_time: ti.f32, rng: ti.template()):
        """Apply mutation to particle attributes"""
        if idx < 0 or idx >= self.max_particles or self.active[idx] == 0:
            return
        
        # Check if enough time has passed since last mutation
        if current_time - self.last_mutation[idx] < 1.0:  # Minimum 1 time unit between mutations
            return
        
        # Mutate attributes
        for j in range(4):
            mutation = rng.next_gaussian(0.0, mutation_strength)
            self.attributes[idx][j] += mutation
            
            # Keep mass positive
            if j == 0 and self.attributes[idx][j] <= 0:
                self.attributes[idx][j] = 0.1
        
        # Update mutation time
        self.last_mutation[idx] = current_time
    
    @ti.kernel
    def get_particle_neighbors(self, idx: ti.i32, positions: ti.template(),
                             radius: ti.f32, neighbors: ti.template()) -> ti.i32:
        """Get neighbors of a particle within radius"""
        if idx < 0 or idx >= self.max_particles or self.active[idx] == 0:
            return 0
        
        pos = self.positions[idx]
        count = 0
        
        for i in range(self.max_particles):
            if i != idx and self.active[i] == 1:
                other_pos = self.positions[i]
                distance = (pos - other_pos).norm()
                
                if distance <= radius and count < neighbors.shape[0]:
                    neighbors[count] = i
                    count += 1
        
        return count
    
    def get_active_particles(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get data for all active particles"""
        positions = self.positions.to_numpy()
        velocities = self.velocities.to_numpy()
        attributes = self.attributes.to_numpy()
        active = self.active.to_numpy()
        
        # Filter only active particles
        active_mask = active == 1
        active_positions = positions[active_mask]
        active_velocities = velocities[active_mask]
        active_attributes = attributes[active_mask]
        
        return active_positions, active_velocities, active_attributes, active_mask
    
    def get_particle_by_index(self, idx: int) -> Optional[Dict]:
        """Get particle data by index"""
        if idx < 0 or idx >= self.max_particles or self.active[idx] == 0:
            return None
        
        return {
            'position': self.positions[idx].to_numpy(),
            'velocity': self.velocities[idx].to_numpy(),
            'attributes': self.attributes[idx].to_numpy(),
            'type_id': self.type_ids[idx],
            'binding_sites': self.binding_sites[idx],
            'binding_strength': self.binding_strength[idx],
            'energy': self.energy[idx],
            'age': self.age[idx]
        }
    
    def get_stats(self) -> Dict:
        """Get particle system statistics"""
        active_count = 0
        total_energy = 0.0
        total_mass = 0.0
        
        for i in range(self.max_particles):
            if self.active[i] == 1:
                active_count += 1
                total_energy += self.energy[i]
                total_mass += self.attributes[i][0]
        
        return {
            'particle_count': active_count,
            'max_particles': self.max_particles,
            'total_energy': total_energy,
            'total_mass': total_mass,
            'registered_types': len(self.type_registry),
            'average_energy': total_energy / max(active_count, 1),
            'average_mass': total_mass / max(active_count, 1)
        }
    
    def get_type_info(self, type_id: int) -> Optional[Dict]:
        """Get information about a particle type"""
        for name, info in self.type_registry.items():
            if info['id'] == type_id:
                return info
        return None
