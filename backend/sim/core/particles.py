"""
Particle management for Live 2.0 simulation
Handles particle properties, attributes, and basic operations
"""

import taichi as ti
import numpy as np
from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)
from ..config import SimulationConfig

# Global variables for Taichi fields (will be initialized later)
positions = None
velocities = None
attributes = None
active = None
particle_count = None
type_ids = None
binding_sites = None
binding_strength = None
energy = None
age = None
last_mutation = None

# Compile-time constant for max particles
MAX_PARTICLES_COMPILE = 10000

def init_particle_fields():
    """Initialize Taichi fields after ti.init() has been called"""
    global positions, velocities, attributes, active, particle_count
    global type_ids, binding_sites, binding_strength
    global energy, age, last_mutation
    
    positions = ti.Vector.field(2, dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    velocities = ti.Vector.field(2, dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    attributes = ti.Vector.field(4, dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    active = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    particle_count = ti.field(dtype=ti.i32, shape=())
    
    type_ids = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    binding_sites = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    binding_strength = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    
    # Dynamic properties
    energy = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    age = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    last_mutation = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))

# Define kernels at module level with compile-time constants
@ti.kernel
def reset_particles_kernel():
    """Reset all particles - module-level kernel"""
    for i in range(MAX_PARTICLES_COMPILE):
        active[i] = 0
        energy[i] = 0.0
        age[i] = 0.0
        last_mutation[i] = 0.0
    particle_count[None] = 0

@ti.kernel
def update_positions_kernel(dt: ti.f32):
    """Update particle positions - module-level kernel"""
    for i in range(MAX_PARTICLES_COMPILE):
        if active[i] == 1:
            positions[i] += velocities[i] * dt
            age[i] += dt

@ti.kernel
def apply_forces_kernel(forces: ti.template(), dt: ti.f32):
    """Apply forces to particles - module-level kernel"""
    for i in range(MAX_PARTICLES_COMPILE):
        if active[i] == 1:
            mass = attributes[i][0]
            if mass > 0:
                acceleration = forces[i] / mass
                velocities[i] += acceleration * dt

@ti.kernel
def remove_particle_kernel(idx: ti.i32):
    """Remove a particle - module-level kernel"""
    if idx >= 0 and idx < MAX_PARTICLES_COMPILE:
        active[idx] = 0

@ti.kernel
def add_energy_kernel(energy_amount: ti.template()):
    """Add energy to particles - module-level kernel"""
    for i in range(MAX_PARTICLES_COMPILE):
        if active[i] == 1:
            energy[i] += energy_amount[i]

@ti.kernel
def decay_energy_kernel(decay_rate: ti.f32):
    """Decay particle energy - module-level kernel"""
    for i in range(MAX_PARTICLES_COMPILE):
        if active[i] == 1:
            energy[i] *= decay_rate

@ti.kernel
def mutate_particle_kernel(idx: ti.i32, mutation_strength: ti.f32, 
                           current_time: ti.f32, rng: ti.template()):
    """Apply mutation to particle - module-level kernel"""
    if idx >= 0 and idx < MAX_PARTICLES_COMPILE and active[idx] == 1:
        # Check if enough time has passed since last mutation
        if current_time - last_mutation[idx] >= 1.0:
            # Mutate attributes
            for j in range(4):
                mutation = rng.next_gaussian(0.0, mutation_strength)
                attributes[idx][j] += mutation
                
                # Keep mass positive
                if j == 0 and attributes[idx][j] <= 0:
                    attributes[idx][j] = 0.1
            
            # Update mutation time
            last_mutation[idx] = current_time

@ti.kernel
def get_particle_neighbors_kernel(idx: ti.i32, positions_field: ti.template(),
                                 radius: ti.f32, neighbors: ti.template()) -> ti.i32:
    """Get neighbors of a particle - module-level kernel"""
    count = 0
    if idx >= 0 and idx < MAX_PARTICLES_COMPILE and active[idx] == 1:
        pos = positions[idx]
        for i in range(MAX_PARTICLES_COMPILE):
            if active[i] == 1 and i != idx:
                diff = positions[i] - pos
                dist_sq = diff[0] * diff[0] + diff[1] * diff[1]
                if dist_sq < radius * radius:
                    if count < neighbors.shape[0]:
                        neighbors[count] = i
                        count += 1
    return count

@ti.data_oriented
class ParticleSystem:
    """Manages particles and their properties in the simulation"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.max_particles = config.max_particles
        
        # Initialize Taichi fields if not already done
        if positions is None:
            init_particle_fields()
        
        # Use module-level fields
        self.positions = positions
        self.velocities = velocities
        self.attributes = attributes
        self.active = active
        self.particle_count = particle_count
        
        # Particle types and properties
        self.type_ids = type_ids
        self.binding_sites = binding_sites
        self.binding_strength = binding_strength
        
        # Dynamic properties
        self.energy = energy
        self.age = age
        self.last_mutation = last_mutation
        
        # Pre-allocated Taichi fields for efficient data transfer
        self._numpy_positions_taichi = ti.field(dtype=ti.f32, shape=(self.max_particles, 2))
        self._numpy_velocities_taichi = ti.field(dtype=ti.f32, shape=(self.max_particles, 2))
        self._numpy_attributes_taichi = ti.field(dtype=ti.f32, shape=(self.max_particles, 4))
        self._numpy_active_taichi = ti.field(dtype=ti.i32, shape=(self.max_particles))
        
        # Type registry for tracking particle types
        self.type_registry = {}
        self.next_type_id = 0
        
        # Initialize
        self.reset()
    
    def reset(self):
        """Reset particle system to initial state"""
        reset_particles_kernel()
    
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
    
    def add_particle_py(self, pos, vel, attributes, type_id: int,
                        binding_sites: int, binding_strength: float) -> int:
        """Add a particle from Python scope to avoid Taichi kernel return constraints."""
        idx = int(self.particle_count[None])
        print(f"DEBUG: add_particle_py called, current count={idx}, max={self.max_particles}")
        if idx >= int(self.max_particles):
            print(f"DEBUG: add_particle_py failed - too many particles")
            return -1
        self.positions[idx] = pos
        self.velocities[idx] = vel
        self.attributes[idx] = attributes
        self.type_ids[idx] = int(type_id)
        self.binding_sites[idx] = int(binding_sites)
        self.binding_strength[idx] = float(binding_strength)
        self.active[idx] = 1
        self.energy[idx] = 0.0
        self.age[idx] = 0.0
        self.last_mutation[idx] = 0.0
        self.particle_count[None] = idx + 1
        print(f"DEBUG: add_particle_py success, new count={self.particle_count[None]}")
        return idx
    
    def remove_particle(self, idx: int):
        """Remove particle by index"""
        remove_particle_kernel(idx)
    
    def update_positions(self, dt: float):
        """Update particle positions using velocity"""
        update_positions_kernel(dt)
    
    def apply_forces(self, forces, dt: float):
        """Apply forces to particles and update velocities"""
        apply_forces_kernel(forces, dt)
    
    def add_energy(self, energy_amount):
        """Add energy to particles"""
        # Debug energy addition
        total_energy_added = 0.0
        for i in range(min(10, self.particle_count[None])):  # Check first 10 particles
            if self.active[i] == 1:
                total_energy_added += energy_amount[i]
        logger.info(f"DEBUG: add_energy called, total energy to add={total_energy_added:.4f}")
        
        add_energy_kernel(energy_amount)
    
    def decay_energy(self, decay_rate: float):
        """Apply energy decay to all particles"""
        decay_energy_kernel(decay_rate)
    
    def mutate_particle(self, idx: int, mutation_strength: float, 
                       current_time: float, rng):
        """Apply mutation to particle attributes"""
        # Call module-level kernel
        mutate_particle_kernel(idx, mutation_strength, current_time, rng)
    
    def get_particle_neighbors(self, idx: int, positions,
                             radius: float, neighbors) -> int:
        """Get neighbors of a particle within radius"""
        return get_particle_neighbors_kernel(idx, positions, radius, neighbors)
    
    def get_active_particles(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get data for all active particles - OPTIMIZED VERSION with energy"""
        # Use single kernel to copy data to Taichi fields more efficiently
        self._copy_particles_to_taichi()
        
        # Convert to numpy and filter only active particles
        positions = self._numpy_positions_taichi.to_numpy()
        velocities = self._numpy_velocities_taichi.to_numpy()
        attributes = self._numpy_attributes_taichi.to_numpy()
        active = self._numpy_active_taichi.to_numpy()
        energies = self.energy.to_numpy()
        
        # Filter only active particles
        active_mask = active == 1
        active_positions = positions[active_mask]
        active_velocities = velocities[active_mask]
        active_attributes = attributes[active_mask]
        active_energies = energies[active_mask]
        
        return active_positions, active_velocities, active_attributes, active_mask, active_energies
    
    @ti.kernel
    def _copy_particles_to_taichi(self):
        """Copy particle data to Taichi fields using single kernel"""
        for i in range(self.max_particles):
            # Copy vector fields to 2D fields
            self._numpy_positions_taichi[i, 0] = self.positions[i][0]
            self._numpy_positions_taichi[i, 1] = self.positions[i][1]
            self._numpy_velocities_taichi[i, 0] = self.velocities[i][0]
            self._numpy_velocities_taichi[i, 1] = self.velocities[i][1]
            self._numpy_attributes_taichi[i, 0] = self.attributes[i][0]
            self._numpy_attributes_taichi[i, 1] = self.attributes[i][1]
            self._numpy_attributes_taichi[i, 2] = self.attributes[i][2]
            self._numpy_attributes_taichi[i, 3] = self.attributes[i][3]
            self._numpy_active_taichi[i] = self.active[i]
    
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
    
    @ti.kernel
    def thermal_kick(self, vmax: float, k_sigma: float, energy_field: ti.template()):
        """Add thermal noise to particle velocities based on local energy"""
        for p in range(self.max_particles):
            if self.active[p] == 1:
                # Get grid position
                i = ti.cast(self.positions[p].x, ti.i32) % energy_field.shape[0]
                j = ti.cast(self.positions[p].y, ti.i32) % energy_field.shape[1]
                E = energy_field[i, j]
                
                # Sigma grows with energy (clipped)
                sigma = k_sigma * ti.min(1.0, E)
                
                # Add thermal noise
                self.velocities[p].x += (ti.random(ti.f32) - 0.5) * sigma
                self.velocities[p].y += (ti.random(ti.f32) - 0.5) * sigma
                
                # Clamp velocity
                v2 = self.velocities[p].x * self.velocities[p].x + self.velocities[p].y * self.velocities[p].y
                if v2 > vmax * vmax:
                    s = vmax / ti.sqrt(v2)
                    self.velocities[p].x *= s
                    self.velocities[p].y *= s