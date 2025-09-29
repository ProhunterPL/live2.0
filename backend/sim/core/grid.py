"""
Grid management for Live 2.0 simulation
Handles 2D periodic grid with particle positioning and neighbor finding
"""

import taichi as ti
import numpy as np
from typing import Tuple, List
from ..config import SimulationConfig

# Compile-time constants
MAX_PARTICLES_COMPILE = 10000
GRID_WIDTH_COMPILE = 256
GRID_HEIGHT_COMPILE = 256
GRID_CELLS_X_COMPILE = 128  # GRID_WIDTH_COMPILE / 2.0
GRID_CELLS_Y_COMPILE = 128  # GRID_HEIGHT_COMPILE / 2.0
MAX_PARTICLES_PER_CELL_COMPILE = 32

# Global Taichi fields
particle_positions_field = None
particle_attributes_field = None
particle_active_field = None
particle_count_field = None
spatial_hash_field = None
cell_counts_field = None
energy_field_global = None

def init_taichi_fields():
    """Initialize global Taichi fields for grid"""
    global particle_positions_field, particle_attributes_field, particle_active_field
    global particle_count_field, spatial_hash_field, cell_counts_field, energy_field_global
    
    particle_positions_field = ti.Vector.field(2, dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    particle_attributes_field = ti.Vector.field(4, dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    particle_active_field = ti.field(dtype=ti.i32, shape=(MAX_PARTICLES_COMPILE,))
    particle_count_field = ti.field(dtype=ti.i32, shape=())
    
    spatial_hash_field = ti.field(dtype=ti.i32, shape=(GRID_CELLS_X_COMPILE, GRID_CELLS_Y_COMPILE, MAX_PARTICLES_PER_CELL_COMPILE))
    cell_counts_field = ti.field(dtype=ti.i32, shape=(GRID_CELLS_X_COMPILE, GRID_CELLS_Y_COMPILE))
    
    energy_field_global = ti.field(dtype=ti.f32, shape=(GRID_WIDTH_COMPILE, GRID_HEIGHT_COMPILE))

# Module-level kernels
@ti.kernel
def reset_grid_kernel():
    """Reset grid - module-level kernel"""
    # Reset all particles
    for i in range(MAX_PARTICLES_COMPILE):
        particle_active_field[i] = 0
    
    particle_count_field[None] = 0
    
    # Clear spatial hash
    for i, j in ti.ndrange(GRID_CELLS_X_COMPILE, GRID_CELLS_Y_COMPILE):
        cell_counts_field[i, j] = 0
        for k in range(MAX_PARTICLES_PER_CELL_COMPILE):
            spatial_hash_field[i, j, k] = -1
    
    # Clear energy field
    for i, j in ti.ndrange(GRID_WIDTH_COMPILE, GRID_HEIGHT_COMPILE):
        energy_field_global[i, j] = 0.0

@ti.kernel
def update_spatial_hash_kernel():
    """Update spatial hash table - module-level kernel"""
    # Clear hash table
    for i, j in ti.ndrange(GRID_CELLS_X_COMPILE, GRID_CELLS_Y_COMPILE):
        cell_counts_field[i, j] = 0
        for k in range(MAX_PARTICLES_PER_CELL_COMPILE):
            spatial_hash_field[i, j, k] = -1
    
    # Rebuild hash table
    for i in range(MAX_PARTICLES_COMPILE):
        if particle_active_field[i] == 1:
            pos = particle_positions_field[i]
            cell_x = int(pos[0] / 2.0)  # cell_size = 2.0
            cell_y = int(pos[1] / 2.0)
            
            # Wrap to grid bounds
            cell_x = cell_x % GRID_CELLS_X_COMPILE
            cell_y = cell_y % GRID_CELLS_Y_COMPILE
            
            # Add to hash table
            count = cell_counts_field[cell_x, cell_y]
            if count < MAX_PARTICLES_PER_CELL_COMPILE:
                spatial_hash_field[cell_x, cell_y, count] = i
                cell_counts_field[cell_x, cell_y] = count + 1

@ti.kernel
def apply_periodic_boundary_kernel():
    """Apply periodic boundary conditions - module-level kernel"""
    for i in range(MAX_PARTICLES_COMPILE):
        if particle_active_field[i] == 1:
            pos = particle_positions_field[i]
            # Wrap positions to grid bounds
            pos[0] = pos[0] % GRID_WIDTH_COMPILE
            pos[1] = pos[1] % GRID_HEIGHT_COMPILE
            particle_positions_field[i] = pos

@ti.kernel
def decay_energy_field_kernel(decay_rate: ti.f32):
    """Decay energy field - module-level kernel"""
    for i, j in ti.ndrange(GRID_WIDTH_COMPILE, GRID_HEIGHT_COMPILE):
        energy_field_global[i, j] *= decay_rate

@ti.data_oriented
class Grid:
    """2D periodic grid for particle simulation"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.height = config.grid_height
        self.width = config.grid_width
        self.max_particles = config.max_particles
        
        # Spatial hash parameters
        self.cell_size = 2.0  # Should be > 2 * particle_radius
        self.grid_cells_x = int(np.ceil(self.width / self.cell_size))
        self.grid_cells_y = int(np.ceil(self.height / self.cell_size))
        self.max_particles_per_cell = 32

        # Initialize global fields if not already done
        if particle_positions_field is None:
            init_taichi_fields()
        
        # Use global fields
        self.particle_positions = particle_positions_field
        self.particle_attributes = particle_attributes_field
        self.particle_active = particle_active_field
        self.particle_count = particle_count_field
        
        # Hash table: cell -> list of particle indices
        self.spatial_hash = spatial_hash_field
        self.cell_counts = cell_counts_field
        
        # Energy field
        self.energy_field = energy_field_global
        
        # Initialize
        self.reset()
    
    def reset(self):
        """Reset grid to initial state"""
        reset_grid_kernel()
    
    @ti.kernel
    def add_particle(self, pos: ti.template(), attributes: ti.template()) -> ti.i32:
        """Add a particle to the grid. Returns particle index or -1 if failed"""
        if self.particle_count[None] >= self.max_particles:  # Use configured max_particles
            return -1
        
        idx = self.particle_count[None]
        self.particle_positions[idx] = pos
        self.particle_attributes[idx] = attributes
        self.particle_active[idx] = 1
        self.particle_count[None] += 1
        
        return idx
    
    @ti.kernel
    def remove_particle(self, idx: ti.i32):
        """Remove particle by index"""
        if idx < 0 or idx >= self.particle_positions.shape[0]:
            return
        
        self.particle_active[idx] = 0
    
    def update_spatial_hash(self):
        """Update spatial hash table for neighbor finding"""
        update_spatial_hash_kernel()
    
    @ti.kernel
    def get_neighbors(self, pos: ti.template(), radius: ti.f32, neighbors: ti.template()) -> ti.i32:
        """Get neighbors within radius. Returns count of neighbors found"""
        count = 0
        wrapped_x = pos[0] % self.width
        wrapped_y = pos[1] % self.height
        
        # Calculate cell range
        cell_radius = int(radius / self.cell_size) + 1
        
        for dx in range(-cell_radius, cell_radius + 1):
            for dy in range(-cell_radius, cell_radius + 1):
                cell_x = int(wrapped_x / self.cell_size) + dx
                cell_y = int(wrapped_y / self.cell_size) + dy
                
                # Wrap cell coordinates
                cell_x = cell_x % self.grid_cells_x
                cell_y = cell_y % self.grid_cells_y
                
                # Check particles in this cell
                for k in range(self.cell_counts[cell_x, cell_y]):
                    p_idx = self.spatial_hash[cell_x, cell_y, k]
                    if p_idx >= 0 and self.particle_active[p_idx] == 1:
                        p_pos = self.particle_positions[p_idx]
                        
                        # Calculate distance with periodic boundary conditions
                        dx_pos = p_pos[0] - pos[0]
                        dy_pos = p_pos[1] - pos[1]
                        
                        # Apply periodic boundary conditions
                        if dx_pos > self.width / 2:
                            dx_pos -= self.width
                        elif dx_pos < -self.width / 2:
                            dx_pos += self.width
                        
                        if dy_pos > self.height / 2:
                            dy_pos -= self.height
                        elif dy_pos < -self.height / 2:
                            dy_pos += self.height
                        
                        distance = ti.sqrt(dx_pos * dx_pos + dy_pos * dy_pos)
                        
                        if distance <= radius and count < neighbors.shape[0]:
                            neighbors[count] = p_idx
                            count += 1
        
        return count
    
    def apply_periodic_boundary(self):
        """Apply periodic boundary conditions to all particles"""
        apply_periodic_boundary_kernel()
    
    def decay_energy_field(self, decay_rate: float):
        """Apply energy decay to the energy field"""
        decay_energy_field_kernel(decay_rate)
    
    def get_particle_data(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get particle data as numpy arrays for external use"""
        positions = self.particle_positions.to_numpy()
        attributes = self.particle_attributes.to_numpy()
        active = self.particle_active.to_numpy()
        
        # Filter only active particles
        active_mask = active == 1
        active_positions = positions[active_mask]
        active_attributes = attributes[active_mask]
        
        return active_positions, active_attributes, active_mask
    
    def get_energy_field(self) -> np.ndarray:
        """Get energy field as numpy array"""
        return self.energy_field.to_numpy()
    
    def get_stats(self) -> dict:
        """Get grid statistics"""
        return {
            "particle_count": self.particle_count[None],
            "grid_size": (self.width, self.height),
            "cell_size": self.cell_size,
            "max_particles": self.particle_positions.shape[0]
        }
