"""
Grid management for Live 2.0 simulation
Handles 2D periodic grid with particle positioning and neighbor finding
"""

import taichi as ti
import numpy as np
from typing import Tuple, List
from ..config import SimulationConfig

# Global variables for Taichi fields (will be initialized later)
particle_positions = None
particle_attributes = None
particle_active = None
particle_count = None
spatial_hash = None
cell_counts = None
energy_field = None

def init_taichi_fields(max_particles: int, grid_cells_x: int, grid_cells_y: int,
                       energy_width: int, energy_height: int, max_particles_per_cell: int = 32):
    """Initialize Taichi fields after ti.init() has been called"""
    global particle_positions, particle_attributes, particle_active, particle_count
    global spatial_hash, cell_counts, energy_field
    
    particle_positions = ti.Vector.field(2, dtype=ti.f32, shape=(max_particles,))
    particle_attributes = ti.Vector.field(4, dtype=ti.f32, shape=(max_particles,))
    particle_active = ti.field(dtype=ti.i32, shape=(max_particles,))
    particle_count = ti.field(dtype=ti.i32, shape=())
    
    spatial_hash = ti.field(dtype=ti.i32, shape=(grid_cells_x, grid_cells_y, max_particles_per_cell))
    cell_counts = ti.field(dtype=ti.i32, shape=(grid_cells_x, grid_cells_y))
    
    energy_field = ti.field(dtype=ti.f32, shape=(energy_width, energy_height))

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

        # Initialize Taichi fields if not already done
        if particle_positions is None:
            init_taichi_fields(
                max_particles=self.max_particles,
                grid_cells_x=self.grid_cells_x,
                grid_cells_y=self.grid_cells_y,
                energy_width=self.width,
                energy_height=self.height,
                max_particles_per_cell=self.max_particles_per_cell,
            )
        
        # Use module-level fields
        self.particle_positions = particle_positions
        self.particle_attributes = particle_attributes
        self.particle_active = particle_active
        self.particle_count = particle_count
        
        # Hash table: cell -> list of particle indices
        self.spatial_hash = spatial_hash
        self.cell_counts = cell_counts
        
        # Energy field
        self.energy_field = energy_field
        
        # Initialize
        self.reset()
    
    @ti.kernel
    def reset(self):
        """Reset grid to initial state"""
        # Reset all particles
        for i in range(self.max_particles):
            self.particle_active[i] = 0
        
        self.particle_count[None] = 0
        
        # Clear spatial hash
        for i, j in ti.ndrange(self.grid_cells_x, self.grid_cells_y):
            self.cell_counts[i, j] = 0
            for k in range(self.max_particles_per_cell):
                self.spatial_hash[i, j, k] = -1
        
        # Clear energy field
        for i, j in ti.ndrange(self.width, self.height):
            self.energy_field[i, j] = 0.0
    
    @ti.kernel
    def add_particle(self, pos: ti.template(), attributes: ti.template()) -> ti.i32:
        """Add a particle to the grid. Returns particle index or -1 if failed"""
        if self.particle_count[None] >= 10000:  # Use fixed max_particles
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
    
    @ti.kernel
    def update_spatial_hash(self):
        """Update spatial hash table for neighbor finding"""
        # Clear hash table
        for i, j in ti.ndrange(self.grid_cells_x, self.grid_cells_y):
            self.cell_counts[i, j] = 0
            for k in range(self.max_particles_per_cell):
                self.spatial_hash[i, j, k] = -1
        
        # Fill hash table
        for p in range(self.particle_positions.shape[0]):
            if self.particle_active[p] == 1:
                pos = self.particle_positions[p]
                
                # Wrap position to grid bounds
                wrapped_x = pos[0] % self.width
                wrapped_y = pos[1] % self.height
                
                # Calculate cell coordinates
                cell_x = int(wrapped_x / self.cell_size)
                cell_y = int(wrapped_y / self.cell_size)
                
                # Clamp to valid range
                cell_x = ti.max(0, ti.min(cell_x, self.grid_cells_x - 1))
                cell_y = ti.max(0, ti.min(cell_y, self.grid_cells_y - 1))
                
                # Add to hash table
                count = self.cell_counts[cell_x, cell_y]
                if count < self.max_particles_per_cell:
                    self.spatial_hash[cell_x, cell_y, count] = p
                    self.cell_counts[cell_x, cell_y] += 1
    
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
    
    @ti.kernel
    def apply_periodic_boundary(self):
        """Apply periodic boundary conditions to all particles"""
        for i in range(self.particle_positions.shape[0]):
            if self.particle_active[i] == 1:
                pos = self.particle_positions[i]
                self.particle_positions[i] = ti.Vector([
                    pos[0] % self.width,
                    pos[1] % self.height
                ])
    
    @ti.kernel
    def decay_energy_field(self, decay_rate: ti.f32):
        """Apply energy decay to the energy field"""
        for i, j in ti.ndrange(self.width, self.height):
            self.energy_field[i, j] *= decay_rate
    
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
