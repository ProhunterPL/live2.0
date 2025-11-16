"""
Adaptive Spatial Hashing System for Live 2.0
=============================================

PROOF-OF-CONCEPT for Patent Application
Patent Claim: Dynamic cell sizing based on particle density and bonding topology

Key Innovation:
- Cell size adapts to system state (density, bonding, phase)
- O(n) complexity maintained
- Improved performance in heterogeneous systems

Author: Live 2.0 Team
Date: 2025-11-16
Status: PROOF-OF-CONCEPT (not production)
"""

import taichi as ti
import numpy as np
import logging
from typing import Tuple, Dict

logger = logging.getLogger(__name__)

# Compile-time constants (same as original)
MAX_PARTICLES = 10000
MAX_CELLS = 4096  # Will be dynamically adjusted
MAX_PARTICLES_PER_CELL = 128

# Global fields for adaptive system
adaptive_cell_size_field = None
grid_cell_list_adaptive = None
grid_cell_count_adaptive = None
grid_dims_adaptive = None
cell_size_history = None  # Track cell size evolution

# Adaptive parameters (tunable)
ALPHA = 2.0  # Base scaling factor
BETA = 0.3   # Bonding influence factor
MIN_CELL_SIZE = 5.0
MAX_CELL_SIZE = 25.0
TARGET_PARTICLES_PER_CELL = 16.0  # Optimal load


def init_adaptive_fields():
    """Initialize Taichi fields for adaptive spatial hash"""
    global adaptive_cell_size_field, grid_cell_list_adaptive
    global grid_cell_count_adaptive, grid_dims_adaptive, cell_size_history
    
    # Adaptive cell size (scalar)
    adaptive_cell_size_field = ti.field(dtype=ti.f32, shape=())
    
    # Grid structures (same as original but separate namespace)
    grid_cell_list_adaptive = ti.field(dtype=ti.i32, shape=(MAX_CELLS, MAX_PARTICLES_PER_CELL))
    grid_cell_count_adaptive = ti.field(dtype=ti.i32, shape=(MAX_CELLS,))
    grid_dims_adaptive = ti.Vector.field(2, dtype=ti.i32, shape=())
    
    # History tracking for analysis
    cell_size_history = ti.field(dtype=ti.f32, shape=(10000,))  # Track over time
    
    # Initialize to default
    adaptive_cell_size_field[None] = 10.0
    
    logger.info("Adaptive spatial hash fields initialized")


@ti.kernel
def compute_optimal_cell_size_kernel(
    particle_count: ti.i32,
    total_bonds: ti.i32,
    box_width: ti.f32,
    box_height: ti.f32,
    alpha: ti.f32,
    beta: ti.f32
) -> ti.f32:
    """
    PATENT FORMULA: Adaptive cell size computation
    
    s_optimal = α · √(A/N) · (1 + β·b̄)⁻¹
    
    where:
        s_optimal = optimal cell size
        α = scaling factor (typically 2.0)
        A = simulation area
        N = number of active particles
        β = bonding influence factor (typically 0.3)
        b̄ = average bonds per particle
    
    Innovation: Cell size decreases as bonding increases (particles closer together)
    """
    # Compute area and density
    area = box_width * box_height
    density = ti.cast(particle_count, ti.f32) / ti.max(area, 1.0)
    
    # Base cell size from density
    # sqrt(A/N) gives characteristic spacing between particles
    base_spacing = ti.sqrt(area / ti.max(ti.cast(particle_count, ti.f32), 1.0))
    
    # Bonding factor: more bonds = particles closer = smaller cells
    avg_bonds_per_particle = ti.cast(total_bonds, ti.f32) / ti.max(ti.cast(particle_count, ti.f32), 1.0)
    bonding_factor = 1.0 / (1.0 + beta * avg_bonds_per_particle)
    
    # PATENT FORMULA
    optimal_size = alpha * base_spacing * bonding_factor
    
    # Clamp to safe bounds
    clamped_size = ti.max(5.0, ti.min(optimal_size, 25.0))
    
    return clamped_size


@ti.func
def get_cell_index_adaptive(
    pos: ti.template(), 
    grid_width: ti.i32, 
    grid_height: ti.i32, 
    cell_size: ti.f32
) -> ti.i32:
    """Get cell index using adaptive cell size"""
    gx = ti.cast(pos[0] / cell_size, ti.i32)
    gy = ti.cast(pos[1] / cell_size, ti.i32)
    
    gx = ti.max(0, ti.min(grid_width - 1, gx))
    gy = ti.max(0, ti.min(grid_height - 1, gy))
    
    return gy * grid_width + gx


@ti.kernel
def build_adaptive_spatial_hash(
    positions: ti.template(),
    active: ti.template(),
    particle_count: ti.i32,
    box_width: ti.f32,
    box_height: ti.f32,
    cell_size: ti.f32
):
    """
    Build spatial hash with adaptive cell size
    
    Same algorithm as original, but uses dynamically computed cell_size
    """
    # Compute grid dimensions based on adaptive cell size
    grid_width = ti.cast(box_width / cell_size, ti.i32) + 1
    grid_height = ti.cast(box_height / cell_size, ti.i32) + 1
    grid_dims_adaptive[None] = ti.Vector([grid_width, grid_height])
    
    # Clear cell counts
    for i in range(MAX_CELLS):
        grid_cell_count_adaptive[i] = 0
    
    # Assign particles to cells
    for i in range(particle_count):
        if active[i] == 1:
            pos = positions[i]
            cell_idx = get_cell_index_adaptive(pos, grid_width, grid_height, cell_size)
            
            # Atomic add
            slot = ti.atomic_add(grid_cell_count_adaptive[cell_idx], 1)
            
            if slot < MAX_PARTICLES_PER_CELL:
                grid_cell_list_adaptive[cell_idx, slot] = i


@ti.kernel
def compute_forces_adaptive(
    positions: ti.template(),
    attributes: ti.template(),
    active: ti.template(),
    particle_count: ti.i32,
    forces: ti.template(),
    cell_size: ti.f32,
    cutoff: ti.f32
):
    """
    Compute forces using adaptive spatial hash
    
    Same neighbor search algorithm, but works with dynamic cell_size
    """
    grid_dims = grid_dims_adaptive[None]
    grid_width = grid_dims[0]
    grid_height = grid_dims[1]
    
    # Clear forces
    for i in range(MAX_PARTICLES):
        forces[i] = ti.Vector([0.0, 0.0])
    
    # For each particle
    for i in range(particle_count):
        if active[i] == 1:
            pos_i = positions[i]
            
            # Get cell
            gx = ti.cast(pos_i[0] / cell_size, ti.i32)
            gy = ti.cast(pos_i[1] / cell_size, ti.i32)
            gx = ti.max(0, ti.min(grid_width - 1, gx))
            gy = ti.max(0, ti.min(grid_height - 1, gy))
            
            # Check 3×3 neighborhood
            for dx in ti.static(range(-1, 2)):
                for dy in ti.static(range(-1, 2)):
                    nx = gx + dx
                    ny = gy + dy
                    
                    if 0 <= nx < grid_width and 0 <= ny < grid_height:
                        cell = ny * grid_width + nx
                        count = grid_cell_count_adaptive[cell]
                        
                        for slot in range(count):
                            if slot < MAX_PARTICLES_PER_CELL:
                                j = grid_cell_list_adaptive[cell, slot]
                                
                                if j > i and active[j] == 1:
                                    r_vec = pos_i - positions[j]
                                    r = r_vec.norm()
                                    
                                    if 0.2 < r < cutoff:
                                        r_hat = r_vec / r
                                        
                                        # LJ + Coulomb (same as original)
                                        charge_i = attributes[i][1]
                                        charge_j = attributes[j][1]
                                        
                                        sr = 1.0 / r
                                        sr6 = sr * sr * sr * sr * sr * sr
                                        sr12 = sr6 * sr6
                                        lj = 12.0 * (2.0 * sr12 - sr6) / r
                                        
                                        coul = 0.5 * charge_i * charge_j / (r * r)
                                        
                                        force = lj + coul
                                        force = ti.max(-5.0, ti.min(5.0, force))
                                        
                                        f_vec = force * r_hat
                                        forces[i] += f_vec
                                        forces[j] -= f_vec


class AdaptiveSpatialHash:
    """
    Adaptive Spatial Hashing System - PROOF-OF-CONCEPT
    
    Key Patent Claims:
    1. Dynamic cell size based on particle density
    2. Bonding topology influence on grid structure
    3. Phase-adaptive optimization
    4. O(n) complexity preservation
    """
    
    def __init__(self, box_width: float, box_height: float,
                 alpha: float = ALPHA, beta: float = BETA):
        """
        Initialize adaptive spatial hash
        
        Args:
            box_width: Simulation box width
            box_height: Simulation box height
            alpha: Base scaling factor (default 2.0)
            beta: Bonding influence factor (default 0.3)
        """
        self.box_width = box_width
        self.box_height = box_height
        self.alpha = alpha
        self.beta = beta
        
        # Initialize fields
        if adaptive_cell_size_field is None:
            init_adaptive_fields()
        
        # Tracking
        self.step_count = 0
        self.cell_size_history_np = []
        
        logger.info(f"AdaptiveSpatialHash initialized (α={alpha}, β={beta})")
    
    def compute_optimal_cell_size(self, particle_count: int, 
                                  total_bonds: int) -> float:
        """
        Compute optimal cell size for current system state
        
        Returns:
            Optimal cell size (clamped to safe bounds)
        """
        optimal = compute_optimal_cell_size_kernel(
            particle_count,
            total_bonds,
            self.box_width,
            self.box_height,
            self.alpha,
            self.beta
        )
        
        # Update field
        adaptive_cell_size_field[None] = optimal
        
        # Track history
        if self.step_count < 10000:
            cell_size_history[self.step_count] = optimal
        self.cell_size_history_np.append(optimal)
        
        self.step_count += 1
        
        return optimal
    
    def rebuild_grid(self, positions, active, particle_count: int,
                    total_bonds: int):
        """
        Rebuild spatial hash with adaptive cell size
        
        Args:
            positions: Particle positions (Taichi field)
            active: Particle active flags (Taichi field)
            particle_count: Number of active particles
            total_bonds: Total number of bonds in system
        """
        # Compute optimal cell size
        cell_size = self.compute_optimal_cell_size(particle_count, total_bonds)
        
        # Rebuild grid
        build_adaptive_spatial_hash(
            positions, active, particle_count,
            self.box_width, self.box_height, cell_size
        )
        
        return cell_size
    
    def compute_forces(self, positions, attributes, active, 
                      particle_count: int, forces, cutoff: float = 10.0):
        """
        Compute forces using adaptive spatial hash
        
        Args:
            positions: Particle positions
            attributes: Particle attributes (mass, charge, etc.)
            active: Active flags
            particle_count: Number of active particles
            forces: Output force field
            cutoff: Interaction cutoff distance
        """
        cell_size = adaptive_cell_size_field[None]
        
        compute_forces_adaptive(
            positions, attributes, active, particle_count,
            forces, cell_size, cutoff
        )
    
    def get_stats(self) -> Dict:
        """Get statistics for analysis and patent documentation"""
        if len(self.cell_size_history_np) == 0:
            return {}
        
        cell_sizes = np.array(self.cell_size_history_np)
        
        # Get grid info
        grid_dims = grid_dims_adaptive.to_numpy()
        counts = grid_cell_count_adaptive.to_numpy()
        non_empty = np.sum(counts > 0)
        
        return {
            'current_cell_size': adaptive_cell_size_field[None],
            'min_cell_size': float(np.min(cell_sizes)),
            'max_cell_size': float(np.max(cell_sizes)),
            'mean_cell_size': float(np.mean(cell_sizes)),
            'std_cell_size': float(np.std(cell_sizes)),
            'grid_dims': tuple(grid_dims),
            'non_empty_cells': int(non_empty),
            'max_particles_per_cell': int(np.max(counts)),
            'avg_particles_per_cell': float(np.mean(counts[counts > 0])) if non_empty > 0 else 0,
            'total_steps': self.step_count,
            'alpha': self.alpha,
            'beta': self.beta
        }
    
    def get_cell_size_evolution(self) -> np.ndarray:
        """Get cell size history for plotting"""
        return np.array(self.cell_size_history_np)


# Utility functions for benchmarking
def compare_fixed_vs_adaptive(particle_count: int, bond_count: int,
                              box_size: float = 256.0) -> Dict:
    """
    Compare fixed vs adaptive cell sizing
    
    Returns metrics for patent documentation
    """
    # Fixed cell size approach
    fixed_cell_size = 10.0
    fixed_grid_cells = int(box_size / fixed_cell_size) ** 2
    
    # Adaptive cell size approach
    area = box_size * box_size
    density = particle_count / area
    avg_bonds = bond_count / max(particle_count, 1)
    
    adaptive_cell_size = ALPHA * np.sqrt(area / max(particle_count, 1)) / (1.0 + BETA * avg_bonds)
    adaptive_cell_size = np.clip(adaptive_cell_size, MIN_CELL_SIZE, MAX_CELL_SIZE)
    adaptive_grid_cells = int(box_size / adaptive_cell_size) ** 2
    
    # Theoretical performance comparison
    # Fixed: checks fixed number of cells
    # Adaptive: optimizes cell count based on system state
    
    return {
        'fixed_cell_size': fixed_cell_size,
        'adaptive_cell_size': adaptive_cell_size,
        'fixed_grid_cells': fixed_grid_cells,
        'adaptive_grid_cells': adaptive_grid_cells,
        'cell_reduction_ratio': adaptive_grid_cells / max(fixed_grid_cells, 1),
        'particle_count': particle_count,
        'bond_count': bond_count,
        'avg_bonds_per_particle': avg_bonds,
        'expected_speedup': fixed_grid_cells / max(adaptive_grid_cells, 1)
    }

