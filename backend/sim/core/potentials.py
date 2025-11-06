"""
Potential functions for Live 2.0 simulation
Handles particle-particle interactions and binding potentials

PERFORMANCE: Now uses spatial hashing for O(n) force computation!
"""

import taichi as ti
import numpy as np
from typing import Tuple, Dict, Optional
from pathlib import Path
import logging
from ..config import SimulationConfig

logger = logging.getLogger(__name__)

# Import spatial hashing
try:
    from .spatial_hash import (
        init_spatial_hash_fields, 
        build_spatial_hash, 
        compute_forces_spatial_simple,
        get_stats as get_spatial_stats
    )
    SPATIAL_HASH_AVAILABLE = True
    logger.info("Spatial hashing enabled - O(n) force computation!")
except ImportError as e:
    logger.warning(f"Spatial hashing not available: {e}")
    SPATIAL_HASH_AVAILABLE = False

# Compile-time constants
MAX_PARTICLES_COMPILE = 10000

# Global Taichi fields
forces_field = None
binding_matrix_field = None
potential_strength_field = None
potential_range_field = None
binding_threshold_field = None
unbinding_threshold_field = None

def init_potential_fields():
    """Initialize global Taichi fields for potentials"""
    global forces_field, binding_matrix_field
    global potential_strength_field, potential_range_field
    global binding_threshold_field, unbinding_threshold_field
    
    forces_field = ti.Vector.field(2, dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE,))
    binding_matrix_field = ti.field(dtype=ti.f32, shape=(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE))
    potential_strength_field = ti.field(dtype=ti.f32, shape=())
    potential_range_field = ti.field(dtype=ti.f32, shape=())
    binding_threshold_field = ti.field(dtype=ti.f32, shape=())
    unbinding_threshold_field = ti.field(dtype=ti.f32, shape=())

# Module-level kernels
@ti.kernel
def reset_forces_kernel():
    """Reset forces - module-level kernel"""
    for i in range(MAX_PARTICLES_COMPILE):
        forces_field[i] = ti.Vector([0.0, 0.0])

@ti.kernel
def reset_binding_matrix_kernel():
    """Reset binding matrix - module-level kernel"""
    for i, j in ti.ndrange(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE):
        binding_matrix_field[i, j] = 0.0

@ti.kernel
def compute_forces_kernel(positions: ti.template(), attributes: ti.template(),
                         active: ti.template(), particle_count: ti.i32):
    """Compute forces between all particle pairs - module-level kernel with stability improvements"""
    # Clear forces
    for i in range(MAX_PARTICLES_COMPILE):
        forces_field[i] = ti.Vector([0.0, 0.0])
    
    # Compute pairwise forces with improved stability
    for i in range(particle_count):
        if active[i] == 1:
            for j in range(i + 1, particle_count):
                if active[j] == 1:
                    pos_i = positions[i]
                    pos_j = positions[j]
                    
                    # Calculate distance vector
                    r_vec = pos_i - pos_j
                    r = r_vec.norm()
                    
                    # Improved distance handling with soft cutoff
                    min_distance = 0.2  # Increased from 0.1 for stability
                    max_distance = 10.0  # Add maximum interaction distance
                    
                    if r > min_distance and r < max_distance:
                        # Normalize distance vector
                        r_hat = r_vec / r
                        
                        # Get particle attributes
                        mass_i = attributes[i][0]
                        mass_j = attributes[j][0]
                        charge_i = attributes[i][1]  # x-component of charge vector
                        charge_j = attributes[j][1]  # x-component of charge vector
                        
                        # Compute forces with reduced strength for stability
                        # Lennard-Jones force (reduced strength)
                        lj_force = lennard_jones_force_func(r, 0.5, 1.0)  # Reduced epsilon
                        
                        # Coulomb force (reduced strength)
                        coulomb_force = coulomb_force_func(r, charge_i, charge_j, 0.5)  # Reduced strength
                        
                        # Total force magnitude with cap
                        total_force = lj_force + coulomb_force
                        max_force = 5.0  # Cap maximum force magnitude
                        total_force = ti.max(-max_force, ti.min(max_force, total_force))
                        
                        # Apply force to both particles (Newton's third law)
                        force_vec = total_force * r_hat
                        
                        forces_field[i] += force_vec
                        forces_field[j] -= force_vec

@ti.func
def lennard_jones_force_func(r: ti.f32, epsilon: ti.f32 = 1.0, sigma: ti.f32 = 1.0) -> ti.f32:
    """Force magnitude from Lennard-Jones potential"""
    if r < 0.1:
        r = 0.1
    
    sr6 = (sigma / r) ** 6
    sr12 = sr6 * sr6
    return 24.0 * epsilon * (2.0 * sr12 - sr6) / r

@ti.func
def coulomb_force_func(r: ti.f32, q1: ti.f32, q2: ti.f32, k: ti.f32 = 1.0) -> ti.f32:
    """Force magnitude from Coulomb potential"""
    if r < 0.1:
        r = 0.1
    
    return k * q1 * q2 / (r * r)

@ti.kernel
def update_binding_matrix_kernel(positions: ti.template(), attributes: ti.template(),
                               active: ti.template(), particle_count: ti.i32):
    """Update binding matrix - module-level kernel"""
    # Clear binding matrix
    for i, j in ti.ndrange(MAX_PARTICLES_COMPILE, MAX_PARTICLES_COMPILE):
        binding_matrix_field[i, j] = 0.0
    
    # Update binding strengths
    for i in range(particle_count):
        if active[i] == 1:
            for j in range(i + 1, particle_count):
                if active[j] == 1:
                    pos_i = positions[i]
                    pos_j = positions[j]
                    
                    r_vec = pos_i - pos_j
                    r = r_vec.norm()
                    
                    if r < potential_range_field[None]:
                        # Compute binding strength based on distance and particle properties
                        binding_strength = compute_binding_strength_func(i, j, r, attributes)
                        
                        # Update binding matrix (symmetric)
                        binding_matrix_field[i, j] = binding_strength
                        binding_matrix_field[j, i] = binding_strength

@ti.func
def compute_binding_strength_func(i: ti.i32, j: ti.i32, r: ti.f32, 
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
    
    distance_factor = ti.exp(-r / potential_range_field[None])
    mass_factor = ti.min(mass_i, mass_j) / ti.max(mass_i, mass_j)
    charge_factor = 1.0 - ti.abs(charge_i - charge_j) / ti.max(ti.abs(charge_i), ti.abs(charge_j), 1.0)
    
    binding_strength = distance_factor * mass_factor * charge_factor * potential_strength_field[None]
    
    return binding_strength

@ti.kernel
def apply_binding_forces_kernel(positions: ti.template(), attributes: ti.template(),
                              active: ti.template(), particle_count: ti.i32):
    """Apply binding forces - module-level kernel"""
    for i in range(particle_count):
        if active[i] == 1:
            for j in range(i + 1, particle_count):
                if active[j] == 1:
                    binding_strength = binding_matrix_field[i, j]
                    
                    if binding_strength > binding_threshold_field[None]:
                        pos_i = positions[i]
                        pos_j = positions[j]
                        
                        r_vec = pos_i - pos_j
                        r = r_vec.norm()
                        
                        if r > 0.1:
                            r_hat = r_vec / r
                            
                            # Binding force (attractive)
                            binding_force = binding_force_func(r, binding_strength, 1.0)
                            force_vec = binding_force * r_hat
                            
                            forces_field[i] += force_vec
                            forces_field[j] -= force_vec

@ti.func
def binding_force_func(r: ti.f32, binding_strength: ti.f32,
                      equilibrium_distance: ti.f32 = 1.0) -> ti.f32:
    """Force magnitude from binding potential"""
    if r < 0.1:
        r = 0.1
    
    dr = r - equilibrium_distance
    return binding_strength * (2.0 * dr - 1.5 * dr * dr / equilibrium_distance)

@ti.data_oriented
class PotentialSystem:
    """Manages potential functions and particle interactions"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        
        # Initialize global fields if not already done
        if forces_field is None:
            init_potential_fields()
        
        # Use global fields
        self.forces = forces_field
        self.binding_matrix = binding_matrix_field
        self.potential_strength = potential_strength_field
        self.potential_range = potential_range_field
        self.binding_threshold = binding_threshold_field
        self.unbinding_threshold = unbinding_threshold_field
        
        # Initialize parameters
        self.potential_strength[None] = 1.0
        self.potential_range[None] = 2.0
        self.binding_threshold[None] = config.binding_threshold
        self.unbinding_threshold[None] = config.unbinding_threshold
        
        # Initialize spatial hashing (O(n) performance!) - check config first
        self.use_spatial_hash = SPATIAL_HASH_AVAILABLE and config.use_spatial_hash
        self.spatial_cell_size = 10.0  # Default
        if self.use_spatial_hash:
            try:
                # Check config for cell size
                if hasattr(config, 'spatial_hash_cell_size'):
                    self.spatial_cell_size = config.spatial_hash_cell_size
                elif hasattr(config, 'spatial_cell_size'):
                    self.spatial_cell_size = config.spatial_cell_size
                
                init_spatial_hash_fields(cell_size=self.spatial_cell_size)
                logger.info(f"Spatial hashing initialized (cell_size={self.spatial_cell_size})")
            except Exception as e:
                logger.error(f"Failed to initialize spatial hashing: {e}")
                self.use_spatial_hash = False
        
        # Load Physics Database (PHASE 1 WEEK 2: Literature parameters)
        self.physics_db = None
        self.use_physics_db = config.use_physics_db
        
        if self.use_physics_db:
            try:
                from .physics_db import PhysicsDatabase
                db_path = Path(config.physics_db_path)
                
                # Fix relative paths - try multiple locations
                import os
                
                # List of possible paths to try
                possible_paths = [
                    db_path,  # Original path
                    Path(config.physics_db_path),  # Config path
                    Path.cwd() / config.physics_db_path,  # CWD relative
                    Path(__file__).parent.parent.parent / "data" / "physics_parameters.json",  # Backend relative
                    Path(__file__).parent.parent.parent.parent / "data" / "physics_parameters.json",  # Backend relative (one more level)
                    Path.cwd() / "data" / "physics_parameters.json",  # Root/data
                    Path.cwd().parent / "data" / "physics_parameters.json",  # Parent/data
                    Path.home() / "live2.0" / "data" / "physics_parameters.json",  # Home directory
                ]
                
                # Also try to find project root by looking for common markers
                current = Path.cwd()
                for _ in range(5):  # Go up max 5 levels
                    test_path = current / "data" / "physics_parameters.json"
                    if test_path.exists():
                        possible_paths.insert(0, test_path)  # Add to front of list
                        break
                    if (current / "backend").exists() and (current / "scripts").exists():
                        # Found project root
                        test_path = current / "data" / "physics_parameters.json"
                        if test_path.exists():
                            possible_paths.insert(0, test_path)
                        break
                    current = current.parent
                    if current == current.parent:  # Reached filesystem root
                        break
                
                # Find first existing path
                for possible_path in possible_paths:
                    try:
                        if possible_path.exists():
                            db_path = possible_path
                            logger.info(f"Found PhysicsDatabase at: {db_path}")
                            break
                    except:
                        continue
                
                if db_path.exists():
                    self.physics_db = PhysicsDatabase(str(db_path))
                    logger.info(f"Loaded PhysicsDatabase from {db_path}")
                    
                    stats = self.physics_db.get_statistics()
                    logger.info(f"  Bond parameters: {stats['total_bonds']}")
                    logger.info(f"  VDW parameters: {stats['total_vdw']}")
                    logger.info(f"  Citations: {stats['unique_citations']}")
                else:
                    logger.error(f"PhysicsDatabase not found, using fallback parameters")
                    self.use_physics_db = False
                    
            except Exception as e:
                logger.error(f"Failed to load PhysicsDatabase: {e}")
                logger.info("Using fallback parameters from config")
                self.use_physics_db = False
        
        # Store fallback parameters for easy access
        self.default_epsilon = config.default_epsilon
        self.default_sigma = config.default_sigma
        self.default_bond_D_e = config.default_bond_D_e
        self.default_bond_r_e = config.default_bond_r_e
        self.default_bond_a = config.default_bond_a
        
        # Initialize
        self.reset()
    
    def get_vdw_parameters(self, atom_a: str = 'C', atom_b: str = 'C') -> Tuple[float, float]:
        """
        Get Van der Waals parameters (epsilon, sigma) for atom pair.
        
        Returns:
            (epsilon, sigma) in (kJ/mol, Angstrom)
        """
        if self.use_physics_db and self.physics_db is not None:
            try:
                epsilon, sigma = self.physics_db.get_vdw_parameters(atom_a, atom_b)
                return (epsilon, sigma)
            except Exception as e:
                logger.warning(f"Failed to get VDW parameters for {atom_a}-{atom_b}: {e}")
        
        # Fallback
        return (self.default_epsilon, self.default_sigma)
    
    def get_bond_parameters(self, atom_a: str = 'C', atom_b: str = 'C', order: int = 1) -> Tuple[float, float, float]:
        """
        Get bond parameters (D_e, r_e, a) for atom pair.
        
        Returns:
            (D_e, r_e, a) in (kJ/mol, Angstrom, 1/Angstrom)
        """
        if self.use_physics_db and self.physics_db is not None:
            try:
                params = self.physics_db.get_bond_parameters(atom_a, atom_b, order)
                if params is not None:
                    return (params.D_e, params.r_e, params.a)
            except Exception as e:
                logger.warning(f"Failed to get bond parameters for {atom_a}-{atom_b} order {order}: {e}")
        
        # Fallback
        return (self.default_bond_D_e, self.default_bond_r_e, self.default_bond_a)
    
    def reset(self):
        """Reset potential system"""
        reset_forces_kernel()
        reset_binding_matrix_kernel()
    
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
    def morse_potential(self, r: ti.f32, D_e: ti.f32, r_e: ti.f32, a: ti.f32) -> ti.f32:
        """
        Morse potential: V(r) = D_e * (1 - exp(-a*(r - r_e)))²
        
        Used for chemical bonds with literature parameters.
        
        Args:
            r: Distance (Angstrom)
            D_e: Dissociation energy (kJ/mol)
            r_e: Equilibrium distance (Angstrom)
            a: Width parameter (1/Angstrom)
        """
        if r < 0.1:
            r = 0.1
        
        exp_term = ti.exp(-a * (r - r_e))
        return D_e * (1.0 - exp_term) ** 2
    
    @ti.func
    def morse_force(self, r: ti.f32, D_e: ti.f32, r_e: ti.f32, a: ti.f32) -> ti.f32:
        """
        Force magnitude from Morse potential: F = -dV/dr
        
        F = 2*a*D_e*(1 - exp(-a*(r-r_e)))*exp(-a*(r-r_e))
        """
        if r < 0.1:
            r = 0.1
        
        exp_term = ti.exp(-a * (r - r_e))
        return 2.0 * a * D_e * (1.0 - exp_term) * exp_term
    
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
    
    def compute_forces(self, positions, attributes, active, particle_count: int):
        """
        Compute forces between particles.
        
        Now uses spatial hashing for O(n) performance!
        Falls back to O(n²) if spatial hash not available.
        """
        if self.use_spatial_hash:
            # O(n) spatial hashing approach!
            try:
                # Build spatial hash grid
                box_width = self.config.grid_width
                box_height = self.config.grid_height
                build_spatial_hash(positions, active, particle_count, box_width, box_height)
                
                # Compute forces using spatial hash
                compute_forces_spatial_simple(positions, attributes, active, particle_count, self.forces)
            except Exception as e:
                logger.error(f"Spatial hash failed: {e}, falling back to O(n²)")
                compute_forces_kernel(positions, attributes, active, particle_count)
        else:
            # O(n²) fallback
            compute_forces_kernel(positions, attributes, active, particle_count)
    
    def update_binding_matrix(self, positions, attributes, active, particle_count: int):
        """Update binding matrix based on particle distances and properties"""
        update_binding_matrix_kernel(positions, attributes, active, particle_count)
    
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
    
    def apply_binding_forces(self, positions, attributes, active, particle_count: int):
        """Apply binding forces based on current binding matrix"""
        apply_binding_forces_kernel(positions, attributes, active, particle_count)
    
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
