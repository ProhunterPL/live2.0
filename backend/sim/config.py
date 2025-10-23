# Live 2.0 Configuration
from typing import Dict, Any, Optional
from pydantic import BaseModel, Field
import numpy as np

class SimulationConfig(BaseModel):
    """Main configuration for Live 2.0 simulation"""
    
    model_config = {"extra": "ignore"}  # Ignore extra fields from frontend
    
    # Grid settings
    grid_height: int = Field(default=256, ge=64, le=1024)
    grid_width: int = Field(default=256, ge=64, le=1024)
    
    # Simulation mode
    mode: str = Field(default="open_chemistry", pattern="^(preset_prebiotic|open_chemistry)$")
    
    # Time settings - STABILIZED for particle retention
    dt: float = Field(default=0.005, gt=0, le=1.0)  # Even smaller timestep for stability
    max_time: float = Field(default=1000.0, gt=0)
    
    # Energy settings - SCIENTIFICALLY CALIBRATED for prebiotic chemistry
    # Based on Miller-Urey (1953): ~1-10 eV per molecule in discharge zone
    # FIXED: Reduced pulse amplitude and frequency for realistic energy control
    energy_decay: float = Field(default=0.90, gt=0, lt=1)  # Faster decay for more controlled energy
    energy_threshold: float = Field(default=0.1, gt=0)
    pulse_every: int = Field(default=150, gt=0)  # Less frequent pulses (was 80) - Miller-Urey continuous but controlled
    pulse_radius: float = Field(default=15.0, gt=0)  # Larger pulse radius (was 12.0)
    pulse_amplitude: float = Field(default=0.8, gt=0)  # FURTHER REDUCED from 1.8 - energy too high (1.82M vs target ~150k)
    diffuse_D: float = Field(default=0.25, gt=0)  # Faster diffusion for more encounters
    target_energy: float = Field(default=0.3, gt=0)  # Lower background energy (was 0.5) - more realistic
    thermostat_alpha: float = Field(default=0.15, gt=0)  # VERY STRONG thermostat (was 0.05) - energy 11× too high, need aggressive control
    
    # Particle settings
    max_particles: int = Field(default=10000, gt=0, le=100000)
    particle_radius: float = Field(default=0.5, gt=0)
    
    # Binding settings - SCIENTIFICALLY CALIBRATED based on bond energies
    # Literature: vdW bonds 2-10 kJ/mol, H-bonds 10-40 kJ/mol, covalent 300-400 kJ/mol
    # FIXED: More permissive thresholds for realistic cluster formation (Miller-Urey experiments)
    binding_threshold: float = Field(default=0.4, gt=0, le=1)  # REDUCED from 0.7 - easier bond formation
    unbinding_threshold: float = Field(default=0.15, gt=0, le=1)  # REDUCED from 0.2 - more stable bonds  # Kept for compatibility
    
    # Novelty detection
    novelty_window: int = Field(default=100, gt=0)
    min_cluster_size: int = Field(default=2, ge=1)  # CHANGED back to 2 - only detect real clusters (not single particles)
    novelty_check_interval: int = Field(default=500, gt=0)  # How often to check for novel substances (steps) - BALANCED for performance
    
    # Visualization
    vis_frequency: int = Field(default=5, gt=0)  # Zmniejszone z 3 dla stabilności
    log_frequency: int = Field(default=100, gt=0)
    
    # Thermodynamic validation (SCIENTIFIC RIGOR WITH SAFETY)
    validate_every_n_steps: int = Field(default=300, gt=0)  # REDUCED frequency for better performance (was 150)
    enable_thermodynamic_validation: bool = Field(default=True)  # Enabled by default
    energy_tolerance: float = Field(default=2e-3, gt=0)  # Balanced tolerance for scientific accuracy
    momentum_tolerance: float = Field(default=2e-4, gt=0)  # Balanced tolerance for scientific accuracy
    
    # Mutations - INCREASED for more diversity and novelty
    p_mut_base: float = Field(default=1e-3, gt=0)  # Increased mutations (2x more)
    p_mut_gain: float = Field(default=30.0, gt=0)  # Increased from 20.0
    attr_sigma: float = Field(default=0.15, gt=0)  # Stronger perturbation
    
    # Performance optimization parameters
    energy_update_interval: int = Field(default=5, gt=0, description="Update energy every N steps")
    metrics_update_interval: int = Field(default=1, gt=0, description="Update metrics every N steps")
    
    # Diagnostics
    enable_diagnostics: bool = Field(default=True)
    diagnostics_dir: str = Field(default="diagnostics")
    diagnostics_frequency: int = Field(default=10, gt=0, description="Log diagnostics every N steps")
    
    # Memory management (based on NAMD/GROMACS best practices)
    memory_check_interval: int = Field(default=1000, gt=0, description="Check memory usage every N steps")
    gc_interval: int = Field(default=5000, gt=0, description="Run Python garbage collection every N steps")
    
    # Physics Database (PHASE 1 WEEK 2: Literature-based parameters)
    use_physics_db: bool = Field(default=True, description="Use literature parameters from physics database")
    physics_db_path: str = Field(default="../data/physics_parameters.json", description="Path to physics parameters database")
    
    # Fallback parameters (when DB not available or parameter not found)
    default_epsilon: float = Field(default=0.439, gt=0, description="Default LJ epsilon (kJ/mol) - Carbon UFF value")
    default_sigma: float = Field(default=3.431, gt=0, description="Default LJ sigma (Angstrom) - Carbon UFF value")
    default_bond_D_e: float = Field(default=348.0, gt=0, description="Default bond dissociation energy (kJ/mol) - C-C single")
    default_bond_r_e: float = Field(default=1.54, gt=0, description="Default bond equilibrium length (Angstrom) - C-C single")
    default_bond_a: float = Field(default=1.8, gt=0, description="Default Morse width parameter (1/Angstrom) - C-C single")
    
    # Random seed
    seed: Optional[int] = Field(default=42)

class PresetPrebioticConfig(BaseModel):
    """Configuration for Preset Prebiotic mode"""
    
    # Chemical species
    species: Dict[str, float] = Field(default={
        "HCN": 0.1,
        "NH2CHO": 0.0,
        "H2O": 0.5
    })
    
    # Reaction rates
    reaction_rates: Dict[str, float] = Field(default={
        "HCN_to_NH2CHO": 0.01
    })
    
    # Diffusion coefficients
    diffusion_coeffs: Dict[str, float] = Field(default={
        "HCN": 0.1,
        "NH2CHO": 0.05,
        "H2O": 0.2
    })

class BondTypeConfig(BaseModel):
    """Configuration for a bond type"""
    name: str = Field(default="")
    k_spring: float = Field(default=2.0, gt=0, description="Spring constant")
    rest_len_factor: float = Field(default=1.0, gt=0, description="Rest length as factor of formation distance")
    damping: float = Field(default=0.1, ge=0, description="Damping coefficient")
    strength: float = Field(default=5.0, gt=0, description="Breaking strength threshold")
    max_age: float = Field(default=100.0, gt=0, description="Maximum age before probabilistic breaking")

class OpenChemistryConfig(BaseModel):
    """Configuration for Open Chemistry mode"""
    
    # Potential parameters
    potential_strength: float = Field(default=1.0, gt=0)
    potential_range: float = Field(default=2.0, gt=0)
    
    # Mutation settings
    mutation_rate: float = Field(default=0.001, ge=0, le=1)
    mutation_strength: float = Field(default=0.1, gt=0)
    
    # Energy sources
    energy_sources: int = Field(default=3, ge=0)
    energy_intensity: float = Field(default=5.0, gt=0)
    
    # Bond system configuration - SCIENTIFICALLY CALIBRATED
    # Literature: Bond lifetimes 10^-9 s (vdW) to 10^6 s (covalent) at 298K
    # FIXED: More permissive binding conditions for realistic cluster formation
    enable_spring_forces: bool = Field(default=True, description="Enable spring forces from bonds")
    enable_advanced_breaking: bool = Field(default=True, description="Enable overload/aging bond breaking")
    theta_bind: float = Field(default=0.3, gt=0)  # REDUCED from 0.5 - easier bond formation (Miller-Urey)
    theta_break: float = Field(default=3.0, gt=0)  # MORE STABLE bonds (was 2.0) - prevent high energy from breaking clusters
    vmax: float = Field(default=6.0, gt=0)  # Lowered from 8.0 - less violent collisions
    neighbor_radius: float = Field(default=2.5, gt=0)  # Increased from 2.0 - detect more neighbors
    rebuild_neighbors_every: int = Field(default=15, gt=0)  # Less frequent updates (was 10)
    clamp_density_per_cell: int = Field(default=16, gt=0)  # Much lower density limit (was 32)
    bond_types: Dict[int, Dict[str, float]] = Field(default={
        0: {'name': 'vdW', 'k_spring': 3.0, 'rest_len_factor': 1.0, 'damping': 0.15, 'strength': 8.0, 'max_age': 100.0},
        1: {'name': 'covalent', 'k_spring': 15.0, 'rest_len_factor': 0.9, 'damping': 0.25, 'strength': 30.0, 'max_age': 400.0},
        2: {'name': 'H-bond', 'k_spring': 7.0, 'rest_len_factor': 1.1, 'damping': 0.2, 'strength': 15.0, 'max_age': 150.0},
        3: {'name': 'metallic', 'k_spring': 10.0, 'rest_len_factor': 1.0, 'damping': 0.3, 'strength': 20.0, 'max_age': 300.0}
    }, description="Bond type parameters - INCREASED strength and longevity")
    
    # Cluster configuration
    min_cluster_size: int = Field(default=2, ge=1, description="Minimum cluster size to track")
    min_cluster_density: float = Field(default=0.15, ge=0, le=1, description="Minimum graph density to consider valid cluster")
    compute_cluster_R_g: bool = Field(default=True, description="Compute radius of gyration for clusters")

def get_default_config() -> SimulationConfig:
    """Get default simulation configuration"""
    return SimulationConfig()

def get_preset_config() -> PresetPrebioticConfig:
    """Get default preset prebiotic configuration"""
    return PresetPrebioticConfig()

def get_open_chemistry_config() -> OpenChemistryConfig:
    """Get default open chemistry configuration"""
    return OpenChemistryConfig()
