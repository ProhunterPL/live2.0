# Live 2.0 Configuration
from typing import Dict, Optional
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
    pulse_every: int = Field(default=100, gt=0)  # SCIENTIFIC: Based on Miller-Urey continuous discharge (was 50, too frequent)
    pulse_radius: float = Field(default=15.0, gt=0)  # Larger pulse radius (was 12.0)
    pulse_amplitude: float = Field(default=1.8, gt=0)  # SCIENTIFIC: Based on Miller-Urey discharge energy (~1-10 eV per molecule)
    diffuse_D: float = Field(default=0.25, gt=0)  # Faster diffusion for more encounters
    target_energy: float = Field(default=0.5, gt=0)  # BALANCED: Increased from 0.3, but not as aggressive as 0.8 - prevents energy drift
    # AGGRESSIVE OPTION: target_energy=0.8 for maximum cluster formation
    thermostat_alpha: float = Field(default=0.2, gt=0)  # BALANCED: Increased from 0.1 - better energy control, prevents drift
    # AGGRESSIVE OPTION: thermostat_alpha=0.1 for maximum cluster stability
    
    # Particle settings - REDUCED for performance
    max_particles: int = Field(default=500, gt=0, le=100000)  # REDUCED from 10000 to 500 for stability
    particle_radius: float = Field(default=0.5, gt=0)
    
    # Binding settings - BALANCED for meaningful clusters
    # Literature: vdW bonds 2-10 kJ/mol, H-bonds 10-40 kJ/mol, covalent 300-400 kJ/mol
    # FIXED: Balanced thresholds for good cluster formation
    binding_threshold: float = Field(default=0.5, gt=0, le=1)  # BALANCED: Medium permissiveness
    unbinding_threshold: float = Field(default=0.18, gt=0, le=1)  # BALANCED: Stable bonds but not too restrictive
    
    # Novelty detection - BALANCED ANTI-BURNOUT settings
    novelty_window: int = Field(default=500, gt=0)  # INCREASED from 100 - longer memory for novelty detection
    min_cluster_size: int = Field(default=3, ge=1)  # INCREASED from 2 - detect meaningful clusters (3+ particles)
    novelty_check_interval: int = Field(default=500, gt=0)  # How often to check for novel substances (steps) - BALANCED for performance
    detect_novel_substances: bool = Field(default=True, description="Enable/disable novelty detection during simulation")
    
    # Visualization - OPTIMIZED for performance
    vis_frequency: int = Field(default=10, gt=0)  # INCREASED from 5 - less frequent updates
    log_frequency: int = Field(default=200, gt=0)  # INCREASED from 100 - less frequent logging
    
    # Thermodynamic validation (SCIENTIFIC RIGOR WITH SAFETY)
    validate_every_n_steps: int = Field(default=300, gt=0)  # REDUCED frequency for better performance (was 150)
    enable_thermodynamic_validation: bool = Field(default=True)  # Enabled by default
    energy_tolerance: float = Field(default=2e-3, gt=0)  # Balanced tolerance for scientific accuracy
    momentum_tolerance: float = Field(default=2e-4, gt=0)  # Balanced tolerance for scientific accuracy
    
    # Mutations - BALANCED ANTI-BURNOUT settings for sustained novelty
    p_mut_base: float = Field(default=3e-3, gt=0)  # BALANCED: Increased from 1e-3, but less aggressive than 5e-3
    p_mut_gain: float = Field(default=30.0, gt=0)  # Increased from 20.0
    attr_sigma: float = Field(default=0.15, gt=0)  # Stronger perturbation
    mutation_interval: int = Field(default=500, gt=0)  # More frequent mutations (was 2000)
    enable_mutations: bool = Field(default=True, description="Enable/disable mutations (disable on CPU/Windows to avoid LLVM errors)")
    
    # Performance optimization parameters - OPTIMIZED for stability
    energy_update_interval: int = Field(default=10, gt=0, description="Update energy every N steps")  # INCREASED from 5
    metrics_update_interval: int = Field(default=5, gt=0, description="Update metrics every N steps")  # INCREASED from 1
    
    # Hybrid GPU+CPU mode settings
    chemistry_snapshot_interval: int = Field(default=100, gt=0, description="Send snapshot to CPU for chemistry analysis every N steps (hybrid mode)")
    
    # Spatial hashing - DISABLED on Windows due to LLVM compilation errors
    use_spatial_hash: bool = Field(default=False, description="Use O(n) spatial hashing instead of O(n²) (disabled on Windows)")
    spatial_hash_cell_size: float = Field(default=10.0, gt=0, description="Cell size for spatial hashing")
    
    # Diagnostics
    enable_diagnostics: bool = Field(default=True)
    diagnostics_dir: str = Field(default="diagnostics")
    diagnostics_frequency: int = Field(default=10, gt=0, description="Log diagnostics every N steps")
    
    # Memory management (based on NAMD/GROMACS best practices)
    memory_check_interval: int = Field(default=1000, gt=0, description="Check memory usage every N steps")
    gc_interval: int = Field(default=5000, gt=0, description="Run Python garbage collection every N steps")
    
    # Physics Database (PHASE 1 WEEK 2: Literature-based parameters)
    use_physics_db: bool = Field(default=True, description="Use literature parameters from physics database")
    physics_db_path: str = Field(default="data/physics_parameters.json", description="Path to physics parameters database (relative to project root)")
    
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
    
    # Bond system configuration - MORE RESTRICTIVE for performance
    # Literature: Bond lifetimes 10^-9 s (vdW) to 10^6 s (covalent) at 298K
    # FIXED: More restrictive binding conditions to reduce computational load
    enable_spring_forces: bool = Field(default=True, description="Enable spring forces from bonds")
    enable_advanced_breaking: bool = Field(default=True, description="Enable overload/aging bond breaking")
    theta_bind: float = Field(default=0.4, gt=0)  # INCREASED from 0.3 - more restrictive bond formation
    theta_break: float = Field(default=1.8, gt=0)  # INCREASED from 1.5 - stronger bonds required
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
