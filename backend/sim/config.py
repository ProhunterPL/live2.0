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
    
    # Time settings
    dt: float = Field(default=0.01, gt=0, le=1.0)  # Bardzo mały krok dla stabilności numerycznej
    max_time: float = Field(default=1000.0, gt=0)
    
    # Energy settings
    energy_decay: float = Field(default=0.96, gt=0, lt=1)  # Wolniejsze wygaszanie z 0.95
    energy_threshold: float = Field(default=0.1, gt=0)
    pulse_every: int = Field(default=24, gt=0)  # Częstsze pulsy (co 24 kroki)
    pulse_radius: float = Field(default=32.0, gt=0)  # Większy promień (32 jednostki)
    pulse_amplitude: float = Field(default=8.0, gt=0)  # Większa amplituda (8.0)
    diffuse_D: float = Field(default=0.5, gt=0)  # Szybsza dyfuzja energii
    target_energy: float = Field(default=0.5, gt=0)  # Wyższe tło energii
    thermostat_alpha: float = Field(default=0.01, gt=0)  # Silniejszy termostat
    
    # Particle settings
    max_particles: int = Field(default=10000, gt=0, le=100000)
    particle_radius: float = Field(default=0.5, gt=0)
    
    # Binding settings - OPTIMIZED for active bond formation
    binding_threshold: float = Field(default=0.25, gt=0, le=1)  # NAPRAWIONE: Obniżone z 0.5 dla łatwiejszego tworzenia wiązań
    unbinding_threshold: float = Field(default=0.15, gt=0, le=1)  # NAPRAWIONE: Zwiększone z 0.1 dla stabilnych wiązań
    
    # Novelty detection
    novelty_window: int = Field(default=100, gt=0)
    min_cluster_size: int = Field(default=1, ge=1)  # Obniżone z 2 dla wykrywania pojedynczych wiązań
    
    # Visualization
    vis_frequency: int = Field(default=5, gt=0)  # Zmniejszone z 3 dla stabilności
    log_frequency: int = Field(default=100, gt=0)
    
    # Mutations
    p_mut_base: float = Field(default=5e-4, gt=0)  # Zwiększone mutacje (5x więcej)
    p_mut_gain: float = Field(default=20.0, gt=0)  # Zwiększone z 14.0
    attr_sigma: float = Field(default=0.12, gt=0)  # Silniejsza perturbacja
    
    # Performance optimization parameters
    energy_update_interval: int = Field(default=5, gt=0, description="Update energy every N steps")
    metrics_update_interval: int = Field(default=1, gt=0, description="Update metrics every N steps")
    
    # Diagnostics
    enable_diagnostics: bool = Field(default=True)
    diagnostics_dir: str = Field(default="diagnostics")
    diagnostics_frequency: int = Field(default=10, gt=0, description="Log diagnostics every N steps")
    
    # Random seed
    seed: Optional[int] = Field(default=None)

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
    
    # Bond system configuration
    enable_spring_forces: bool = Field(default=True, description="Enable spring forces from bonds")
    enable_advanced_breaking: bool = Field(default=True, description="Enable overload/aging bond breaking")
    theta_bind: float = Field(default=0.25, gt=0)  # NAPRAWIONE: obniżone z 0.3 dla łatwiejszego tworzenia wiązań
    theta_break: float = Field(default=1.5, gt=0)  # Podniesione z 1.15 dla stabilniejszych wiązań
    vmax: float = Field(default=8.0, gt=0)  # NOWE: zwiększone z 4.5
    neighbor_radius: float = Field(default=3.2, gt=0)  # NOWE
    rebuild_neighbors_every: int = Field(default=8, gt=0)  # NOWE
    clamp_density_per_cell: int = Field(default=64, gt=0)  # NOWE
    bond_types: Dict[int, Dict[str, float]] = Field(default={
        0: {'name': 'vdW', 'k_spring': 2.0, 'rest_len_factor': 1.0, 'damping': 0.1, 'strength': 5.0, 'max_age': 50.0},
        1: {'name': 'covalent', 'k_spring': 10.0, 'rest_len_factor': 0.9, 'damping': 0.2, 'strength': 20.0, 'max_age': 200.0},
        2: {'name': 'H-bond', 'k_spring': 5.0, 'rest_len_factor': 1.1, 'damping': 0.15, 'strength': 10.0, 'max_age': 80.0},
        3: {'name': 'metallic', 'k_spring': 7.0, 'rest_len_factor': 1.0, 'damping': 0.25, 'strength': 15.0, 'max_age': 150.0}
    }, description="Bond type parameters")
    
    # Cluster configuration
    min_cluster_size: int = Field(default=4, ge=2, description="Minimum cluster size to track")
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
