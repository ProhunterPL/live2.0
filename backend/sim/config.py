# Live 2.0 Configuration
from typing import Dict, Any, Optional
from pydantic import BaseModel, Field
import numpy as np

class SimulationConfig(BaseModel):
    """Main configuration for Live 2.0 simulation"""
    
    # Grid settings
    grid_height: int = Field(default=256, ge=64, le=1024)
    grid_width: int = Field(default=256, ge=64, le=1024)
    
    # Simulation mode
    mode: str = Field(default="open_chemistry", pattern="^(preset_prebiotic|open_chemistry)$")
    
    # Time settings
    dt: float = Field(default=0.01, gt=0, le=1.0)
    max_time: float = Field(default=1000.0, gt=0)
    
    # Energy settings
    energy_decay: float = Field(default=0.95, gt=0, lt=1)
    energy_threshold: float = Field(default=0.1, gt=0)
    
    # Particle settings
    max_particles: int = Field(default=10000, gt=0, le=100000)
    particle_radius: float = Field(default=0.5, gt=0)
    
    # Binding settings
    binding_threshold: float = Field(default=0.8, gt=0, le=1)
    unbinding_threshold: float = Field(default=0.2, gt=0, le=1)
    
    # Novelty detection
    novelty_window: int = Field(default=100, gt=0)
    min_cluster_size: int = Field(default=2, ge=1)
    
    # Visualization
    vis_frequency: int = Field(default=10, gt=0)
    log_frequency: int = Field(default=100, gt=0)
    
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

def get_default_config() -> SimulationConfig:
    """Get default simulation configuration"""
    return SimulationConfig()

def get_preset_config() -> PresetPrebioticConfig:
    """Get default preset prebiotic configuration"""
    return PresetPrebioticConfig()

def get_open_chemistry_config() -> OpenChemistryConfig:
    """Get default open chemistry configuration"""
    return OpenChemistryConfig()
