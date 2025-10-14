"""
Phase 2 Configuration Extension
================================

Adds Phase 2-specific parameters to Live 2.0 simulation system.
Supports Miller-Urey, Hydrothermal, and Formamide scenarios.
"""

from typing import Dict, List, Optional, Any
from pydantic import BaseModel, Field


class InitialMolecule(BaseModel):
    """Configuration for initial molecule species"""
    name: str = Field(..., description="Molecule name (e.g., 'methane', 'water')")
    formula: str = Field(..., description="Chemical formula (e.g., 'CH4', 'H2O')")
    count: int = Field(..., gt=0, description="Number of molecules to add")
    
    # Optional: Pre-bonded structure
    atoms: Optional[List[str]] = Field(default=None, description="List of atom symbols")
    bonds: Optional[List[tuple]] = Field(default=None, description="List of (atom1_idx, atom2_idx, order)")


class EnergyInjectionConfig(BaseModel):
    """Configuration for energy injection (electrical discharge, UV, etc.)"""
    enabled: bool = Field(default=False)
    type: str = Field(default="none", pattern="^(none|electrical|electrical_discharge|uv|thermal)$")
    
    # Pulse parameters
    pulse_interval: int = Field(default=1000, gt=0, description="Steps between pulses")
    pulse_energy: float = Field(default=50.0, gt=0, description="Energy per pulse (kJ/mol)")
    pulse_duration: int = Field(default=10, gt=0, description="Pulse duration (steps)")
    
    # Spatial distribution
    pulse_radius: float = Field(default=10.0, gt=0, description="Pulse radius (Angstrom)")
    pulse_shape: str = Field(default="spherical", pattern="^(spherical|planar|uniform)$")


class CatalystConfig(BaseModel):
    """Configuration for catalytic surfaces/particles"""
    name: str = Field(..., description="Catalyst name (e.g., 'FeS', 'TiO2')")
    concentration: float = Field(..., ge=0, le=1, description="Fraction of particles as catalyst")
    effect: str = Field(..., description="Catalytic effect type")
    boost_factor: float = Field(default=1.0, gt=0, description="Reaction rate boost")


class Phase2Config(BaseModel):
    """Phase 2 prebiotic scenario configuration"""
    
    # Scenario identification
    scenario_name: str = Field(..., description="Scenario name (miller_urey, hydrothermal, formamide)")
    description: str = Field(default="", description="Scenario description")
    
    # Initial composition
    initial_molecules: List[InitialMolecule] = Field(default_factory=list)
    
    # Energy injection (electrical discharge, UV, etc.)
    energy_injection: EnergyInjectionConfig = Field(default_factory=EnergyInjectionConfig)
    
    # Catalysts
    catalysts: List[CatalystConfig] = Field(default_factory=list)
    
    # Environmental conditions
    temperature: float = Field(default=298.0, gt=0, description="Temperature (Kelvin)")
    temperature_control: bool = Field(default=True)
    temperature_range: Optional[tuple] = Field(default=None, description="(min, max) for fluctuations")
    
    pH: Optional[float] = Field(default=7.0, ge=0, le=14, description="pH value")
    
    # Phase 2 specific tracking
    track_known_molecules: bool = Field(default=True)
    expected_products: List[str] = Field(default_factory=list, description="Expected molecules to track")
    
    # Output
    output_base_dir: str = Field(default="results/phase2")
    save_snapshots: bool = Field(default=True)
    snapshot_interval: int = Field(default=50000, description="Steps between snapshots")
    
    # Simulation parameters (from YAML simulation section)
    n_particles: Optional[int] = Field(default=None, description="Max particles override")
    dt: Optional[float] = Field(default=None, description="Timestep override")
    box_size: Optional[float] = Field(default=None, description="Box size override")
    max_steps: Optional[int] = Field(default=None, description="Max steps override")
    
    # Performance parameters
    enable_thermodynamic_validation: Optional[bool] = Field(default=None)
    validate_every_n_steps: Optional[int] = Field(default=None)
    energy_update_interval: Optional[int] = Field(default=None)
    metrics_update_interval: Optional[int] = Field(default=None)
    enable_diagnostics: Optional[bool] = Field(default=None)
    diagnostics_frequency: Optional[int] = Field(default=None)


def load_phase2_config_from_yaml(yaml_path: str) -> Phase2Config:
    """
    Load Phase 2 configuration from YAML file
    
    Args:
        yaml_path: Path to YAML configuration file
    
    Returns:
        Phase2Config object
    """
    import yaml
    
    with open(yaml_path, 'r') as f:
        config_dict = yaml.safe_load(f)
    
    # Extract Phase 2 relevant sections
    sim_config = config_dict.get('simulation', {})
    initial_mols = config_dict.get('initial_molecules', [])
    energy_inj = config_dict.get('energy_injection', {})
    catalysts = config_dict.get('catalysis', {}).get('catalysts', [])
    chemistry = config_dict.get('chemistry', {})
    benchmark = config_dict.get('benchmark', {})
    output = config_dict.get('output', {})
    
    # Convert to Phase2Config
    phase2_config = Phase2Config(
        scenario_name=sim_config.get('name', 'unknown'),
        description=sim_config.get('description', ''),
        
        # Initial molecules
        initial_molecules=[
            InitialMolecule(
                name=mol.get('name', 'unknown'),
                formula=mol.get('formula', ''),
                count=mol.get('count', 0)
            )
            for mol in initial_mols
        ],
        
        # Energy injection
        energy_injection=EnergyInjectionConfig(
            enabled=energy_inj.get('enabled', False),
            type=energy_inj.get('type', 'none'),
            pulse_interval=energy_inj.get('pulse_interval', 1000),
            pulse_energy=energy_inj.get('pulse_energy', 50.0),
            pulse_duration=energy_inj.get('pulse_duration', 10),
            pulse_radius=energy_inj.get('pulse_radius', 10.0),
            pulse_shape=energy_inj.get('pulse_shape', 'spherical')
        ),
        
        # Catalysts
        catalysts=[
            CatalystConfig(
                name=cat.get('name', 'unknown'),
                concentration=cat.get('concentration', 0.01),
                effect=cat.get('effect', 'none'),
                boost_factor=cat.get('boost_factor', 1.0)
            )
            for cat in catalysts
        ],
        
        # Conditions
        temperature=sim_config.get('target_temperature', 298.0),
        temperature_control=sim_config.get('temperature_control', True),
        pH=chemistry.get('pH', 7.0),
        
        # Tracking
        track_known_molecules=benchmark.get('track_known_molecules', True),
        expected_products=benchmark.get('expected_products', []),
        
        # Output
        output_base_dir=output.get('base_dir', 'results/phase2'),
        save_snapshots=output.get('save_snapshots', True),
        snapshot_interval=output.get('snapshot_interval', 50000),
        
        # Simulation parameters from YAML
        n_particles=sim_config.get('n_particles'),
        dt=sim_config.get('dt'),
        box_size=sim_config.get('box_size'),
        max_steps=sim_config.get('max_steps'),
        
        # Performance parameters from YAML
        enable_thermodynamic_validation=sim_config.get('enable_thermodynamic_validation'),
        validate_every_n_steps=sim_config.get('validate_every_n_steps'),
        energy_update_interval=sim_config.get('energy_update_interval'),
        metrics_update_interval=sim_config.get('metrics_update_interval'),
        enable_diagnostics=sim_config.get('enable_diagnostics'),
        diagnostics_frequency=sim_config.get('diagnostics_frequency')
    )
    
    return phase2_config


def create_miller_urey_config() -> Phase2Config:
    """Create default Miller-Urey configuration"""
    return Phase2Config(
        scenario_name="miller_urey",
        description="Miller-Urey 1953 - Reducing atmosphere with electrical discharge",
        
        initial_molecules=[
            InitialMolecule(name="methane", formula="CH4", count=500),
            InitialMolecule(name="ammonia", formula="NH3", count=400),
            InitialMolecule(name="water", formula="H2O", count=800),
            InitialMolecule(name="hydrogen", formula="H2", count=300),
        ],
        
        energy_injection=EnergyInjectionConfig(
            enabled=True,
            type="electrical",
            pulse_interval=1000,
            pulse_energy=50.0,
            pulse_duration=10,
            pulse_radius=10.0,
            pulse_shape="spherical"
        ),
        
        temperature=298.0,
        pH=7.0,
        
        expected_products=["glycine", "alanine", "formaldehyde", "HCN", "formic_acid"]
    )


def create_hydrothermal_config() -> Phase2Config:
    """Create default hydrothermal vent configuration"""
    return Phase2Config(
        scenario_name="hydrothermal",
        description="Alkaline hydrothermal vent conditions",
        
        initial_molecules=[
            InitialMolecule(name="hydrogen", formula="H2", count=600),
            InitialMolecule(name="hydrogen_sulfide", formula="H2S", count=300),
            InitialMolecule(name="carbon_dioxide", formula="CO2", count=400),
            InitialMolecule(name="water", formula="H2O", count=600),
            InitialMolecule(name="ammonia", formula="NH3", count=100),
        ],
        
        catalysts=[
            CatalystConfig(name="FeS", concentration=0.01, effect="bond_formation_boost", boost_factor=2.0),
            CatalystConfig(name="FeS2", concentration=0.005, effect="electron_transfer", boost_factor=1.5),
        ],
        
        temperature=373.0,  # 100°C
        temperature_control=True,
        temperature_range=(323.0, 423.0),  # 50-150°C
        pH=10.0,  # Alkaline
        
        expected_products=["formic_acid", "acetic_acid", "pyruvic_acid", "thiols"]
    )


def create_formamide_config() -> Phase2Config:
    """Create default formamide-rich configuration"""
    return Phase2Config(
        scenario_name="formamide",
        description="Formamide-rich with UV radiation and mineral catalysts",
        
        initial_molecules=[
            InitialMolecule(name="formamide", formula="HCONH2", count=1200),
            InitialMolecule(name="water", formula="H2O", count=400),
            InitialMolecule(name="ammonia", formula="NH3", count=200),
            InitialMolecule(name="hydrogen_cyanide", formula="HCN", count=200),
        ],
        
        energy_injection=EnergyInjectionConfig(
            enabled=True,
            type="uv",
            pulse_interval=500,
            pulse_energy=20.0,
            pulse_duration=20,
            pulse_radius=15.0,
            pulse_shape="uniform"
        ),
        
        catalysts=[
            CatalystConfig(name="TiO2", concentration=0.02, effect="uv_enhancement", boost_factor=3.0),
            CatalystConfig(name="ZnS", concentration=0.01, effect="bond_formation_boost", boost_factor=2.0),
        ],
        
        temperature=323.0,  # 50°C
        pH=7.0,
        
        expected_products=["adenine", "guanine", "cytosine", "uracil", "glycine", "glycolaldehyde"]
    )


# Example usage
if __name__ == "__main__":
    # Example: Load from YAML
    config = load_phase2_config_from_yaml("configs/phase2_miller_urey.yaml")
    print(f"Loaded Phase 2 config: {config.scenario_name}")
    print(f"Initial molecules: {len(config.initial_molecules)}")
    print(f"Energy injection: {config.energy_injection.enabled}")
    
    # Example: Create programmatically
    miller_urey = create_miller_urey_config()
    print(f"\nMiller-Urey config: {miller_urey.description}")

