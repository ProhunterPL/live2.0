"""
Phase 2 Molecule Initializer
=============================

Initializes prebiotic molecules for Phase 2 simulations.
Converts molecular formulas to atoms and places them in the simulation.
"""

import taichi as ti
import numpy as np
import logging
from typing import List, Dict, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

from .phase2_config import Phase2Config, InitialMolecule
from .core.particles import ParticleSystem
from .config import SimulationConfig


# Atomic data
ATOM_DATA = {
    'H': {'mass': 1.008, 'type': 0, 'symbol': 'H'},
    'C': {'mass': 12.011, 'type': 1, 'symbol': 'C'},
    'N': {'mass': 14.007, 'type': 2, 'symbol': 'N'},
    'O': {'mass': 15.999, 'type': 3, 'symbol': 'O'},
    'S': {'mass': 32.06, 'type': 4, 'symbol': 'S'},
    'P': {'mass': 30.974, 'type': 5, 'symbol': 'P'},
}

# Molecular geometries (simplified 2D projections)
# Format: {formula: [(atom_symbol, (x_offset, y_offset)), ...]}
MOLECULAR_GEOMETRIES = {
    # Simple molecules
    'H2': [('H', (0.0, 0.0)), ('H', (0.74, 0.0))],
    'H2O': [('O', (0.0, 0.0)), ('H', (0.96, 0.0)), ('H', (0.48, 0.83))],
    'NH3': [('N', (0.0, 0.0)), ('H', (1.0, 0.0)), ('H', (-0.5, 0.87)), ('H', (-0.5, -0.87))],
    'CH4': [('C', (0.0, 0.0)), ('H', (1.09, 0.0)), ('H', (-1.09, 0.0)), ('H', (0.0, 1.09)), ('H', (0.0, -1.09))],
    'CO2': [('C', (0.0, 0.0)), ('O', (-1.16, 0.0)), ('O', (1.16, 0.0))],
    'H2S': [('S', (0.0, 0.0)), ('H', (1.34, 0.0)), ('H', (0.67, 1.16))],
    'HCN': [('H', (0.0, 0.0)), ('C', (1.07, 0.0)), ('N', (2.22, 0.0))],
    
    # Organics
    'HCONH2': [  # Formamide
        ('H', (0.0, 0.0)),
        ('C', (1.09, 0.0)),
        ('O', (1.59, 1.22)),
        ('N', (1.59, -1.22)),
        ('H', (2.59, -1.22)),
        ('H', (1.09, -2.12))
    ],
    'CH3OH': [  # Methanol
        ('C', (0.0, 0.0)),
        ('O', (1.43, 0.0)),
        ('H', (2.03, 0.0)),
        ('H', (-0.5, 0.87)),
        ('H', (-0.5, -0.87)),
        ('H', (-1.0, 0.0))
    ],
    'CH2O': [  # Formaldehyde
        ('C', (0.0, 0.0)),
        ('O', (1.22, 0.0)),
        ('H', (-0.6, 0.87)),
        ('H', (-0.6, -0.87))
    ],
}


class Phase2Initializer:
    """Initializes Phase 2 prebiotic scenarios"""
    
    def __init__(self, sim_config: SimulationConfig, phase2_config: Phase2Config):
        """
        Initialize Phase 2 molecule placer
        
        Args:
            sim_config: Main simulation configuration
            phase2_config: Phase 2 specific configuration
        """
        self.sim_config = sim_config
        self.phase2_config = phase2_config
        
        # Simulation box size
        self.box_width = float(sim_config.grid_width)
        self.box_height = float(sim_config.grid_height)
        
        # Placement grid to avoid overlaps
        self.grid_spacing = 5.0  # Angstroms between molecules
        
        logger.info(f"Phase 2 Initializer created for scenario: {phase2_config.scenario_name}")
        logger.info(f"Box size: {self.box_width} x {self.box_height}")
    
    def initialize_molecules(self, particle_system: ParticleSystem):
        """
        Initialize all molecules from Phase 2 config
        
        Args:
            particle_system: ParticleSystem to add molecules to
        """
        logger.info("=" * 70)
        logger.info("INITIALIZING PHASE 2 MOLECULES")
        logger.info("=" * 70)
        
        total_molecules = 0
        total_atoms = 0
        
        for mol_config in self.phase2_config.initial_molecules:
            n_placed, n_atoms = self._place_molecules(particle_system, mol_config)
            total_molecules += n_placed
            total_atoms += n_atoms
            
            logger.info(
                f"  {mol_config.name} ({mol_config.formula}): "
                f"{n_placed}/{mol_config.count} molecules, {n_atoms} atoms"
            )
        
        logger.info("=" * 70)
        logger.info(f"TOTAL: {total_molecules} molecules, {total_atoms} atoms")
        logger.info("=" * 70)
        
        return total_molecules, total_atoms
    
    def _place_molecules(self, particle_system: ParticleSystem, 
                        mol_config: InitialMolecule) -> Tuple[int, int]:
        """
        Place multiple copies of a molecule
        
        Returns:
            (molecules_placed, atoms_placed)
        """
        formula = mol_config.formula
        count = mol_config.count
        
        # Get molecular geometry
        if formula not in MOLECULAR_GEOMETRIES:
            logger.warning(f"No geometry for {formula}, using fallback")
            geometry = self._fallback_geometry(formula)
        else:
            geometry = MOLECULAR_GEOMETRIES[formula]
        
        molecules_placed = 0
        atoms_placed = 0
        
        # Place molecules in grid pattern to avoid overlaps
        for i in range(count):
            # Random position in box
            center_x = np.random.uniform(10, self.box_width - 10)
            center_y = np.random.uniform(10, self.box_height - 10)
            
            # Random rotation
            angle = np.random.uniform(0, 2 * np.pi)
            cos_a = np.cos(angle)
            sin_a = np.sin(angle)
            
            # Place all atoms in this molecule
            molecule_atoms = []
            for atom_symbol, (dx, dy) in geometry:
                # Rotate and translate
                x = center_x + dx * cos_a - dy * sin_a
                y = center_y + dx * sin_a + dy * cos_a
                
                # Ensure in box
                x = np.clip(x, 2.0, self.box_width - 2.0)
                y = np.clip(y, 2.0, self.box_height - 2.0)
                
                # Get atom data
                atom_data = ATOM_DATA.get(atom_symbol, ATOM_DATA['C'])
                
                # Create particle
                pos = ti.Vector([x, y])
                
                # Initial velocity (Maxwell-Boltzmann at temperature)
                temp = self.phase2_config.temperature
                vel_scale = np.sqrt(temp / 298.0) * 0.5  # Scale by sqrt(T)
                vx = np.random.normal(0, vel_scale)
                vy = np.random.normal(0, vel_scale)
                vel = ti.Vector([vx, vy])
                
                # Attributes: [mass, charge_x, charge_y, charge_z]
                attributes = ti.Vector([
                    atom_data['mass'],
                    0.0,  # charge - would need proper calculation
                    0.0,
                    0.0
                ])
                
                # Add particle
                try:
                    particle_system.add_particle_py(
                        pos=pos,
                        vel=vel,
                        attributes=attributes,
                        type_id=atom_data['type'],
                        binding_sites=4,  # Allow bonding
                        binding_strength=1.0
                    )
                    molecule_atoms.append(particle_system.particle_count[None] - 1)
                    atoms_placed += 1
                except Exception as e:
                    logger.error(f"Failed to add atom {atom_symbol}: {e}")
                    break
            
            if len(molecule_atoms) == len(geometry):
                molecules_placed += 1
                # TODO: Create bonds between atoms in molecule
                # Would need to call binding system here
        
        return molecules_placed, atoms_placed
    
    def _fallback_geometry(self, formula: str) -> List[Tuple[str, Tuple[float, float]]]:
        """
        Create fallback geometry for unknown formulas
        Parses simple formulas like "C2H6"
        """
        # Simple parser
        atoms = []
        i = 0
        while i < len(formula):
            if formula[i].isupper():
                symbol = formula[i]
                i += 1
                
                # Check for lowercase (e.g., Cl)
                if i < len(formula) and formula[i].islower():
                    symbol += formula[i]
                    i += 1
                
                # Check for count
                count_str = ''
                while i < len(formula) and formula[i].isdigit():
                    count_str += formula[i]
                    i += 1
                
                count = int(count_str) if count_str else 1
                atoms.extend([symbol] * count)
            else:
                i += 1
        
        # Create linear geometry
        geometry = []
        spacing = 1.5
        for idx, symbol in enumerate(atoms):
            x = idx * spacing
            y = 0.0
            geometry.append((symbol, (x, y)))
        
        return geometry
    
    def setup_energy_injection(self, stepper):
        """
        Configure energy injection for the simulation
        
        Args:
            stepper: SimulationStepper instance
        """
        if not self.phase2_config.energy_injection.enabled:
            logger.info("Energy injection disabled")
            return
        
        inj_config = self.phase2_config.energy_injection
        
        logger.info(f"Energy injection enabled:")
        logger.info(f"  Type: {inj_config.type}")
        logger.info(f"  Interval: {inj_config.pulse_interval} steps")
        logger.info(f"  Energy: {inj_config.pulse_energy} kJ/mol")
        logger.info(f"  Radius: {inj_config.pulse_radius} Angstrom")
        
        # TODO: Configure stepper to apply energy pulses
        # This would modify the energy_manager or grid to inject energy
        # at specified intervals
        
        # For now, adjust pulse parameters in config
        stepper.config.pulse_every = inj_config.pulse_interval
        stepper.config.pulse_amplitude = inj_config.pulse_energy / 50.0  # Scale
        stepper.config.pulse_radius = inj_config.pulse_radius
    
    def setup_catalysts(self, particle_system: ParticleSystem):
        """
        Add catalyst particles to the simulation
        
        Args:
            particle_system: ParticleSystem to add catalysts to
        """
        if not self.phase2_config.catalysts:
            logger.info("No catalysts specified")
            return
        
        logger.info("Adding catalysts:")
        
        for catalyst in self.phase2_config.catalysts:
            # Number of catalyst particles
            n_catalyst = int(self.sim_config.max_particles * catalyst.concentration)
            
            logger.info(
                f"  {catalyst.name}: {n_catalyst} particles "
                f"({catalyst.concentration*100:.1f}%)"
            )
            
            # Place catalyst particles
            for i in range(n_catalyst):
                x = np.random.uniform(0, self.box_width)
                y = np.random.uniform(0, self.box_height)
                
                pos = ti.Vector([x, y])
                vel = ti.Vector([0.0, 0.0])  # Stationary catalysts
                
                # Heavy, neutral particles
                attributes = ti.Vector([100.0, 0.0, 0.0, 0.0])
                
                # Special type for catalyst
                catalyst_type = 10 + len(self.phase2_config.catalysts)  # Unique type
                
                try:
                    particle_system.add_particle_py(
                        pos=pos,
                        vel=vel,
                        attributes=attributes,
                        type_id=catalyst_type,
                        binding_sites=8,  # High binding sites
                        binding_strength=catalyst.boost_factor
                    )
                except Exception as e:
                    logger.error(f"Failed to add catalyst: {e}")
                    break


# Convenience function
def initialize_phase2_simulation(sim_config: SimulationConfig,
                                 phase2_config: Phase2Config,
                                 particle_system: ParticleSystem,
                                 stepper) -> Dict:
    """
    Complete Phase 2 initialization
    
    Args:
        sim_config: Main simulation configuration
        phase2_config: Phase 2 specific configuration
        particle_system: ParticleSystem to initialize
        stepper: SimulationStepper instance
    
    Returns:
        Initialization statistics
    """
    initializer = Phase2Initializer(sim_config, phase2_config)
    
    # Add molecules
    n_molecules, n_atoms = initializer.initialize_molecules(particle_system)
    
    # Setup energy injection
    initializer.setup_energy_injection(stepper)
    
    # Add catalysts
    initializer.setup_catalysts(particle_system)
    
    stats = {
        'scenario': phase2_config.scenario_name,
        'molecules_placed': n_molecules,
        'atoms_placed': n_atoms,
        'catalysts': len(phase2_config.catalysts),
        'energy_injection_enabled': phase2_config.energy_injection.enabled,
        'temperature': phase2_config.temperature
    }
    
    logger.info("Phase 2 initialization complete!")
    logger.info(f"Statistics: {stats}")
    
    return stats

