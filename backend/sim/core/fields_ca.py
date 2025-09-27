"""
Cellular Automata fields for Live 2.0 Preset Prebiotic mode
Handles concentration fields, diffusion, and chemical reactions
"""

import taichi as ti
import numpy as np
from typing import Dict, List, Tuple, Optional
from ..config import PresetPrebioticConfig

@ti.data_oriented
class ConcentrationFields:
    """Manages concentration fields for chemical species"""
    
    def __init__(self, config: PresetPrebioticConfig, width: int, height: int):
        self.config = config
        self.width = width
        self.height = height
        
        # Species mapping
        self.species_names = list(config.species.keys())
        self.num_species = len(self.species_names)
        self.species_to_index = {name: i for i, name in enumerate(self.species_names)}
        
        # Concentration fields
        self.concentrations = ti.field(dtype=ti.f32, 
                                     shape=(self.num_species, self.width, self.height))
        
        # Diffusion coefficients
        self.diffusion_coeffs = ti.field(dtype=ti.f32, shape=(self.num_species,))
        
        # Reaction rates
        self.reaction_rates = ti.field(dtype=ti.f32, shape=(self.num_species, self.num_species))
        
        # Initialize fields
        self.initialize_fields()
    
    def initialize_fields(self):
        """Initialize concentration fields with initial values"""
        # Set diffusion coefficients
        for i, species_name in enumerate(self.species_names):
            if species_name in self.config.diffusion_coeffs:
                self.diffusion_coeffs[i] = self.config.diffusion_coeffs[species_name]
            else:
                self.diffusion_coeffs[i] = 0.1  # Default diffusion coefficient
        
        # Set reaction rates
        for i in range(self.num_species):
            for j in range(self.num_species):
                self.reaction_rates[i, j] = 0.0
        
        # Set specific reaction rates
        for reaction_name, rate in self.config.reaction_rates.items():
            if reaction_name == "HCN_to_NH2CHO":
                hcn_idx = self.species_to_index.get("HCN", -1)
                nh2cho_idx = self.species_to_index.get("NH2CHO", -1)
                if hcn_idx >= 0 and nh2cho_idx >= 0:
                    self.reaction_rates[hcn_idx, nh2cho_idx] = rate
        
        # Initialize concentrations
        self.set_initial_concentrations()
    
    @ti.kernel
    def set_initial_concentrations(self):
        """Set initial concentration values"""
        for i, j in ti.ndrange(self.width, self.height):
            for s in range(self.num_species):
                species_name = self.species_names[s]
                initial_conc = self.config.species[species_name]
                self.concentrations[s, i, j] = initial_conc
    
    @ti.kernel
    def apply_diffusion(self, dt: ti.f32):
        """Apply diffusion using discrete Laplacian"""
        # Create temporary field
        temp_field = ti.field(dtype=ti.f32, 
                            shape=(self.num_species, self.width, self.height))
        
        # Copy current concentrations
        for s, i, j in ti.ndrange(self.num_species, self.width, self.height):
            temp_field[s, i, j] = self.concentrations[s, i, j]
        
        # Apply diffusion
        for s, i, j in ti.ndrange(self.num_species, self.width, self.height):
            diffusion_coeff = self.diffusion_coeffs[s]
            
            # Get neighbors with periodic boundary conditions
            left = (i - 1) % self.width
            right = (i + 1) % self.width
            up = (j - 1) % self.height
            down = (j + 1) % self.height
            
            # Discrete Laplacian
            laplacian = (temp_field[s, left, j] + temp_field[s, right, j] + 
                        temp_field[s, i, up] + temp_field[s, i, down] - 
                        4.0 * temp_field[s, i, j])
            
            # Update concentration
            self.concentrations[s, i, j] += diffusion_coeff * laplacian * dt
    
    @ti.kernel
    def apply_reactions(self, dt: ti.f32, energy_field: ti.template()):
        """Apply chemical reactions"""
        for i, j in ti.ndrange(self.width, self.height):
            energy = energy_field[i, j]
            energy_factor = 1.0 + energy  # Higher energy increases reaction rates
            
            # Apply reactions
            for s1 in range(self.num_species):
                for s2 in range(self.num_species):
                    if s1 != s2 and self.reaction_rates[s1, s2] > 0:
                        # Reaction: s1 -> s2
                        rate = self.reaction_rates[s1, s2] * energy_factor
                        
                        # First-order reaction
                        conc_s1 = self.concentrations[s1, i, j]
                        reaction_rate = rate * conc_s1
                        
                        # Update concentrations
                        self.concentrations[s1, i, j] -= reaction_rate * dt
                        self.concentrations[s2, i, j] += reaction_rate * dt
    
    @ti.kernel
    def add_source(self, species_idx: ti.i32, pos: ti.template(), 
                  intensity: ti.f32, radius: ti.f32):
        """Add concentration source"""
        for i, j in ti.ndrange(self.width, self.height):
            dx = i - pos[0]
            dy = j - pos[1]
            
            # Handle periodic boundary conditions
            if dx > self.width / 2:
                dx -= self.width
            elif dx < -self.width / 2:
                dx += self.width
            
            if dy > self.height / 2:
                dy -= self.height
            elif dy < -self.height / 2:
                dy += self.height
            
            distance = ti.sqrt(dx * dx + dy * dy)
            
            if distance <= radius:
                energy_factor = 1.0 - (distance / radius)
                energy_factor = ti.max(energy_factor, 0.0)
                self.concentrations[species_idx, i, j] += intensity * energy_factor
    
    @ti.kernel
    def apply_boundary_conditions(self):
        """Apply boundary conditions (periodic)"""
        # Concentrations are already periodic due to modulo operations
        pass
    
    def get_concentration(self, species_name: str) -> np.ndarray:
        """Get concentration field for a specific species"""
        if species_name not in self.species_to_index:
            raise ValueError(f"Unknown species: {species_name}")
        
        species_idx = self.species_to_index[species_name]
        return self.concentrations[species_idx, :, :].to_numpy()
    
    def get_all_concentrations(self) -> Dict[str, np.ndarray]:
        """Get all concentration fields"""
        concentrations = {}
        for species_name in self.species_names:
            concentrations[species_name] = self.get_concentration(species_name)
        return concentrations
    
    def set_concentration(self, species_name: str, field: np.ndarray):
        """Set concentration field for a specific species"""
        if species_name not in self.species_to_index:
            raise ValueError(f"Unknown species: {species_name}")
        
        species_idx = self.species_to_index[species_name]
        self.concentrations[species_idx, :, :].from_numpy(field)
    
    def get_total_mass(self) -> Dict[str, float]:
        """Get total mass of each species"""
        masses = {}
        for species_name in self.species_names:
            field = self.get_concentration(species_name)
            masses[species_name] = float(np.sum(field))
        return masses
    
    def get_stats(self) -> Dict:
        """Get field statistics"""
        stats = {}
        for species_name in self.species_names:
            field = self.get_concentration(species_name)
            stats[species_name] = {
                'total_mass': float(np.sum(field)),
                'max_concentration': float(np.max(field)),
                'mean_concentration': float(np.mean(field)),
                'std_concentration': float(np.std(field))
            }
        return stats

@ti.data_oriented
class ReactionSystem:
    """Manages chemical reactions in preset mode"""
    
    def __init__(self, config: PresetPrebioticConfig):
        self.config = config
        
        # Reaction definitions
        self.reactions = []
        self.reaction_rates = {}
        
        # Initialize reactions
        self.initialize_reactions()
    
    def initialize_reactions(self):
        """Initialize reaction system"""
        # Define reactions
        self.reactions = [
            {
                'name': 'HCN_to_NH2CHO',
                'reactants': ['HCN'],
                'products': ['NH2CHO'],
                'rate': self.config.reaction_rates.get('HCN_to_NH2CHO', 0.01),
                'stoichiometry': {'HCN': -1, 'NH2CHO': 1}
            }
        ]
        
        # Set reaction rates
        for reaction in self.reactions:
            self.reaction_rates[reaction['name']] = reaction['rate']
    
    @ti.kernel
    def compute_reaction_rates(self, concentrations: ti.template(), 
                             energy_field: ti.template(), 
                             width: ti.i32, height: ti.i32,
                             reaction_rates: ti.template()) -> ti.f32:
        """Compute reaction rates at each grid point"""
        total_rate = 0.0
        
        for i, j in ti.ndrange(width, height):
            energy = energy_field[i, j]
            energy_factor = 1.0 + energy
            
            # HCN -> NH2CHO reaction
            hcn_conc = concentrations[0, i, j]  # Assuming HCN is index 0
            rate = reaction_rates[0] * energy_factor * hcn_conc
            total_rate += rate
        
        return total_rate
    
    def get_reaction_info(self) -> List[Dict]:
        """Get information about all reactions"""
        return self.reactions.copy()
    
    def get_reaction_stats(self) -> Dict:
        """Get reaction statistics"""
        return {
            'num_reactions': len(self.reactions),
            'reaction_rates': self.reaction_rates.copy()
        }

class PresetPrebioticSimulator:
    """Main simulator for preset prebiotic mode"""
    
    def __init__(self, config: PresetPrebioticConfig, width: int, height: int):
        self.config = config
        self.width = width
        self.height = height
        
        # Initialize components
        self.concentration_fields = ConcentrationFields(config, width, height)
        self.reaction_system = ReactionSystem(config)
        
        # Energy field (shared with main simulation)
        self.energy_field = ti.field(dtype=ti.f32, shape=(width, height))
        
        # Time tracking
        self.current_time = 0.0
        self.dt = 0.01
        
        # Initialize
        self.reset()
    
    def reset(self):
        """Reset simulator to initial state"""
        self.current_time = 0.0
        self.concentration_fields.set_initial_concentrations()
        
        # Clear energy field
        for i, j in ti.ndrange(self.width, self.height):
            self.energy_field[i, j] = 0.0
    
    def step(self, dt: float = None):
        """Perform one simulation step"""
        if dt is None:
            dt = self.dt
        
        # Apply diffusion
        self.concentration_fields.apply_diffusion(dt)
        
        # Apply reactions
        self.concentration_fields.apply_reactions(dt, self.energy_field)
        
        # Apply boundary conditions
        self.concentration_fields.apply_boundary_conditions()
        
        # Update time
        self.current_time += dt
    
    def add_energy_source(self, position: Tuple[float, float], intensity: float, radius: float):
        """Add energy source to the field"""
        pos_ti = ti.Vector(position)
        
        for i, j in ti.ndrange(self.width, self.height):
            dx = i - pos_ti[0]
            dy = j - pos_ti[1]
            
            # Handle periodic boundary conditions
            if dx > self.width / 2:
                dx -= self.width
            elif dx < -self.width / 2:
                dx += self.width
            
            if dy > self.height / 2:
                dy -= self.height
            elif dy < -self.height / 2:
                dy += self.height
            
            distance = ti.sqrt(dx * dx + dy * dy)
            
            if distance <= radius:
                energy_factor = 1.0 - (distance / radius)
                energy_factor = ti.max(energy_factor, 0.0)
                self.energy_field[i, j] += intensity * energy_factor
    
    def decay_energy_field(self, decay_rate: float):
        """Apply energy decay"""
        for i, j in ti.ndrange(self.width, self.height):
            self.energy_field[i, j] *= decay_rate
    
    def get_visualization_data(self) -> Dict:
        """Get data for visualization"""
        concentrations = self.concentration_fields.get_all_concentrations()
        energy_field = self.energy_field.to_numpy()
        
        return {
            'concentrations': {k: v.tolist() for k, v in concentrations.items()},
            'energy_field': energy_field.tolist(),
            'time': self.current_time,
            'stats': self.concentration_fields.get_stats()
        }
    
    def get_stats(self) -> Dict:
        """Get simulation statistics"""
        return {
            'current_time': self.current_time,
            'concentration_stats': self.concentration_fields.get_stats(),
            'reaction_stats': self.reaction_system.get_reaction_stats(),
            'total_energy': float(np.sum(self.energy_field.to_numpy()))
        }
