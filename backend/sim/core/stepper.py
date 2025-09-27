"""
Main simulation stepper for Live 2.0
Coordinates all simulation components and manages time evolution
"""

import taichi as ti
import numpy as np
import time
from typing import Dict, List, Optional, Any
from .config import SimulationConfig, PresetPrebioticConfig, OpenChemistryConfig
from .grid import Grid
from .particles import ParticleSystem
from .potentials import PotentialSystem
from .binding import BindingSystem
from .graphs import GraphProcessor, MolecularGraph
from .catalog import SubstanceCatalog
from .metrics import MetricsCollector, NoveltyTracker, MetricsAggregator
from .energy import EnergyManager
from .rng import RNG

class SimulationStepper:
    """Main simulation coordinator"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.current_time = 0.0
        self.step_count = 0
        
        # Initialize components
        self.grid = Grid(config)
        self.particles = ParticleSystem(config)
        self.potentials = PotentialSystem(config)
        self.binding = BindingSystem(config)
        self.graphs = GraphProcessor(config.max_particles)
        self.catalog = SubstanceCatalog()
        self.metrics = MetricsCollector(config.max_particles)
        self.novelty_tracker = NoveltyTracker()
        self.energy_manager = EnergyManager(config)
        self.rng = RNG(config.seed)
        
        # Metrics aggregator
        self.aggregator = MetricsAggregator()
        self.aggregator.set_components(self.metrics, self.novelty_tracker, self.catalog)
        
        # Mode-specific configurations
        self.preset_config = PresetPrebioticConfig()
        self.open_chemistry_config = OpenChemistryConfig()
        
        # Simulation state
        self.is_running = False
        self.is_paused = False
        
        # Initialize simulation
        self.initialize_simulation()
    
    def initialize_simulation(self):
        """Initialize simulation with initial conditions"""
        if self.config.mode == "preset_prebiotic":
            self.initialize_preset_mode()
        else:
            self.initialize_open_chemistry_mode()
        
        # Initialize metrics
        self.metrics.record_metrics()
        self.aggregator.update_aggregated_stats()
    
    def initialize_preset_mode(self):
        """Initialize preset prebiotic mode"""
        # Add initial particles based on preset configuration
        species = self.preset_config.species
        
        for species_name, concentration in species.items():
            if concentration > 0:
                # Register particle type
                type_id = self.particles.register_particle_type(
                    name=species_name,
                    mass=1.0,
                    charge=(0.0, 0.0, 0.0)
                )
                
                # Add particles randomly distributed
                num_particles = int(concentration * self.config.max_particles)
                for _ in range(num_particles):
                    pos = self.rng.next_vector2(0, self.config.grid_width)
                    vel = self.rng.next_vector2(-1, 1)
                    attributes = ti.Vector([1.0, 0.0, 0.0, 0.0])  # mass, charge_x, charge_y, charge_z
                    
                    self.particles.add_particle(pos, vel, attributes, type_id, 2, 1.0)
        
        # Add energy sources
        for _ in range(self.open_chemistry_config.energy_sources):
            pos = self.rng.next_vector2(0, self.config.grid_width)
            intensity = self.open_chemistry_config.energy_intensity
            radius = 10.0
            
            self.energy_manager.add_energy_source(
                pos.to_numpy(), intensity, radius, duration=100.0
            )
    
    def initialize_open_chemistry_mode(self):
        """Initialize open chemistry mode"""
        # Add initial particles with random properties
        num_initial_particles = min(100, self.config.max_particles)
        
        for i in range(num_initial_particles):
            # Random position
            pos = self.rng.next_vector2(0, self.config.grid_width)
            
            # Random velocity
            vel = self.rng.next_vector2(-0.5, 0.5)
            
            # Random attributes
            mass = self.rng.next_range(0.5, 2.0)
            charge_x = self.rng.next_range(-1.0, 1.0)
            charge_y = self.rng.next_range(-1.0, 1.0)
            charge_z = self.rng.next_range(-1.0, 1.0)
            attributes = ti.Vector([mass, charge_x, charge_y, charge_z])
            
            # Random type
            type_id = self.particles.register_particle_type(
                name=f"particle_{i}",
                mass=mass,
                charge=(charge_x, charge_y, charge_z)
            )
            
            # Add particle
            self.particles.add_particle(pos, vel, attributes, type_id, 2, 1.0)
        
        # Add energy sources
        for _ in range(self.open_chemistry_config.energy_sources):
            pos = self.rng.next_vector2(0, self.config.grid_width)
            intensity = self.open_chemistry_config.energy_intensity
            radius = 15.0
            
            self.energy_manager.add_energy_source(
                pos.to_numpy(), intensity, radius, duration=200.0
            )
    
    def step(self, dt: float = None):
        """Perform one simulation step"""
        if dt is None:
            dt = self.config.dt
        
        if not self.is_running or self.is_paused:
            return
        
        # Update energy system
        self.energy_manager.update(dt)
        
        # Update particle positions
        self.particles.update_positions(dt)
        
        # Update spatial hash
        self.grid.update_spatial_hash()
        
        # Compute forces
        self.potentials.compute_forces(
            self.particles.positions,
            self.particles.attributes,
            self.particles.active,
            self.particles.particle_count[None]
        )
        
        # Apply forces
        self.particles.apply_forces(self.potentials.forces, dt)
        
        # Update binding system
        self.binding.update_bonds(
            self.particles.positions,
            self.particles.attributes,
            self.particles.active,
            self.particles.particle_count[None],
            dt
        )
        
        # Update clusters
        self.binding.update_clusters(
            self.particles.active,
            self.particles.particle_count[None]
        )
        
        # Apply periodic boundary conditions
        self.grid.apply_periodic_boundary()
        
        # Update energy field
        self.grid.decay_energy_field(self.config.energy_decay)
        
        # Add energy from sources to particles
        self.add_energy_to_particles()
        
        # Apply mutations in high energy regions
        self.apply_mutations(dt)
        
        # Update graph representation
        self.update_graph_representation()
        
        # Detect novel substances
        self.detect_novel_substances()
        
        # Update metrics
        self.update_metrics()
        
        # Update time and step count
        self.current_time += dt
        self.step_count += 1
    
    def add_energy_to_particles(self):
        """Add energy from field to particles"""
        energy_amounts = ti.field(dtype=ti.f32, shape=(self.config.max_particles,))
        
        # Clear energy amounts
        for i in range(self.config.max_particles):
            energy_amounts[i] = 0.0
        
        # Add energy from field
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                pos = self.particles.positions[i]
                energy = self.energy_manager.get_energy_at_position(pos.to_numpy())
                energy_amounts[i] = energy
        
        # Apply energy to particles
        self.particles.add_energy(energy_amounts)
    
    def apply_mutations(self, dt: float):
        """Apply mutations to particles in high energy regions"""
        mutation_rate = self.open_chemistry_config.mutation_rate
        mutation_strength = self.open_chemistry_config.mutation_strength
        
        # Find high energy regions
        high_energy_regions = self.energy_manager.get_high_energy_regions()
        
        for region_pos in high_energy_regions:
            # Find particles near this region
            region_ti = ti.Vector(region_pos)
            neighbors = ti.field(dtype=ti.i32, shape=(100,))
            count = self.grid.get_neighbors(region_ti, 5.0, neighbors)
            
            # Apply mutations to nearby particles
            for i in range(count):
                particle_idx = neighbors[i]
                if self.particles.active[particle_idx] == 1:
                    # Check mutation probability
                    if self.rng.next() < mutation_rate * dt:
                        self.particles.mutate_particle(
                            particle_idx, mutation_strength, self.current_time, self.rng
                        )
    
    def update_graph_representation(self):
        """Update graph representation of molecular structures"""
        # Update adjacency matrix from bonds
        self.graphs.update_adjacency_matrix(
            self.binding.bond_matrix,
            self.particles.active,
            self.particles.particle_count[None]
        )
        
        # Compute degrees
        self.graphs.compute_degrees(
            self.particles.active,
            self.particles.particle_count[None]
        )
        
        # Compute graph hash
        self.graphs.compute_graph_hash(
            self.particles.active,
            self.particles.particle_count[None]
        )
    
    def detect_novel_substances(self):
        """Detect and catalog novel molecular structures"""
        # Get clusters
        clusters = self.binding.get_clusters(min_size=self.config.min_cluster_size)
        
        for cluster in clusters:
            if len(cluster) >= self.config.min_cluster_size:
                # Create molecular graph
                bonds = self.binding.get_bonds()
                cluster_bonds = [(i, j) for i, j, _ in bonds if i in cluster and j in cluster]
                
                # Get particle attributes for cluster
                particle_attributes = {}
                for particle_idx in cluster:
                    if self.particles.active[particle_idx] == 1:
                        particle_attributes[particle_idx] = self.particles.attributes[particle_idx].to_numpy()
                
                # Create graph
                graph = MolecularGraph(cluster, cluster_bonds, particle_attributes)
                
                # Calculate complexity
                complexity = self.calculate_complexity(graph, particle_attributes)
                
                # Add to catalog
                is_novel, substance_id = self.catalog.add_substance(graph, self.current_time)
                
                # Update novelty tracker
                self.novelty_tracker.record_discovery(is_novel, complexity)
    
    def calculate_complexity(self, graph: MolecularGraph, 
                           particle_attributes: Dict[int, np.ndarray]) -> float:
        """Calculate complexity of a molecular structure"""
        from .metrics import ComplexityAnalyzer
        
        return ComplexityAnalyzer.calculate_emergent_complexity(graph, particle_attributes)
    
    def update_metrics(self):
        """Update all metrics"""
        # Update particle metrics
        self.metrics.update_particle_metrics(
            self.particles.active,
            self.particles.attributes,
            self.particles.energy,
            self.particles.particle_count[None]
        )
        
        # Update bond metrics
        self.metrics.update_bond_metrics(
            self.binding.bond_matrix,
            self.particles.particle_count[None]
        )
        
        # Update cluster metrics
        self.metrics.update_cluster_metrics(
            self.binding.cluster_sizes,
            self.particles.particle_count[None]
        )
        
        # Update energy field metrics
        energy_field = self.energy_manager.energy_system.energy_field
        self.metrics.update_energy_field_metrics(
            energy_field,
            self.config.grid_width,
            self.config.grid_height
        )
        
        # Record metrics
        additional_metrics = {
            'novelty_rate': self.novelty_tracker.get_novelty_rate(),
            'discovery_rate': self.novelty_tracker.get_discovery_rate(),
            'total_substances': len(self.catalog.substances),
            'health_score': self.aggregator.get_health_score()
        }
        
        self.metrics.record_metrics(additional_metrics)
        self.aggregator.update_aggregated_stats()
    
    def start(self):
        """Start simulation"""
        self.is_running = True
        self.is_paused = False
    
    def pause(self):
        """Pause simulation"""
        self.is_paused = True
    
    def resume(self):
        """Resume simulation"""
        self.is_paused = False
    
    def stop(self):
        """Stop simulation"""
        self.is_running = False
        self.is_paused = False
    
    def reset(self):
        """Reset simulation to initial state"""
        self.stop()
        self.current_time = 0.0
        self.step_count = 0
        
        # Reset all components
        self.grid.reset()
        self.particles.reset()
        self.potentials.reset()
        self.binding.reset()
        self.graphs.reset()
        self.catalog.clear()
        self.metrics.reset()
        self.novelty_tracker = NoveltyTracker()
        self.energy_manager = EnergyManager(self.config)
        self.rng.reset(self.config.seed)
        
        # Reinitialize
        self.initialize_simulation()
    
    def get_simulation_state(self) -> Dict:
        """Get current simulation state"""
        return {
            'current_time': self.current_time,
            'step_count': self.step_count,
            'is_running': self.is_running,
            'is_paused': self.is_paused,
            'mode': self.config.mode,
            'particle_count': self.particles.particle_count[None],
            'bond_count': self.binding.get_bond_stats()['num_bonds'],
            'cluster_count': self.binding.get_cluster_stats()['num_clusters'],
            'novelty_rate': self.novelty_tracker.get_novelty_rate(),
            'total_substances': len(self.catalog.substances),
            'health_score': self.aggregator.get_health_score()
        }
    
    def get_visualization_data(self) -> Dict:
        """Get data for visualization"""
        # Get particle data
        positions, attributes, active_mask = self.particles.get_active_particles()
        
        # Get energy field
        energy_field = self.energy_manager.energy_system.energy_field.to_numpy()
        
        # Get bonds
        bonds = self.binding.get_bonds()
        
        # Get clusters
        clusters = self.binding.get_clusters()
        
        return {
            'particles': {
                'positions': positions.tolist(),
                'attributes': attributes.tolist(),
                'active_mask': active_mask.tolist()
            },
            'energy_field': energy_field.tolist(),
            'bonds': bonds,
            'clusters': clusters,
            'metrics': self.aggregator.get_aggregated_stats()
        }
    
    def get_novel_substances(self, count: int = 10) -> List[Dict]:
        """Get recent novel substances"""
        recent_substances = self.catalog.get_recent_substances(count)
        
        return [substance.to_dict() for substance in recent_substances]
    
    def save_snapshot(self, filename: str):
        """Save simulation snapshot"""
        snapshot = {
            'config': self.config.dict(),
            'current_time': self.current_time,
            'step_count': self.step_count,
            'particles': self.particles.get_active_particles(),
            'bonds': self.binding.get_bonds(),
            'clusters': self.binding.get_clusters(),
            'energy_field': self.energy_manager.energy_system.energy_field.to_numpy().tolist(),
            'catalog': {
                'substances': {k: v.to_dict() for k, v in self.catalog.substances.items()},
                'statistics': self.catalog.get_catalog_stats()
            },
            'metrics': self.aggregator.get_aggregated_stats()
        }
        
        import json
        with open(filename, 'w') as f:
            json.dump(snapshot, f, indent=2)
    
    def load_snapshot(self, filename: str):
        """Load simulation snapshot"""
        import json
        with open(filename, 'r') as f:
            snapshot = json.load(f)
        
        # Load configuration
        self.config = SimulationConfig(**snapshot['config'])
        
        # Load simulation state
        self.current_time = snapshot['current_time']
        self.step_count = snapshot['step_count']
        
        # Load particles
        self.particles.reset()
        # ... (implement particle loading)
        
        # Load other components
        # ... (implement loading for other components)
        
        # Update metrics
        self.update_metrics()
