"""
Main simulation stepper for Live 2.0
Coordinates all simulation components and manages time evolution
"""

import taichi as ti
import numpy as np
import time
from typing import Dict, List, Optional, Any
from ..config import SimulationConfig, PresetPrebioticConfig, OpenChemistryConfig

from .grid import Grid
from .particles import ParticleSystem
from .potentials import PotentialSystem
from .binding import BindingSystem
from .graphs import GraphProcessor, MolecularGraph
from .catalog import SubstanceCatalog
from .metrics import MetricsCollector, NoveltyTracker, MetricsAggregator
from .energy import EnergyManager
from .rng import RNG
from .fields_ca import PresetPrebioticSimulator
from ..io.snapshot import SnapshotSerializer

class SimulationStepper:
    """Main simulation coordinator"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.current_time = 0.0
        self.step_count = 0
        self._current_dt = float(config.dt)
        
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
        self.preset_simulator = None
        
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
        # Initialize the preset fields/chemistry simulator
        self.preset_simulator = PresetPrebioticSimulator(
            self.preset_config,
            self.config.grid_width,
            self.config.grid_height,
        )
        
        # Add energy sources for preset mode as well
        for _ in range(self.open_chemistry_config.energy_sources):
            px, py = self.rng.py_next_vector2(0, float(self.config.grid_width))
            pos = ti.Vector([px, py])
            intensity = self.open_chemistry_config.energy_intensity
            radius = 10.0
            self.energy_manager.add_energy_source(
                (float(pos[0]), float(pos[1])), intensity, radius, duration=100.0
            )
    
    def initialize_open_chemistry_mode(self):
        """Initialize open chemistry mode"""
        # Add initial particles with random properties
        num_initial_particles = min(100, self.config.max_particles)
        
        for i in range(num_initial_particles):
            # Random position
            px, py = self.rng.py_next_vector2(0, float(self.config.grid_width))
            pos = ti.Vector([px, py])
            
            # Random velocity
            vx, vy = self.rng.py_next_vector2(-0.5, 0.5)
            vel = ti.Vector([vx, vy])
            
            # Random attributes
            mass = float(self.rng.py_next_range(0.5, 2.0))
            charge_x = float(self.rng.py_next_range(-1.0, 1.0))
            charge_y = float(self.rng.py_next_range(-1.0, 1.0))
            charge_z = float(self.rng.py_next_range(-1.0, 1.0))
            attributes = ti.Vector([mass, charge_x, charge_y, charge_z])
            
            # Random type
            type_id = self.particles.register_particle_type(
                name=f"particle_{i}",
                mass=mass,
                charge=(charge_x, charge_y, charge_z)
            )
            
            # Add particle
            self.particles.add_particle_py(pos, vel, attributes, type_id, 2, 1.0)
        
        # Add energy sources
        for _ in range(self.open_chemistry_config.energy_sources):
            px, py = self.rng.py_next_vector2(0, float(self.config.grid_width))
            pos = ti.Vector([px, py])
            intensity = self.open_chemistry_config.energy_intensity
            radius = 15.0
            
            self.energy_manager.add_energy_source(
                (float(pos[0]), float(pos[1])), intensity, radius, duration=200.0
            )
    
    def step(self, dt: float = None):
        """Perform one simulation step"""
        if dt is None:
            # Adaptive timestep based on recent energy mean (simple controller)
            recent_mean = self.aggregator.aggregated_stats.get('energy_field_mean', None)
            base_dt = float(self.config.dt)
            if recent_mean is None:
                self._current_dt = base_dt
            else:
                # Scale dt inversely with energy; clamp to [base_dt/4, base_dt]
                # target_mean ~ 0.1
                target = 0.1
                factor = max(0.25, min(1.0, target / max(recent_mean, 1e-6)))
                self._current_dt = base_dt * factor
            dt = self._current_dt
        
        if not self.is_running or self.is_paused:
            return
        
        print(f"STEP {self.step_count + 1}: Starting - sim_time={self.current_time:.6f}, dt={dt:.6f}")
        
        try:
            # Update energy system
            print(f"STEP {self.step_count + 1}: Updating energy system")
            self.energy_manager.update(dt)

            if self.config.mode == "preset_prebiotic":
                # Synchronize preset energy field with main energy system (read-only copy)
                if self.preset_simulator is not None:
                    energy_np = self.energy_manager.energy_system.energy_field.to_numpy()
                    self.preset_simulator.energy_field.from_numpy(energy_np)
                    # Step preset chemistry
                    print(f"STEP {self.step_count + 1}: Stepping preset prebiotic simulator")
                    self.preset_simulator.step(dt)
                # Metrics for preset mode
                print(f"STEP {self.step_count + 1}: Updating metrics (preset)")
                self.update_metrics()
            else:
                # Open chemistry branch
                # Update particle positions
                print(f"STEP {self.step_count + 1}: Updating particle positions")
                self.particles.update_positions(dt)
                
                # Update spatial hash
                print(f"STEP {self.step_count + 1}: Updating spatial hash")
                self.grid.update_spatial_hash()
                
                # Compute forces
                print(f"STEP {self.step_count + 1}: Computing forces")
                self.potentials.compute_forces(
                    self.particles.positions,
                    self.particles.attributes,
                    self.particles.active,
                    self.particles.particle_count[None]
                )
                
                # Apply forces
                print(f"STEP {self.step_count + 1}: Applying forces")
                self.particles.apply_forces(self.potentials.forces, dt)
                
                # Update binding system
                print(f"STEP {self.step_count + 1}: Updating binding system")
                self.binding.update_bonds(
                    self.particles.positions,
                    self.particles.attributes,
                    self.particles.active,
                    self.particles.particle_count[None],
                    dt
                )
                
                # Update clusters
                print(f"STEP {self.step_count + 1}: Updating clusters")
                self.binding.update_clusters(
                    self.particles.active,
                    self.particles.particle_count[None]
                )
                
                # Apply periodic boundary conditions
                print(f"STEP {self.step_count + 1}: Applying periodic boundary conditions")
                self.grid.apply_periodic_boundary()
                
                # Update energy field (legacy decay path)
                print(f"STEP {self.step_count + 1}: Updating energy field (grid decay)")
                self.grid.decay_energy_field(self.config.energy_decay)
                
                # Add energy from sources to particles
                print(f"STEP {self.step_count + 1}: Adding energy to particles")
                self.add_energy_to_particles()
                
                # Apply mutations in high energy regions
                print(f"STEP {self.step_count + 1}: Applying mutations")
                self.apply_mutations(dt)
                
                # Update graph representation
                print(f"STEP {self.step_count + 1}: Updating graph representation")
                self.update_graph_representation()
                
                # Detect novel substances
                print(f"STEP {self.step_count + 1}: Detecting novel substances")
                self.detect_novel_substances()
                
                # Update metrics
                print(f"STEP {self.step_count + 1}: Updating metrics")
                self.update_metrics()
            
        except Exception as e:
            print(f"STEP {self.step_count + 1}: ERROR - {e}")
            import traceback
            traceback.print_exc()
            raise  # Re-raise to see the full error
        
        # Update time and step count
        self.current_time += dt
        self.step_count += 1
        print(f"STEP {self.step_count}: COMPLETED - sim_time={self.current_time:.6f}, step_count={self.step_count}")
    
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
        
        # Throttle mutation application to avoid heavy work every step
        mutation_interval = getattr(self.open_chemistry_config, "mutation_interval", 5)
        if mutation_interval is None or mutation_interval < 1:
            mutation_interval = 5
        if (self.step_count % mutation_interval) != 0:
            print(f"STEP {self.step_count + 1}: Mutations skipped (interval={mutation_interval})")
            return
        
        # Find high energy regions (may be many)
        high_energy_regions = self.energy_manager.get_high_energy_regions()
        if not high_energy_regions:
            print(f"STEP {self.step_count + 1}: No high energy regions for mutations")
            return
        
        # Cap the number of regions processed per step
        max_regions_per_step = 64
        selected_regions = high_energy_regions[:max_regions_per_step]
        
        # Limit neighbors considered per region and total mutations per step
        max_neighbors_per_region = 64
        max_mutations_per_step = 200
        neighbors = ti.field(dtype=ti.i32, shape=(max_neighbors_per_region,))
        total_applied = 0
        
        for region_pos in selected_regions:
            if total_applied >= max_mutations_per_step:
                break
            region_ti = ti.Vector(region_pos)
            count = int(self.grid.get_neighbors(region_ti, 5.0, neighbors))
            count = max(0, min(count, max_neighbors_per_region))
            if count == 0:
                continue
            # Pull neighbors to numpy for safe Python-side iteration
            neighbor_indices = neighbors.to_numpy()[:count]
            
            for particle_idx in neighbor_indices:
                if total_applied >= max_mutations_per_step:
                    break
                if self.particles.active[particle_idx] == 1:
                    if self.rng.py_next() < mutation_rate * dt:
                        self.particles.mutate_particle(
                            int(particle_idx), float(mutation_strength), float(self.current_time), self.rng
                        )
                        total_applied += 1
        
        print(
            f"STEP {self.step_count + 1}: Mutations applied={total_applied}, "
            f"regions_considered={len(selected_regions)}/{len(high_energy_regions)}, "
            f"interval={mutation_interval}"
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
        
        # Convert dict to list for ComplexityAnalyzer
        attributes_list = list(particle_attributes.values())
        return ComplexityAnalyzer.calculate_emergent_complexity(graph, attributes_list)
    
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
        data = {
            'metrics': self.aggregator.get_aggregated_stats()
        }
        
        # Energy field (always present)
        energy_field = self.energy_manager.energy_system.energy_field.to_numpy()
        data['energy_field'] = energy_field.tolist()
        
        if self.config.mode == "preset_prebiotic":
            # Provide concentration fields for preset mode
            if self.preset_simulator is not None:
                # Extract only one channel (selected species index 0 by default) to reduce bandwidth
                try:
                    conc = self.preset_simulator.concentration_fields.get_all_concentrations()
                    # pick first available species
                    first_key = next(iter(conc.keys())) if conc else None
                    if first_key is not None:
                        data['concentration_view'] = conc[first_key].tolist()
                        data['concentration_species'] = list(conc.keys())
                        # Back-compat: also provide single-key concentrations for older clients
                        data['concentrations'] = {first_key: conc[first_key].tolist()}
                except Exception:
                    pass
        else:
            # Provide particle/bond/cluster data for open chemistry
            positions, velocities, attributes, active_mask = self.particles.get_active_particles()
            bonds = self.binding.get_bonds()
            clusters = self.binding.get_clusters()
            data['particles'] = {
                'positions': positions.tolist(),
                'attributes': attributes.tolist(),
                'active_mask': active_mask.tolist()
            }
            data['bonds'] = bonds
            data['clusters'] = clusters
        
        return data
    
    def get_novel_substances(self, count: int = 10) -> List[Dict]:
        """Get recent novel substances"""
        recent_substances = self.catalog.get_recent_substances(count)
        
        return [substance.to_dict() for substance in recent_substances]
    
    def save_snapshot(self, filename: str):
        """Save simulation snapshot using serializer"""
        snapshot = SnapshotSerializer.serialize_simulation(self)
        import json
        with open(filename, 'w') as f:
            json.dump(snapshot, f, indent=2)
    
    def load_snapshot(self, filename: str):
        """Load simulation snapshot using serializer"""
        import json
        with open(filename, 'r') as f:
            snapshot = json.load(f)
        if not SnapshotSerializer.validate_snapshot(snapshot):
            raise ValueError("Invalid snapshot data")
        SnapshotSerializer.deserialize_simulation(snapshot, self)
