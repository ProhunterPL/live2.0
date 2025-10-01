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
from .diagnostics import DiagnosticsLogger
from ..io.snapshot import SnapshotSerializer

class SimulationStepper:
    """Main simulation coordinator"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.current_time = 0.0
        self.step_count = 0
        self._current_dt = float(config.dt)
        
        # Pre-allocate energy amounts field to avoid recreation in each step
        self.energy_amounts = ti.field(dtype=ti.f32, shape=(config.max_particles,))
        
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
        
        # Energy conservation monitoring
        self.initial_energy = 0.0
        self.energy_history = []
        self.max_energy_history = 1000  # Keep last 1000 energy values
        self.energy_conservation_threshold = 0.05  # 5% energy drift threshold
        
        # Metrics aggregator
        self.aggregator = MetricsAggregator()
        self.aggregator.set_components(self.metrics, self.novelty_tracker, self.catalog)
        
        # Diagnostics logger
        diagnostics_enabled = getattr(config, 'enable_diagnostics', True)
        diagnostics_dir = getattr(config, 'diagnostics_dir', 'diagnostics')
        self.diagnostics = DiagnosticsLogger(
            output_dir=diagnostics_dir,
            enabled=diagnostics_enabled
        )
        
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
        
        # Initialize energy conservation monitoring
        self.initial_energy = self._get_total_energy()
        self.energy_history = [self.initial_energy]
    
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
        """Perform one simulation step with adaptive timestep control"""
        if dt is None:
            dt = self._compute_adaptive_timestep()
        
        # Store previous state for error estimation
        prev_energy = self._get_total_energy()
        
        # Perform step with current dt
        self._perform_step(dt)
        
        # Estimate error and adjust timestep if needed
        current_energy = self._get_total_energy()
        energy_error = abs(current_energy - prev_energy) / max(prev_energy, 1e-6)
        
        # Adaptive timestep control based on energy conservation
        if energy_error > 0.01:  # 1% energy error threshold
            self._current_dt *= 0.8  # Reduce timestep
        elif energy_error < 0.001:  # 0.1% energy error threshold
            self._current_dt *= 1.1  # Increase timestep
        
        # Clamp timestep to reasonable bounds
        base_dt = float(self.config.dt)
        self._current_dt = max(base_dt * 0.1, min(base_dt * 2.0, self._current_dt))
    
    def _compute_adaptive_timestep(self):
        """Compute adaptive timestep based on multiple factors"""
        base_dt = float(self.config.dt)
        
        # Factor 1: Energy field stability
        recent_mean = self.aggregator.aggregated_stats.get('energy_field_mean', None)
        if recent_mean is not None:
            target = 0.1
            energy_factor = max(0.25, min(1.0, target / max(recent_mean, 1e-6)))
        else:
            energy_factor = 1.0
        
        # Factor 2: Particle velocity (CFL condition)
        max_velocity = self._get_max_particle_velocity()
        if max_velocity > 0:
            cfl_dt = 0.1 / max_velocity  # CFL number = 0.1
            cfl_factor = min(1.0, cfl_dt / base_dt)
        else:
            cfl_factor = 1.0
        
        # Factor 3: Force magnitude
        max_force = self._get_max_force_magnitude()
        if max_force > 0:
            force_dt = 0.01 / max_force  # Force-based timestep limit
            force_factor = min(1.0, force_dt / base_dt)
        else:
            force_factor = 1.0
        
        # Combine factors (use most restrictive)
        combined_factor = min(energy_factor, cfl_factor, force_factor)
        self._current_dt = base_dt * combined_factor
        
        return self._current_dt
    
    def _perform_step(self, dt: float):
        """Perform the actual simulation step"""
        
        if not self.is_running or self.is_paused:
            return
        
        # Log only every 100 steps to reduce overhead
        if self.step_count % 100 == 0:
            print(f"STEP {self.step_count + 1}: Starting - sim_time={self.current_time:.6f}, dt={dt:.6f}")
        
        try:
            # Update energy system
            self.energy_manager.update(dt)

            if self.config.mode == "preset_prebiotic":
                # Synchronize preset energy field with main energy system (read-only copy)
                if self.preset_simulator is not None:
                    energy_np = self.energy_manager.energy_system.energy_field.to_numpy()
                    self.preset_simulator.energy_field.from_numpy(energy_np)
                    # Step preset chemistry
                    self.preset_simulator.step(dt)
                # Metrics for preset mode
                self.update_metrics()
            else:
                # Open chemistry branch
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
                
                # Apply bond spring forces (if enabled)
                enable_spring_forces = getattr(self.open_chemistry_config, 'enable_spring_forces', True)
                if enable_spring_forces:
                    self.binding.apply_bond_forces(
                        self.particles.positions,
                        self.particles.velocities,
                        self.potentials.forces,  # Add to existing forces
                        self.particles.active,
                        self.particles.particle_count[None]
                    )
                
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
                
                # Update energy field (legacy decay path)
                self.grid.decay_energy_field(self.config.energy_decay)
                
                # Add energy from sources to particles (now optimized)
                self.add_energy_to_particles()
                
                # Re-enable operations one by one for performance testing
                
                # Update metrics (optimized)
                self.update_metrics()
                
                # Log diagnostics (periodically to reduce overhead)
                diag_freq = getattr(self.config, 'diagnostics_frequency', 10)
                if self.step_count % diag_freq == 0:
                    self._log_diagnostics()
                
                # Skip expensive operations for now
                # self.apply_mutations(dt)
                # self.update_graph_representation()  
                # self.detect_novel_substances()
            
        except Exception as e:
            print(f"STEP {self.step_count + 1}: ERROR - {e}")
            import traceback
            traceback.print_exc()
            raise  # Re-raise to see the full error
        
        # Update time and step count
        self.current_time += dt
        self.step_count += 1
        
        # Update energy conservation monitoring
        self._update_energy_conservation()
        
        # Log completion only every 100 steps
        if self.step_count % 100 == 0:
            energy_drift = self._get_energy_drift()
            print(f"STEP {self.step_count}: COMPLETED - sim_time={self.current_time:.6f}, step_count={self.step_count}, energy_drift={energy_drift:.4f}%")
    
    def _get_total_energy(self):
        """Calculate total energy of the system"""
        total_energy = 0.0
        
        # Energy field energy
        energy_field = self.energy_manager.energy_system.energy_field.to_numpy()
        total_energy += float(energy_field.sum())
        
        # Particle kinetic energy
        if self.config.mode == "open_chemistry":
            velocities = self.particles.velocities.to_numpy()
            attributes = self.particles.attributes.to_numpy()
            active = self.particles.active.to_numpy()
            
            for i in range(len(velocities)):
                if active[i] == 1:
                    mass = attributes[i][0]
                    vx, vy = velocities[i]
                    kinetic = 0.5 * mass * (vx * vx + vy * vy)
                    total_energy += kinetic
        
        return total_energy
    
    def _get_max_particle_velocity(self):
        """Get maximum particle velocity for CFL condition"""
        if self.config.mode != "open_chemistry":
            return 0.0
        
        velocities = self.particles.velocities.to_numpy()
        active = self.particles.active.to_numpy()
        
        max_vel = 0.0
        for i in range(len(velocities)):
            if active[i] == 1:
                vx, vy = velocities[i]
                vel_mag = (vx * vx + vy * vy) ** 0.5
                max_vel = max(max_vel, vel_mag)
        
        return max_vel
    
    def _get_max_force_magnitude(self):
        """Get maximum force magnitude for timestep control"""
        if self.config.mode != "open_chemistry":
            return 0.0
        
        forces = self.potentials.forces.to_numpy()
        active = self.particles.active.to_numpy()
        
        max_force = 0.0
        for i in range(len(forces)):
            if active[i] == 1:
                fx, fy = forces[i]
                force_mag = (fx * fx + fy * fy) ** 0.5
                max_force = max(max_force, force_mag)
        
        return max_force
    
    def _update_energy_conservation(self):
        """Update energy conservation history"""
        current_energy = self._get_total_energy()
        
        # Initialize initial energy on first call
        if self.initial_energy == 0.0 and current_energy > 0.0:
            self.initial_energy = current_energy
        
        # Add to history
        self.energy_history.append(current_energy)
        
        # Limit history size
        if len(self.energy_history) > self.max_energy_history:
            self.energy_history.pop(0)
        
        # Check for energy drift violation
        drift = self._get_energy_drift()
        if drift > self.energy_conservation_threshold * 100.0:  # Convert threshold to percentage
            # Log warning only every 1000 steps to avoid spam
            if self.step_count % 1000 == 0:
                print(f"WARNING: Energy drift ({drift:.2f}%) exceeds threshold ({self.energy_conservation_threshold*100.0:.2f}%)")
    
    def _get_energy_drift(self):
        """Calculate energy drift percentage from initial energy"""
        if self.initial_energy == 0.0:
            return 0.0
        
        current_energy = self._get_total_energy()
        drift = abs(current_energy - self.initial_energy) / self.initial_energy * 100.0
        return drift
    
    def add_energy_to_particles(self):
        """Add energy from field to particles"""
        # Use pre-allocated field instead of creating new one
        
        # Clear energy amounts
        for i in range(self.config.max_particles):
            self.energy_amounts[i] = 0.0
        
        # Add energy from field
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                pos = self.particles.positions[i]
                energy = self.energy_manager.get_energy_at_position(pos.to_numpy())
                self.energy_amounts[i] = energy
        
        # Apply energy to particles
        self.particles.add_energy(self.energy_amounts)
    
    def apply_mutations(self, dt: float):
        """Apply mutations to particles in high energy regions"""
        mutation_rate = self.open_chemistry_config.mutation_rate
        mutation_strength = self.open_chemistry_config.mutation_strength
        
        # Throttle mutation application to avoid heavy work every step
        mutation_interval = getattr(self.open_chemistry_config, "mutation_interval", 10)
        if mutation_interval is None or mutation_interval < 1:
            mutation_interval = 10
        if (self.step_count % mutation_interval) != 0:
            return
        
        # Find high energy regions (may be many)
        high_energy_regions = self.energy_manager.get_high_energy_regions()
        if not high_energy_regions:
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
        
        # Log mutations only every 50 steps
        if self.step_count % 50 == 0:
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
    
    def _log_diagnostics(self):
        """Log diagnostics data for current step"""
        if not self.diagnostics.enabled:
            return
        
        try:
            # Prepare simulation data for diagnostics
            simulation_data = {}
            
            # Get particle positions
            particle_count = self.particles.particle_count[None]
            if particle_count > 0:
                positions = self.particles.positions.to_numpy()[:particle_count]
                simulation_data['positions'] = positions
                
                # Get particle attributes
                attributes = self.particles.attributes.to_numpy()[:particle_count]
                simulation_data['attributes'] = attributes
                
                # Get particle energies
                particle_energies = self.particles.energy.to_numpy()[:particle_count]
                simulation_data['particle_energies'] = particle_energies
            else:
                simulation_data['positions'] = np.array([])
                simulation_data['attributes'] = np.array([])
                simulation_data['particle_energies'] = np.array([])
            
            # Get bonds
            bonds = []
            bond_matrix = self.binding.bond_matrix.to_numpy()
            for i in range(particle_count):
                for j in range(i + 1, particle_count):
                    strength = bond_matrix[i, j]
                    if strength > 0:
                        bonds.append((i, j, strength))
            simulation_data['bonds'] = bonds
            
            # Get clusters
            clusters = {}
            cluster_stats = self.binding.get_cluster_stats()
            if 'clusters' in cluster_stats:
                clusters = cluster_stats['clusters']
            simulation_data['clusters'] = clusters
            
            # Log to diagnostics
            self.diagnostics.log_step(
                step=self.step_count,
                sim_time=self.current_time,
                simulation_data=simulation_data
            )
            
        except Exception as e:
            import logging
            logger = logging.getLogger(__name__)
            logger.warning(f"Failed to log diagnostics: {e}")
    
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
        
        # Close diagnostics logger
        if hasattr(self, 'diagnostics'):
            self.diagnostics.close()
    
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
            'health_score': self.aggregator.get_health_score(),
            'energy_drift': self._get_energy_drift(),
            'total_energy': self._get_total_energy(),
            'adaptive_dt': self._current_dt
        }
    
    def get_visualization_data(self) -> Dict:
        """Get data for visualization"""
        data = {
            'metrics': self.aggregator.get_aggregated_stats()
        }
        
        # Energy field (always present)
        energy_field = self.energy_manager.energy_system.energy_field.to_numpy()
        data['energy_field'] = energy_field.tolist()
        
        import logging
        logger = logging.getLogger(__name__)
        logger.info(f"DEBUG: get_visualization_data - mode: {self.config.mode}")
        logger.info(f"DEBUG: energy_field shape: {energy_field.shape}, max: {energy_field.max()}, min: {energy_field.min()}")
        
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
                        logger.info(f"DEBUG: preset mode - concentration field shape: {conc[first_key].shape}")
                except Exception as e:
                    logger.info(f"DEBUG: preset mode error: {e}")
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
            
            logger.info(f"DEBUG: open chemistry mode - particles: {len(positions)}, bonds: {len(bonds)}")
            if len(positions) > 0:
                logger.info(f"DEBUG: sample particle position: {positions[0]}")
        
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
