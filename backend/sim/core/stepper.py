"""
Main simulation stepper for Live 2.0
Coordinates all simulation components and manages time evolution
"""

import taichi as ti
import numpy as np
import time
import logging
import sys
from typing import Dict, List, Optional, Any
from .memory_manager import register_taichi_field, get_memory_stats, optimize_memory, should_cleanup_memory

# Set stepper logger to INFO level
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
from ..config import SimulationConfig, PresetPrebioticConfig, OpenChemistryConfig

logger = logging.getLogger(__name__)

from .grid import Grid
from .particles import ParticleSystem
from .potentials import PotentialSystem
from .binding import BindingSystem
from .graphs import GraphProcessor, MolecularGraph
from .catalog import SubstanceCatalog
from .metrics import MetricsCollector, NoveltyTracker, MetricsAggregator, PerformanceMonitor
from .energy import EnergyManager
from .rng import RNG
from .fields_ca import PresetPrebioticSimulator
from .diagnostics import DiagnosticsLogger
from ..io.snapshot import SnapshotManager, SnapshotSerializer
from .integrators import SymplecticIntegrators, AdaptiveIntegrator
from .thermodynamics import ThermodynamicValidator

@ti.data_oriented
class SimulationStepper:
    """Main simulation coordinator"""
    
    def __init__(self, config: SimulationConfig):
        # print(f"DEBUG: SimulationStepper.__init__ called with mode={config.mode}")
        self.config = config
        self.current_time = 0.0
        self.step_count = 0
        self._current_dt = float(config.dt)
        
        # Pre-allocate energy amounts field to avoid recreation in each step
        self.energy_amounts = ti.field(dtype=ti.f32, shape=(config.max_particles,))
        self.debug_energy_total = ti.field(dtype=ti.f32, shape=())
        self.debug_particles_with_energy = ti.field(dtype=ti.i32, shape=())
        
        # Initialize components
        self.grid = Grid(config)
        self.particles = ParticleSystem(config)
        self.potentials = PotentialSystem(config)
        self.binding = BindingSystem(config)
        self.graphs = GraphProcessor(config.max_particles)
        self.catalog = SubstanceCatalog()
        self.metrics = MetricsCollector(config.max_particles)
        self.novelty_tracker = NoveltyTracker(window_size=config.novelty_window)
        self.performance_monitor = PerformanceMonitor()
        self.energy_manager = EnergyManager(config)
        self.rng = RNG(config.seed)
        
        # Energy conservation monitoring
        self.initial_energy = 0.0
        from collections import deque
        self.energy_history = deque(maxlen=1000)  # MEMORY FIX: Use deque instead of list
        self.max_energy_history = 1000  # Keep last 1000 energy values
        self.energy_conservation_threshold = 0.05  # 5% energy drift threshold
        
        # Metrics aggregator
        self.aggregator = MetricsAggregator()
        self.aggregator.set_components(self.metrics, self.novelty_tracker, self.catalog)
        
        # Energy field cache for visualization (reduces GPU->CPU transfers)
        self._energy_field_cache = None
        self._energy_field_cache_step = -999
        self._energy_field_cache_interval = 10
        
        # Initialize particle cache to avoid empty visualization
        self._cached_particles_data = None
        self._cached_bonds_clusters_data = None  # Update every 10 steps
        
        # Diagnostics logger
        diagnostics_enabled = getattr(config, 'enable_diagnostics', True)
        diagnostics_dir = getattr(config, 'diagnostics_dir', 'diagnostics')
        self.diagnostics = DiagnosticsLogger(
            output_dir=diagnostics_dir,
            enabled=diagnostics_enabled
        )
        
        # Snapshot manager with image generation - use absolute path
        import os
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        snapshots_path = os.path.join(project_root, "snapshots")
        self.snapshot_manager = SnapshotManager(snapshot_dir=snapshots_path)
        
        # Thermodynamic validator
        self.enable_validation = getattr(config, 'enable_thermodynamic_validation', False)
        self.validator = ThermodynamicValidator(config) if self.enable_validation else None
        self.validation_interval = getattr(config, 'validate_every_n_steps', 10000)
        self.validation_log = []
        
        # MEMORY MONITORING: Register Taichi fields for tracking
        register_taichi_field('particles_positions', self.particles.positions)
        register_taichi_field('particles_velocities', self.particles.velocities)
        register_taichi_field('particles_attributes', self.particles.attributes)
        register_taichi_field('particles_active', self.particles.active)
        register_taichi_field('energy_field', self.energy_manager.energy_system.energy_field)
        
        # Symplectic integrators for better energy conservation
        self.integrator = AdaptiveIntegrator(config.max_particles, method="verlet")
        self.use_symplectic = getattr(config, 'use_symplectic_integrator', True)
        
        # Memory management
        self.last_cleanup_time = time.time()
        self.cleanup_interval = 60  # Clean up every 1 minute (reduced from 5 minutes)
        
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
        self.metrics.reset()
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
        # Add initial particles with random properties - INCREASED for larger clusters
        num_initial_particles = min(1000, self.config.max_particles)  # INCREASED from 500 to 1000
        
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
        
        # Debug print removed for performance
        
        # Initialize visualization cache with initial particles
        positions, velocities, attributes, active_mask, energies = self.particles.get_active_particles()
        self._cached_particles_data = {
            'positions': positions,
            'attributes': attributes,
            'active_mask': active_mask,
            'energies': energies
        }
        
        # Initialize bonds/clusters cache
        bonds = self.binding.get_bonds()
        clusters = self.binding.get_clusters()
        self._cached_bonds_clusters_data = {
            'bonds': bonds,
            'clusters': clusters
        }
        
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
        # MEMORY MONITORING: Check memory usage every 1000 steps (reduced frequency)
        memory_check_interval = getattr(self.config, 'memory_check_interval', 1000)
        if memory_check_interval > 0 and self.step_count % memory_check_interval == 0:
            memory_stats = get_memory_stats()
            if memory_stats['current_memory_mb'] > 5000:  # Warn only if over 5GB
                logger.warning(f"High memory usage at step {self.step_count}: {memory_stats['current_memory_mb']:.1f} MB")
            
            # Perform memory cleanup only if really needed
            if should_cleanup_memory():
                optimize_memory()
        
        # Garbage collection every 5000 steps (less frequent)
        gc_interval = getattr(self.config, 'gc_interval', 5000)
        if gc_interval > 0 and self.step_count % gc_interval == 0:
            import gc
            gc.collect()
        
        # print("STEP FUNCTION CALLED")  # Disabled to prevent infinite loop
        # Debug prints disabled to prevent infinite loop
        # print(f"STEP ENTRY: step_count={self.step_count}, is_running={self.is_running}, is_paused={self.is_paused}")
        # sys.stdout.flush()
        # print(f"DEBUG: dt param = {dt}")
        # sys.stdout.flush()

        # Debug prints disabled to prevent infinite loop
        # if self.step_count < 5:
        #     print(f"STEP {self.step_count + 1}: step() called with dt={dt}")
        #     sys.stdout.flush()

        if dt is None:
            # print(f"DEBUG: Computing adaptive timestep")
            # sys.stdout.flush()
            try:
                dt = self._compute_adaptive_timestep()
                # print(f"DEBUG: Computed dt = {dt}")
                # sys.stdout.flush()
            except Exception as e:
                logger.error(f"Adaptive timestep computation failed: {e}")
                sys.stdout.flush()
                return

        # if self.step_count < 5:
        #     print(f"STEP {self.step_count + 1}: computed dt={dt:.6f}")
        #     sys.stdout.flush()

        # Store previous state for thermodynamic validation
        state_before = self._create_state_snapshot()
        energy_injected = 0.0
        energy_dissipated = 0.0

        # Store previous state for error estimation
        # print(f"DEBUG: About to get prev_energy")
        # sys.stdout.flush()
        try:
            prev_energy = self._get_total_energy()
            # print(f"DEBUG: prev_energy = {prev_energy}")
        except Exception as e:
            logger.error(f"Energy calculation failed: {e}")
            return

        # Perform step with current dt
        # print(f"DEBUG: About to call _perform_step with dt={dt}")
        self._perform_step(dt)
        # print(f"DEBUG: _perform_step returned")

        # print(f"DEBUG: _perform_step completed, about to get current_energy")

        # Estimate error and adjust timestep if needed
        current_energy = self._get_total_energy()
        # print(f"DEBUG: current_energy = {current_energy}, prev_energy = {prev_energy}")
        energy_error = abs(current_energy - prev_energy) / max(prev_energy, 1e-6)
        
        # Adaptive timestep control based on energy conservation
        if energy_error > 0.01:  # 1% energy error threshold
            self._current_dt *= 0.8  # Reduce timestep
        elif energy_error < 0.001:  # 0.1% energy error threshold
            self._current_dt *= 1.1  # Increase timestep
        
        # Clamp timestep to reasonable bounds
        base_dt = float(self.config.dt)
        self._current_dt = max(base_dt * 0.1, min(base_dt * 2.0, self._current_dt))
        
        # Thermodynamic validation (every N steps) - SIMPLIFIED to prevent hangs
        if self.enable_validation and self.validator is not None and self.step_count % self.validation_interval == 0:
            validation_start = time.time()
            
            try:
                # SMART VALIDATION: Different tests at different frequencies (GROMACS/NAMD best practices)
                # Energy+Momentum: every call (~2ms), M-B: every 20k steps, Entropy: every 50k steps
                state_after = self._create_state_snapshot()
                validation_results = self.validator.validate_smart(
                    state_before, state_after, energy_injected, energy_dissipated, 
                    self.step_count
                )
                
                validation_time = time.time() - validation_start
                
                # Log timing only for full validation (M-B or Entropy, >100ms)
                if validation_time > 0.1:
                    logger.info(f"Full thermodynamic validation at step {self.step_count}: {validation_time*1000:.1f}ms")
                
                # Log validation results if successful
                if validation_results and validation_time < 5.0:  # Increased timeout for full validation
                    self.validator.log_validation_results(validation_results)
                
            except Exception as e:
                logger.error(f"Thermodynamic validation failed at step {self.step_count}: {e}")
                # Continue simulation even if validation fails
        
        # Periodic memory cleanup (more aggressive)
        current_time = time.time()
        if current_time - self.last_cleanup_time > self.cleanup_interval:
            # MEMORY FIX: More aggressive cleanup to prevent accumulation
            self.catalog.cleanup_old_data(max_age_hours=0.25)  # Keep only last 15 minutes (was 30)
            self.last_cleanup_time = current_time
            
            # Force garbage collection every cleanup
            import gc
            gc.collect()
            
            # MEMORY FIX: Clear old validation log entries
            if hasattr(self, 'validator') and self.validator is not None:
                if hasattr(self.validator, 'validation_log') and len(self.validator.validation_log) > 50:
                    self.validator.validation_log = self.validator.validation_log[-50:]  # Keep only last 50
    
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
        
        # Factor 2: Particle velocity (CFL condition) - RELAXED
        max_velocity = self._get_max_particle_velocity()
        if max_velocity > 0:
            cfl_dt = 0.5 / max_velocity  # CFL number = 0.5 (increased from 0.1 for larger timesteps)
            cfl_factor = min(1.0, cfl_dt / base_dt)
        else:
            cfl_factor = 1.0
        
        # Factor 3: Force magnitude - RELAXED
        max_force = self._get_max_force_magnitude()
        if max_force > 0:
            force_dt = 0.05 / max_force  # Force-based timestep limit (increased from 0.01)
            force_factor = min(1.0, force_dt / base_dt)
        else:
            force_factor = 1.0
        
        # Combine factors (use most restrictive) - LESS AGGRESSIVE
        # Don't let dt drop below 20% of base_dt
        combined_factor = max(0.2, min(energy_factor, cfl_factor, force_factor))
        self._current_dt = base_dt * combined_factor
        
        return self._current_dt
    
    def _perform_step(self, dt: float):
        """Perform the actual simulation step"""
        
        # Start performance timing
        self.performance_monitor.start_step_timing()

        if not self.is_running or self.is_paused:
            return

        try:
            self.energy_manager.update(dt)
            self._energy_diffuse(dt)
            self._energy_thermostat()
            
            # Add energy pulses periodically
            pulse_every = getattr(self.config, 'pulse_every', 48)
            if pulse_every and self.step_count % pulse_every == 0:
                self._add_energy_pulse()

            if self.config.mode == "preset_prebiotic":
                # MEMORY LEAK FIX: Use reference instead of to_numpy() copy
                if self.preset_simulator is not None:
                    # Use direct field reference instead of numpy conversion
                    self.preset_simulator.energy_field = self.energy_manager.energy_system.energy_field
                    # Step preset chemistry
                    self.preset_simulator.step(dt)
                # Metrics for preset mode
                self.update_metrics()
            else:
                # Open chemistry branch
                if False:  # Temporarily disable symplectic integrator
                    # Use symplectic integrator for better energy conservation
                    particle_count = self.particles.particle_count[None]
                    
                    if self.step_count < 5:
                        # Debug print removed for performance
                        pass
                    # Update spatial hash first for force computation
                    self.grid.update_spatial_hash()
                    
                    if self.step_count < 5:
                        # Debug print removed for performance
                        pass
                    # Compute forces for Verlet first stage
                    self.potentials.compute_forces(
                        self.particles.positions,
                        self.particles.attributes,
                        self.particles.active,
                        particle_count
                    )
                    
                else:
                    # Original Euler method (fallback)
                    self.particles.update_positions(dt)
                    self.grid.update_spatial_hash()
                    self.potentials.compute_forces(
                        self.particles.positions,
                        self.particles.attributes,
                        self.particles.active,
                        self.particles.particle_count[None]
                    )
                    self.particles.apply_forces(self.potentials.forces, dt)
                
                # Add thermal kick based on local energy
                vmax = getattr(self.open_chemistry_config, 'vmax', 8.0)
                self.particles.thermal_kick(vmax, 0.6, self.energy_manager.energy_system.energy_field)
                
                # PERFORMANCE FIX: These O(n²) operations are killing performance at high particle counts
                # DISABLED for overnight test - too expensive
                # if self.step_count % 50 == 0:
                #     # Add clustering assistance - particles move towards high energy regions
                #     self._assist_clustering()
                #     
                #     # Add strong clustering force - particles move towards center
                #     self._force_clustering_to_center()
                
                # DISABLED: This O(n²) operation is too expensive - causes freeze after 1500 steps
                # if self.step_count % 100 == 0:
                #     self._attract_particles_for_bonding()
                
                # Update metrics for open chemistry mode - HEAVILY THROTTLED
                # DISABLED for overnight test - too expensive with to_numpy() operations
                # Only update at key checkpoints
                # metrics_update_interval = getattr(self.config, 'metrics_update_interval', 300)
                # if self.step_count < 3 or (self.step_count % metrics_update_interval) == 0:
                #     logger.info(f"[STEP {self.step_count}] Calling update_metrics (first call)")
                #     try:
                #         self.update_metrics()
                #         logger.info(f"[STEP {self.step_count}] update_metrics completed")
                #     except Exception as e:
                #         logger.error(f"[STEP {self.step_count}] ERROR in update_metrics: {e}")
                #         import traceback
                #         traceback.print_exc()
                
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
                
                # Update binding system (heavily throttled for performance)
                # OPTIMIZATION: Update bonds every 500 steps (was 150) for overnight test
                if self.step_count % 500 == 0:
                    self.binding.update_bonds(
                        self.particles.positions,
                        self.particles.attributes,
                        self.particles.active,
                        self.particles.particle_count[None],
                        dt
                    )
                
                # OPTIMIZATION: Update clusters every 1000 steps (was 300) for overnight test
                if self.step_count % 1000 == 0:
                    self.binding.update_clusters(
                        self.particles.positions,
                        self.particles.active,
                        self.particles.particle_count[None]
                    )
                
                # Apply periodic boundary conditions
                self.grid.apply_periodic_boundary()
                
                # Additional stability: ensure particles stay within grid bounds
                self._stabilize_particle_positions()
                
                # Update energy field (legacy decay path)
                self.grid.decay_energy_field(self.config.energy_decay)
                
                # Add energy from sources to particles (heavily optimized)
                # Add energy from field to particles (always for first few steps, then heavily throttled)
                energy_update_interval = getattr(self.config, 'energy_update_interval', 50)
                if energy_update_interval is None or energy_update_interval < 1:
                    energy_update_interval = 50
                # Always update energy for first 5 steps to avoid getting stuck
                if self.step_count < 5 or (self.step_count % energy_update_interval) == 0:
                    # logger.info(f"DEBUG: Energy update condition met at step {self.step_count} (step < 10: {self.step_count < 10}, step % {energy_update_interval} == 0: {(self.step_count % energy_update_interval) == 0})")
                    if self.step_count < 5:
                        # Debug print removed for performance
                        pass
                    # Controlled energy transfer - move energy from field to particles (no duplication)
                    self._transfer_energy_to_particles()
                
                # Re-enable operations one by one for performance testing
                
                # Update metrics (OPTIMIZED: configurable interval, default 10000)
                # to_numpy() is very expensive, so minimize calls
                metrics_update_interval = getattr(self.config, 'metrics_update_interval', 10000)
                if metrics_update_interval is None or metrics_update_interval < 1:
                    metrics_update_interval = 10000
                # Skip metrics for first 3 steps, then at specified interval
                if self.step_count > 3 and (self.step_count % metrics_update_interval) == 0:
                    try:
                        self.update_metrics()
                    except Exception as e:
                        logger.error(f"ERROR in update_metrics: {e}")
                
                # Log diagnostics (OPTIMIZED: configurable)
                diag_freq = getattr(self.config, 'diagnostics_frequency', 5000)
                if diag_freq > 0 and self.step_count % diag_freq == 0:
                    self._log_diagnostics()
                
                # Mutations (OPTIMIZED: configurable interval)
                mutation_interval = getattr(self.config, 'mutation_interval', 2000)
                if mutation_interval > 0 and mutation_interval < 999999 and self.step_count % mutation_interval == 0:
                    self.apply_mutations(dt)
                
                # Novelty detection (OPTIMIZED: configurable interval)
                novelty_interval = getattr(self.config, 'novelty_check_interval', 500)  # INCREASED from 100 to 500 for better performance
                if novelty_interval > 0 and novelty_interval < 999999 and self.step_count % novelty_interval == 0:
                    self.detect_novel_substances()
            
            
        except Exception as e:
            logger.error(f"[STEP {self.step_count}] Step execution failed: {e}")
            import traceback
            traceback.print_exc()
            raise  # Re-raise to see the full error
        
        # Update time and step count
        self.current_time += dt
        self.step_count += 1
        
        # End performance timing
        step_time = self.performance_monitor.end_step_timing()
        perf_log_interval = getattr(self.config, 'performance_log_interval', 1000)
        if perf_log_interval > 0 and self.step_count % perf_log_interval == 0:
            logger.info(f"Step {self.step_count} completed in {step_time*1000:.1f}ms")

        # Update energy conservation monitoring
        self._update_energy_conservation()
        
        # Log completion only at specified interval
        if perf_log_interval > 0 and self.step_count % perf_log_interval == 0:
            energy_drift = self._get_energy_drift()
            print(f"Step {self.step_count}: sim_time={self.current_time:.2f}s, energy_drift={energy_drift:.4f}%")
    
    def _create_state_snapshot(self):
        """Create a lightweight snapshot of current simulation state for validation"""
        class StateSnapshot:
            def __init__(self, stepper):
                # CRITICAL: Avoid copying large arrays - use references instead
                # This prevents memory spikes that can cause system hangs
                self.positions = stepper.particles.positions
                self.velocities = stepper.particles.velocities
                self.attributes = stepper.particles.attributes
                self.active = stepper.particles.active
                self.energy_field = stepper.energy_manager.energy_system.energy_field
                self.bond_energy = 0.0  # Placeholder - would need to compute from bonds
                
                # Add particle count to avoid expensive counting operations
                # Use Taichi kernel to count particles instead of to_numpy()
                self.particle_count = stepper._count_active_particles()
        
        return StateSnapshot(self)
    
    @ti.kernel
    def _count_active_particles(self) -> int:
        """Count active particles using Taichi kernel (faster than to_numpy())"""
        count = 0
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                count += 1
        return count
    
    @ti.kernel
    def _compute_energy_field_sum(self) -> float:
        """Compute sum of energy field using Taichi kernel"""
        total = 0.0
        for i, j in ti.ndrange(self.energy_manager.energy_system.energy_field.shape[0], 
                              self.energy_manager.energy_system.energy_field.shape[1]):
            total += self.energy_manager.energy_system.energy_field[i, j]
        return total
    
    @ti.kernel
    def _compute_kinetic_energy(self) -> float:
        """Compute kinetic energy using Taichi kernel"""
        total = 0.0
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                mass = self.particles.attributes[i][0]
                vx, vy = self.particles.velocities[i][0], self.particles.velocities[i][1]
                kinetic = 0.5 * mass * (vx * vx + vy * vy)
                total += kinetic
        return total
    
    @ti.kernel
    def _compute_max_velocity(self) -> float:
        """Compute maximum particle velocity using Taichi kernel"""
        max_vel = 0.0
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                vx, vy = self.particles.velocities[i][0], self.particles.velocities[i][1]
                vel_mag = ti.sqrt(vx * vx + vy * vy)
                if vel_mag > max_vel:
                    max_vel = vel_mag
        return max_vel
    
    def _get_total_energy(self):
        """Calculate total energy of the system"""
        total_energy = 0.0
        
        # MEMORY LEAK FIX: Use Taichi kernel instead of to_numpy()
        # Energy field energy
        total_energy += self._compute_energy_field_sum()
        
        # Particle kinetic energy
        if self.config.mode == "open_chemistry":
            total_energy += self._compute_kinetic_energy()
        
        return total_energy
    
    def _get_max_particle_velocity(self):
        """Get maximum particle velocity for CFL condition - MEMORY LEAK FIX"""
        if self.config.mode != "open_chemistry":
            return 0.0
        
        # MEMORY LEAK FIX: Use Taichi kernel instead of to_numpy()
        return self._compute_max_velocity()
    
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
    
    def _symplectic_step(self, dt: float) -> tuple[bool, float]:
        """Perform symplectic integration step with error control"""
        particle_count = self.particles.particle_count[None]
        
        try:
            # Perform adaptive symplectic step
            actual_dt, success = self.integrator.step(
                self.particles.positions,
                self.particles.velocities,
                self.potentials.forces,
                self.particles.attributes,
                dt,
                self.particles.active,
                particle_count
            )
            
            if success and actual_dt > 0:
                # Successful step - recompute forces with new positions
                self.potentials.compute_forces(
                    self.particles.positions,
                    self.particles.attributes,
                    self.particles.active,
                    particle_count
                )
                
                # Complete Verlet correction if needed
                self.integrator.integrator.verlet_correction(
                    self.particles.velocities,
                    self.potentials.forces,
                    self.particles.attributes,
                    actual_dt,
                    self.particles.active,
                    particle_count
                )
                
                return True, actual_dt
            else:
                # Failed step
                return False, 0.0
                
        except Exception as e:
            logger.warning(f"Symplectic integration failed: {e}")
            return False, 0.0

    def get_integration_stats(self) -> Dict[str, Any]:
        """Get integration quality statistics"""
        integration_stats = self.integrator.get_stats()
        
        # Add energy conservation metrics
        total_energy = self._get_total_energy()
        if hasattr(self, 'initial_energy') and self.initial_energy > 0:
            energy_relative_error = abs(total_energy - self.initial_energy) / self.initial_energy
            integration_stats['energy_conservation'] = {
                'total_energy': total_energy,
                'initial_energy': self.initial_energy,
                'relative_error': energy_relative_error
            }
        
        return integration_stats
    
    def _update_energy_conservation(self):
        """Update energy conservation history"""
        current_energy = self._get_total_energy()
        
        # Initialize initial energy on first call
        if self.initial_energy == 0.0 and current_energy > 0.0:
            self.initial_energy = current_energy
        
        # Add to history (deque automatically handles max length)
        self.energy_history.append(current_energy)
        
        # Check for energy drift violation
        # Only warn for drift in closed systems or very large drift in open systems
        drift = self._get_energy_drift()
        
        # Use different thresholds for open vs closed systems
        if self.config.mode == "open_chemistry":
            # For open systems, only warn if short-term drift is very high (10%)
            threshold = 10.0
        else:
            # For closed systems, use configured threshold (5%)
            threshold = self.energy_conservation_threshold * 100.0
        
        if drift > threshold:
            # Log warning only every 1000 steps to avoid spam
            if self.step_count % 1000 == 0:
                logger.warning(f"Energy drift ({drift:.2f}%) exceeds threshold ({threshold:.2f}%)")
    
    def _get_energy_drift(self):
        """
        Calculate energy drift percentage.
        For open systems (with energy inputs), this measures drift from recent average.
        For closed systems, it measures drift from initial energy.
        """
        current_energy = self._get_total_energy()
        
        # For open chemistry mode, calculate drift from recent average (last 100 steps)
        # instead of from initial energy, since energy is continuously added
        if self.config.mode == "open_chemistry" and len(self.energy_history) > 100:
            # Use average of last 100 steps as baseline
            recent_avg = sum(self.energy_history[-100:]) / 100.0
            if recent_avg > 0:
                drift = abs(current_energy - recent_avg) / recent_avg * 100.0
                return drift
            return 0.0
        
        # For closed systems (preset mode), use initial energy
        if self.initial_energy == 0.0:
            return 0.0
        
        drift = abs(current_energy - self.initial_energy) / self.initial_energy * 100.0
        return drift
    
    def _transfer_energy_to_particles(self):
        """Transfer energy from field to particles without duplication"""
        self._transfer_energy_kernel()
    
    @ti.kernel
    def _transfer_energy_kernel(self):
        """Transfer energy from field to particles (no duplication)"""
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                pos = self.particles.positions[i]
                x = int(pos[0]) % self.energy_manager.energy_system.width
                y = int(pos[1]) % self.energy_manager.energy_system.height
                
                # Get energy from field
                field_energy = self.energy_manager.energy_system.energy_field[x, y]
                
                # Transfer small amount to particle (for bonding)
                transfer_amount = field_energy * 0.05  # Increased from 1% to 5% for better bonding
                
                # Remove from field and add to particle
                self.energy_manager.energy_system.energy_field[x, y] -= transfer_amount
                self.particles.energy[i] += transfer_amount
    
    @ti.kernel
    def _add_energy_to_particles_kernel(self):
        """Vectorized kernel to add energy from field to all particles"""
        total_energy_available = 0.0
        particles_with_energy = 0
        
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                pos = self.particles.positions[i]
                # Direct field access instead of function call
                x = int(pos[0]) % self.energy_manager.energy_system.width
                y = int(pos[1]) % self.energy_manager.energy_system.height
                energy = self.energy_manager.energy_system.energy_field[x, y]
                self.energy_amounts[i] = energy
                
                if energy > 0.0:
                    total_energy_available += energy
                    particles_with_energy += 1
            else:
                self.energy_amounts[i] = 0.0
        
        # Store debug info in a field that can be read outside kernel
        self.debug_energy_total[None] = total_energy_available
        self.debug_particles_with_energy[None] = particles_with_energy
    
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
        
        # Log mutations only every 200 steps
        if self.step_count % 200 == 0:
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
        
        # Debug logging
        if self.step_count % 500 == 0:  # Log every 500 steps
            logger.info(f"detect_novel_substances: Found {len(clusters)} clusters (min_size={self.config.min_cluster_size})")
            for i, cluster in enumerate(clusters[:3]):  # Show first 3 clusters
                logger.info(f"  Cluster {i}: size={len(cluster)}, particles={cluster[:5]}...")
        
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
                
                # Calculate properties for the substance
                properties = self._calculate_substance_properties(particle_attributes)
                
                # Add to catalog
                is_novel, substance_id = self.catalog.add_substance(graph, self.current_time, properties)
                
                # Update novelty tracker
                self.novelty_tracker.record_discovery(is_novel, complexity)
    
    def _calculate_substance_properties(self, particle_attributes: Dict[int, np.ndarray]) -> Dict[str, Any]:
        """Calculate properties for a substance based on particle attributes"""
        if not particle_attributes:
            return {}
        
        # Convert to numpy array for easier processing
        attributes_list = list(particle_attributes.values())
        attributes_array = np.array(attributes_list)
        
        # Calculate average mass (first component of attributes)
        avg_mass = float(np.mean(attributes_array[:, 0]))
        
        # Calculate average charge (components 1-3 of attributes)
        avg_charge = [float(x) for x in np.mean(attributes_array[:, 1:4], axis=0)]
        
        # Calculate other properties
        mass_std = float(np.std(attributes_array[:, 0]))
        charge_magnitude = float(np.linalg.norm(avg_charge))
        
        return {
            'avg_mass': avg_mass,
            'avg_charge': avg_charge,
            'mass_std': mass_std,
            'charge_magnitude': charge_magnitude,
            'particle_count': len(particle_attributes)
        }
    
    def calculate_complexity(self, graph: MolecularGraph, 
                           particle_attributes: Dict[int, np.ndarray]) -> float:
        """Calculate complexity of a molecular structure"""
        from .metrics import ComplexityAnalyzer
        
        # Convert dict to list for ComplexityAnalyzer
        attributes_list = list(particle_attributes.values())
        return ComplexityAnalyzer.calculate_emergent_complexity(graph, attributes_list)
    
    def update_metrics(self):
        """Update all metrics"""
        import time
        start_time = time.time()
        
        # Update particle metrics (includes mass, stored energy, and kinetic energy)
        particle_count = self.particles.particle_count[None]
        # logger.info(f"DEBUG: update_metrics called, particle_count from ParticleSystem={particle_count}")
        
        particle_start = time.time()
        self.metrics.update_particle_metrics(
            self.particles.active,
            self.particles.attributes,
            self.particles.energy,
            self.particles.velocities,
            particle_count
        )
        particle_time = time.time() - particle_start
        # Debug print removed for performance
        
        # Update bond metrics manually - OPTIMIZED with NumPy
        bond_count = 0
        try:
            bond_start = time.time()
            # MEMORY FIX: Limit to 1000 particles (was 500) to improve cluster detection
            max_check = min(particle_count, 1000)
            if max_check > 0:
                bond_matrix = self.binding.bond_matrix.to_numpy()[:max_check, :max_check]
                bond_active = self.binding.bond_active.to_numpy()[:max_check, :max_check]
                
                # Count active bonds using NumPy (much faster than Python loops)
                import numpy as np
                bond_count = int(np.sum(np.triu(bond_active, k=1)))
                self.metrics.bond_count[None] = bond_count
            else:
                self.metrics.bond_count[None] = 0
            bond_time = time.time() - bond_start
            # Debug print removed for performance
            # logger.info(f"DEBUG: Manual bond calculation - bonds={bond_count}")
        except Exception as e:
            logger.error(f"ERROR in bond calculation: {e}")
            self.metrics.bond_count[None] = 0
        
        # Update cluster metrics manually - OPTIMIZED with simplified algorithm
        cluster_count = 0
        try:
            cluster_start = time.time()
            # MEMORY FIX: Limit to 1000 particles (was 500) to improve cluster detection
            max_check = min(particle_count, 1000)
            
            # Safe array access with bounds checking
            if max_check > 0:
                bond_active = self.binding.bond_active.to_numpy()[:max_check, :max_check]
                particles_active = self.particles.active.to_numpy()[:max_check]
                
                # Count particles that have at least one bond
                import numpy as np
                particles_with_bonds = np.sum(particles_active * np.any(bond_active, axis=1))
                
                # Simple approximation: particles_with_bonds / 2 (assuming average 2 particles per cluster)
                cluster_count = max(1, int(particles_with_bonds / 2))
            else:
                cluster_count = 1
            
            self.metrics.cluster_count[None] = cluster_count
            cluster_time = time.time() - cluster_start
            # Debug print removed for performance
            # logger.info(f"DEBUG: Manual cluster calculation - clusters={cluster_count}")
        except Exception as e:
            logger.error(f"ERROR in cluster calculation: {e}")
            # Fallback to simple count
            self.metrics.cluster_count[None] = max(1, min(particle_count, 100))
        
        # Update energy field metrics
        energy_start = time.time()
        energy_field = self.energy_manager.energy_system.energy_field
        self.metrics.update_energy_field_metrics(
            energy_field,
            self.config.grid_width,
            self.config.grid_height
        )
        energy_time = time.time() - energy_start
        # Debug print removed for performance
        
        # Record metrics
        record_start = time.time()
        additional_metrics = {
            'novelty_rate': self.novelty_tracker.get_novelty_rate(),
            'discovery_rate': self.novelty_tracker.get_discovery_rate(),
            'total_substances': len(self.catalog.substances),
            'health_score': self.aggregator.get_health_score(),
            'simulation_time': self.current_time
        }
        
        self.metrics.record_metrics(additional_metrics)
        self.aggregator.update_aggregated_stats()
        record_time = time.time() - record_start
        # Debug print removed for performance
        
        total_time = time.time() - start_time
        # Debug print removed for performance
    
    def _log_diagnostics(self):
        """Log diagnostics data for current step - OPTIMIZED to reduce memory usage"""
        if not self.diagnostics.enabled:
            return
        
        try:
            # PERFORMANCE FIX: Sample particles to reduce memory usage
            particle_count = self.particles.particle_count[None]
            
            # MEMORY FIX: Limit diagnostics to max 200 particles to prevent memory issues
            max_diag_particles = min(200, particle_count)
            
            # Prepare simulation data for diagnostics
            simulation_data = {}
            
            # Get particle positions (sampled)
            if particle_count > 0:
                if particle_count > max_diag_particles:
                    # Sample particles uniformly
                    indices = np.linspace(0, particle_count-1, max_diag_particles, dtype=int)
                    positions = self.particles.positions.to_numpy()[indices]
                    attributes = self.particles.attributes.to_numpy()[indices]
                    particle_energies = self.particles.energy.to_numpy()[indices]
                else:
                    positions = self.particles.positions.to_numpy()[:particle_count]
                    attributes = self.particles.attributes.to_numpy()[:particle_count]
                    particle_energies = self.particles.energy.to_numpy()[:particle_count]
                
                simulation_data['positions'] = positions
                simulation_data['attributes'] = attributes
                simulation_data['particle_energies'] = particle_energies
            else:
                simulation_data['positions'] = np.array([])
                simulation_data['attributes'] = np.array([])
                simulation_data['particle_energies'] = np.array([])
            
            # Get bonds (sampled to prevent O(n²) operation)
            bonds = []
            if max_diag_particles > 0:
                bond_matrix = self.binding.bond_matrix.to_numpy()[:max_diag_particles, :max_diag_particles]
                for i in range(max_diag_particles):
                    for j in range(i + 1, max_diag_particles):
                        strength = bond_matrix[i, j]
                        if strength > 0:
                            bonds.append((i, j, strength))
            simulation_data['bonds'] = bonds
            
            # Get clusters (simplified) - FIX: Use empty dict instead of note to avoid str/int comparison error
            simulation_data['clusters'] = {}
            
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
        self.novelty_tracker = NoveltyTracker(window_size=self.config.novelty_window)
        self.energy_manager = EnergyManager(self.config)
        self.rng.reset(self.config.seed)
        
        # Reinitialize
        self.initialize_simulation()
        
        # Start simulation automatically after reset
        self.start()
    
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
        """Get data for visualization - SIMPLIFIED VERSION for speed"""
        # Start visualization timing
        self.performance_monitor.start_visualization_timing()
        
        import time
        t_start = time.time()
        
        # Time metrics collection
        t_metrics_start = time.time()
        metrics = self.aggregator.get_aggregated_stats()
        t_metrics_end = time.time()
        
        # Time performance metrics
        t_perf_start = time.time()
        performance_metrics = self.performance_monitor.get_performance_metrics(self.step_count)
        t_perf_end = time.time()

        data = {
            'metrics': metrics,
            'step_count': self.step_count,  # Add step count for frontend
            'current_time': self.current_time,  # Add current time for frontend
            'performance': performance_metrics
        }

        # Debug: Log step_count and current_time
        if self.step_count < 5:
            # Debug print removed for performance
            pass
        
        # Time energy field extraction - OPTIMIZED with caching
        t_energy_start = time.time()
        # Cache energy field to reduce expensive GPU->CPU transfers
        if self.step_count - self._energy_field_cache_step >= self._energy_field_cache_interval:
            self._energy_field_cache = self.energy_manager.energy_system.energy_field.to_numpy()
            self._energy_field_cache_step = self.step_count
        data['energy_field'] = self._energy_field_cache
        t_energy_end = time.time()
        
        if self.config.mode == "preset_prebiotic":
            # Provide concentration fields for preset mode
            if self.preset_simulator is not None:
                # Extract only one channel (selected species index 0 by default) to reduce bandwidth
                try:
                    conc = self.preset_simulator.concentration_fields.get_all_concentrations()
                    # pick first available species
                    first_key = next(iter(conc.keys())) if conc else None
                    if first_key is not None:
                        # OPTIMIZATION: Return NumPy array
                        data['concentration_view'] = conc[first_key]
                        data['concentration_species'] = list(conc.keys())
                        # Back-compat: also provide single-key concentrations for older clients
                        data['concentrations'] = {first_key: conc[first_key]}
                except Exception:
                    pass
        else:
            # Provide particle/bond/cluster data for open chemistry - OPTIMIZED with caching
            # Time particles extraction
            t_particles_start = time.time()
            # OPTIMIZATION: Only get particles every 10 steps to reduce load - PERFORMANCE FIX
            if self.step_count % 10 == 0:
                positions, velocities, attributes, active_mask, energies = self.particles.get_active_particles()
                # Cache the data for intermediate steps
                self._cached_particles_data = {
                    'positions': positions,
                    'attributes': attributes,
                    'active_mask': active_mask,
                    'energies': energies
                }
            else:
                # Use cached data for intermediate steps
                if hasattr(self, '_cached_particles_data'):
                    positions = self._cached_particles_data['positions']
                    attributes = self._cached_particles_data['attributes']
                    active_mask = self._cached_particles_data['active_mask']
                    energies = self._cached_particles_data['energies']
                else:
                    # Fallback: get fresh data if cache doesn't exist
                    positions, velocities, attributes, active_mask, energies = self.particles.get_active_particles()
            t_particles_end = time.time()
            
            # Time bonds/clusters extraction
            t_bonds_start = time.time()
            # OPTIMIZATION: Only get bonds/clusters every 50 steps to reduce load (faster than 200) - PERFORMANCE FIX
            if self.step_count % 50 == 0:
                bonds = self.binding.get_bonds()
                clusters = self.binding.get_clusters()
                # Cache the data for intermediate steps
                self._cached_bonds_clusters_data = {
                    'bonds': bonds,
                    'clusters': clusters
                }
            else:
                # Use cached data for intermediate steps
                if hasattr(self, '_cached_bonds_clusters_data') and self._cached_bonds_clusters_data is not None:
                    bonds = self._cached_bonds_clusters_data['bonds']
                    clusters = self._cached_bonds_clusters_data['clusters']
                else:
                    # Fallback: get fresh data if cache doesn't exist
                    bonds = self.binding.get_bonds()
                    clusters = self.binding.get_clusters()
                    # Initialize cache for next time
                    self._cached_bonds_clusters_data = {
                        'bonds': bonds,
                        'clusters': clusters
                    }
            t_bonds_end = time.time()
            
            # OPTIMIZATION: Return NumPy arrays instead of lists
            data['particles'] = {
                'positions': positions,
                'attributes': attributes,
                'active_mask': active_mask,
                'energies': energies
            }
            data['bonds'] = bonds
            data['clusters'] = clusters
        
        t_end = time.time()
        if (t_end - t_start) > 1.0:
            logger.warning(f"Visualization data took {t_end-t_start:.2f}s")
        
        # Log detailed timing breakdown
        if (t_end - t_start) > 0.1:  # Log if total takes more than 100ms
            logger.warning(f"Slow visualization breakdown:")
            logger.warning(f"  Metrics: {(t_metrics_end-t_metrics_start)*1000:.1f}ms")
            logger.warning(f"  Performance: {(t_perf_end-t_perf_start)*1000:.1f}ms")
            logger.warning(f"  Energy field: {(t_energy_end-t_energy_start)*1000:.1f}ms")
            if 't_particles_end' in locals():
                logger.warning(f"  Particles: {(t_particles_end-t_particles_start)*1000:.1f}ms")
            if 't_bonds_end' in locals():
                logger.warning(f"  Bonds/Clusters: {(t_bonds_end-t_bonds_start)*1000:.1f}ms")
            logger.warning(f"  Total: {(t_end-t_start)*1000:.1f}ms")
        
        # End visualization timing
        viz_time = self.performance_monitor.end_visualization_timing()
        if viz_time > 0.1:  # Log if visualization takes more than 100ms
            logger.warning(f"Slow visualization: {viz_time*1000:.1f}ms")
        
        return data
    
    def get_thermodynamic_validation_summary(self) -> Dict[str, Any]:
        """Get summary of thermodynamic validation results"""
        return self.validator.get_validation_summary()
    
    def get_validation_log(self) -> List[Dict[str, Any]]:
        """Get full validation log"""
        return self.validator.validation_log
    
    @ti.kernel
    def _copy_energy_field_to_taichi(self):
        """Copy energy field to Taichi field using kernel"""
        for i, j in ti.ndrange(self.config.grid_height, self.config.grid_width):
            self._energy_field_taichi[i, j] = self.energy_manager.energy_system.energy_field[i, j]
    
    def get_novel_substances(self, count: int = 10) -> List[Dict]:
        """Get recent novel substances"""
        recent_substances = self.catalog.get_recent_substances(count)
        
        return [substance.to_dict() for substance in recent_substances]
    
    def save_snapshot(self, filename: str, save_images: bool = True):
        """Save simulation snapshot with optional image generation"""
        import os
        # Debug print removed for performance
        # Debug print removed for performance
        
        try:
            # Serialize simulation data
            # Debug print removed for performance
            snapshot_data = SnapshotSerializer.serialize_simulation(self)
            # Debug print removed for performance
            
            # Use SnapshotManager to save with image generation
            name = filename.replace('.json', '')  # Remove .json extension for name
            # Debug print removed for performance
            
            saved_filename = self.snapshot_manager.create_snapshot(
                snapshot_data, 
                name=name, 
                save_images=save_images
            )
            logger.info(f"Snapshot saved: {saved_filename}")
            return saved_filename
            
        except Exception as e:
            logger.error(f"Snapshot save failed: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    @ti.kernel
    def _energy_diffuse(self, dt: float):
        """Energy diffusion kernel"""
        D = self.config.diffuse_D
        H, W = self.energy_manager.energy_system.energy_field.shape
        
        for i, j in ti.ndrange(H, W):
            ip = (i + 1) % H
            im = (i - 1 + H) % H
            jp = (j + 1) % W
            jm = (j - 1 + W) % W
            
            # Laplacian
            lap = (self.energy_manager.energy_system.energy_field[ip, j] + 
                   self.energy_manager.energy_system.energy_field[im, j] +
                   self.energy_manager.energy_system.energy_field[i, jp] + 
                   self.energy_manager.energy_system.energy_field[i, jm] -
                   4.0 * self.energy_manager.energy_system.energy_field[i, j])
            
            self.energy_manager.energy_system.energy_field[i, j] += D * lap * dt
    
    def _add_energy_pulse(self):
        """Add periodic energy pulses to random locations"""
        import random
        
        # Get pulse parameters from config
        pulse_radius = getattr(self.config, 'pulse_radius', 24.0)
        pulse_amplitude = getattr(self.config, 'pulse_amplitude', 5.0)
        
        # Random position in the grid
        x = random.uniform(0, self.config.grid_width)
        y = random.uniform(0, self.config.grid_height)
        
        # Add energy impulse
        self.energy_manager.add_energy_impulse(
            position=(x, y),
            intensity=pulse_amplitude,
            radius=pulse_radius
        )
    
    @ti.kernel
    def _energy_thermostat(self):
        """Energy thermostat kernel - controls both field and particle kinetic energy"""
        target = self.config.target_energy
        alpha = self.config.thermostat_alpha
        H, W = self.energy_manager.energy_system.energy_field.shape
        
        # Control energy field
        for i, j in ti.ndrange(H, W):
            e = self.energy_manager.energy_system.energy_field[i, j]
            # Noisy pull to target
            self.energy_manager.energy_system.energy_field[i, j] = e + alpha * (target - e) + (ti.random(ti.f32) - 0.5) * 0.005
        
        # Control particle kinetic energy - CRITICAL FIX for high energy
        kinetic_target = target * 0.02  # Kinetic energy should be ~2% of field energy (FURTHER REDUCED from 5%)
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                # Calculate current kinetic energy
                mass = self.particles.attributes[i][0]
                vx = self.particles.velocities[i][0]
                vy = self.particles.velocities[i][1]
                current_kinetic = 0.5 * mass * (vx*vx + vy*vy)
                
                # Apply thermostat to kinetic energy - MORE AGGRESSIVE
                if current_kinetic > kinetic_target:
                    # Scale down velocity to match target kinetic energy
                    scale_factor = ti.sqrt(kinetic_target / ti.max(current_kinetic, 1e-6))
                    self.particles.velocities[i][0] *= scale_factor
                    self.particles.velocities[i][1] *= scale_factor
                
                # Additional velocity damping for bond formation - MODERATE
                damping_factor = 0.98  # Reduce velocity by 2% each step (was 15% - too aggressive)
                self.particles.velocities[i][0] *= damping_factor
                self.particles.velocities[i][1] *= damping_factor
    
    def load_snapshot(self, filename: str):
        """Load simulation snapshot using serializer"""
        import json
        with open(filename, 'r') as f:
            snapshot = json.load(f)
        if not SnapshotSerializer.validate_snapshot(snapshot):
            raise ValueError("Invalid snapshot data")
        SnapshotSerializer.deserialize_simulation(snapshot, self)
    
    @ti.kernel
    def _stabilize_particle_positions(self):
        """Additional stabilization to prevent particles from drifting out of grid"""
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                pos = self.particles.positions[i]
                
                # Check if particle is outside grid bounds
                if pos[0] < 0 or pos[0] >= self.config.grid_width or \
                   pos[1] < 0 or pos[1] >= self.config.grid_height:
                    
                    # Reset position to center of grid if particle escaped
                    self.particles.positions[i] = ti.Vector([
                        self.config.grid_width * 0.5,
                        self.config.grid_height * 0.5
                    ])
                    
                    # Reset velocity to prevent immediate re-escape
                    self.particles.velocities[i] = ti.Vector([0.0, 0.0])
    
    @ti.kernel
    def _assist_clustering(self):
        """Assist clustering by moving particles towards high energy regions"""
        # Find high energy regions and move nearby particles towards them
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                pos = self.particles.positions[i]
                
                # Get local energy
                x = int(pos[0]) % self.energy_manager.energy_system.width
                y = int(pos[1]) % self.energy_manager.energy_system.height
                local_energy = self.energy_manager.energy_system.energy_field[x, y]
                
                # If in high energy region, add small attractive force towards center
                if local_energy > 0.5:  # High energy threshold
                    # Calculate direction towards grid center
                    center_x = self.config.grid_width * 0.5
                    center_y = self.config.grid_height * 0.5
                    
                    dx = center_x - pos[0]
                    dy = center_y - pos[1]
                    dist = ti.sqrt(dx * dx + dy * dy)
                    
                    if dist > 0.1:  # Avoid division by zero
                        # Normalize direction
                        dx /= dist
                        dy /= dist
                        
                        # Add small attractive velocity
                        attraction_strength = 0.1 * local_energy
                        self.particles.velocities[i].x += dx * attraction_strength
                        self.particles.velocities[i].y += dy * attraction_strength
    
    @ti.kernel
    def _attract_particles_for_bonding(self):
        """Attract particles to each other to facilitate bonding"""
        # Find nearby particles and add attractive forces
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                pos_i = self.particles.positions[i]
                
                # Look for nearby particles
                for j in range(i + 1, self.config.max_particles):
                    if self.particles.active[j] == 1:
                        pos_j = self.particles.positions[j]
                        
                        # Calculate distance
                        dx = pos_j[0] - pos_i[0]
                        dy = pos_j[1] - pos_i[1]
                        dist = ti.sqrt(dx * dx + dy * dy)
                        
                        # If particles are close but not too close, add attraction
                        if dist > 0.5 and dist < 5.0:  # Attraction range
                            # Normalize direction
                            dx /= dist
                            dy /= dist
                            
                            # Add attractive velocity (particles move towards each other)
                            attraction_strength = 0.05  # Small attraction force
                            self.particles.velocities[i].x += dx * attraction_strength
                            self.particles.velocities[i].y += dy * attraction_strength
                            self.particles.velocities[j].x -= dx * attraction_strength
                            self.particles.velocities[j].y -= dy * attraction_strength
    
    @ti.kernel
    def _force_clustering_to_center(self):
        """Force all particles to move towards grid center for maximum clustering"""
        center_x = self.config.grid_width * 0.5
        center_y = self.config.grid_height * 0.5
        
        for i in range(self.config.max_particles):
            if self.particles.active[i] == 1:
                pos = self.particles.positions[i]
                
                # Calculate direction towards center
                dx = center_x - pos[0]
                dy = center_y - pos[1]
                dist = ti.sqrt(dx * dx + dy * dy)
                
                if dist > 0.1:  # Avoid division by zero
                    # Normalize direction
                    dx /= dist
                    dy /= dist
                    
                    # Add strong attractive velocity towards center
                    center_attraction = 0.2  # Strong attraction force
                    self.particles.velocities[i].x += dx * center_attraction
                    self.particles.velocities[i].y += dy * center_attraction
