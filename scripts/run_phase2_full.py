"""
Phase 2 Full Simulation Runner
===============================

Complete Phase 2 runner using the full Live 2.0 simulation system.
Integrates Phase2Config with SimulationStepper.

Usage:
    python scripts/run_phase2_full.py --config configs/phase2_miller_urey.yaml \\
                                      --output results/miller_urey/run01 \\
                                      --steps 10000000 \\
                                      --seed 42
"""

import sys
import argparse
import time
import logging
import json
import multiprocessing
from pathlib import Path
from datetime import datetime
from typing import List, Dict

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import taichi as ti
import numpy as np

from backend.sim.config import SimulationConfig
from backend.sim.phase2_config import Phase2Config, load_phase2_config_from_yaml
from backend.sim.phase2_initializer import initialize_phase2_simulation
from backend.sim.core.stepper import SimulationStepper


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Phase2FullRunner:
    """Full Phase 2 simulation runner"""
    
    def __init__(self,
                 config_path: str,
                 output_dir: str,
                 max_steps: int = 10000000,
                 seed: int = 42,
                 force_cpu: bool = False):
        """
        Initialize Phase 2 runner
        
        Args:
            config_path: Path to Phase 2 YAML configuration
            output_dir: Output directory for results
            max_steps: Maximum simulation steps
            seed: Random seed
            force_cpu: Force CPU mode (skip GPU detection)
        """
        self.config_path = Path(config_path)
        self.output_dir = Path(output_dir)
        self.max_steps = max_steps
        self.seed = seed
        self.force_cpu = force_cpu
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup file logging
        log_file = self.output_dir / "simulation.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        )
        logging.getLogger().addHandler(file_handler)
        
        logger.info("=" * 70)
        logger.info("PHASE 2 FULL SIMULATION RUNNER")
        logger.info("=" * 70)
        logger.info(f"Config: {self.config_path}")
        logger.info(f"Output: {self.output_dir}")
        logger.info(f"Steps: {self.max_steps:,}")
        logger.info(f"Seed: {self.seed}")
        logger.info("=" * 70)
    
    @staticmethod
    def _json_serializer(obj):
        """Custom JSON serializer for NumPy types"""
        if isinstance(obj, (np.integer, np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
    
    def load_configs(self) -> tuple:
        """
        Load Phase 2 and Simulation configurations
        
        Returns:
            (sim_config, phase2_config)
        """
        logger.info("Loading configurations...")
        
        # Load Phase 2 config from YAML
        phase2_config = load_phase2_config_from_yaml(str(self.config_path))
        
        logger.info(f"  Scenario: {phase2_config.scenario_name}")
        logger.info(f"  Description: {phase2_config.description}")
        logger.info(f"  Molecules: {len(phase2_config.initial_molecules)}")
        logger.info(f"  Temperature: {phase2_config.temperature}K")
        logger.info(f"  Energy injection: {phase2_config.energy_injection.enabled}")
        logger.info(f"  Catalysts: {len(phase2_config.catalysts)}")
        
        # Create SimulationConfig
        # Calculate total particle count
        total_atoms = sum(
            self._count_atoms(mol.formula) * mol.count
            for mol in phase2_config.initial_molecules
        )
        
        # Add catalyst particles
        catalyst_particles = sum(
            int(10000 * cat.concentration)
            for cat in phase2_config.catalysts
        )
        
        max_particles = total_atoms + catalyst_particles + 1000  # Buffer
        max_particles = min(max_particles, 50000)  # Cap at 50k
        
        logger.info(f"  Estimated particles: {total_atoms} atoms + {catalyst_particles} catalysts")
        logger.info(f"  Max particles: {max_particles}")
        
        # Create simulation config
        # Use n_particles from phase2_config if available, otherwise use calculated max_particles
        config_n_particles = getattr(phase2_config, 'n_particles', max_particles)
        if config_n_particles > 0:
            max_particles = max(config_n_particles, max_particles)
        
        # Get novelty detection settings
        detect_enabled = getattr(phase2_config, 'detect_novel_substances', True)
        if detect_enabled is None:
            detect_enabled = True  # Default to enabled
        novelty_interval = getattr(phase2_config, 'novelty_check_interval', 500)
        if novelty_interval is None:
            novelty_interval = 500
        
        if not detect_enabled:
            logger.info("FAST MODE: Novelty detection DISABLED during simulation")
            novelty_interval = 99999999  # Never run
        
        # Get mutations control
        enable_mutations_val = getattr(phase2_config, 'enable_mutations', True)
        if enable_mutations_val is None:
            enable_mutations_val = True  # Default to enabled
        
        # Get grid size from config or use box_size as fallback
        # Try multiple fallbacks in order
        grid_width = getattr(phase2_config, 'grid_width', None)
        if grid_width is None:
            grid_width = getattr(phase2_config, 'box_size', 256)
        
        grid_height = getattr(phase2_config, 'grid_height', None)
        if grid_height is None:
            grid_height = getattr(phase2_config, 'box_size', 256)
        
        # Convert to int and ensure valid
        grid_width = int(grid_width)
        grid_height = int(grid_height)
        
        logger.info(f"Grid size: {grid_width}x{grid_height}")
        
        sim_config = SimulationConfig(
            mode='open_chemistry',
            grid_width=grid_width,
            grid_height=grid_height,
            max_particles=max_particles,
            dt=getattr(phase2_config, 'dt', 0.005),
            max_time=float(self.max_steps) * getattr(phase2_config, 'dt', 0.005),
            temperature=phase2_config.temperature,
            seed=self.seed,
            
            # Use validation settings from phase2_config or defaults
            enable_thermodynamic_validation=phase2_config.enable_thermodynamic_validation if phase2_config.enable_thermodynamic_validation is not None else False,
            validate_every_n_steps=phase2_config.validate_every_n_steps if phase2_config.validate_every_n_steps is not None else 10000,
            
            # Performance settings from phase2_config
            energy_update_interval=phase2_config.energy_update_interval if phase2_config.energy_update_interval is not None else 5,
            metrics_update_interval=phase2_config.metrics_update_interval if phase2_config.metrics_update_interval is not None else 1,
            
            # Enable diagnostics
            enable_diagnostics=phase2_config.enable_diagnostics if phase2_config.enable_diagnostics is not None else True,
            diagnostics_frequency=phase2_config.diagnostics_frequency if phase2_config.diagnostics_frequency is not None else 100,
            
            # Novelty detection control
            detect_novel_substances=detect_enabled,
            novelty_check_interval=novelty_interval,
            
            # Mutations control (disable on CPU to avoid LLVM errors)
            enable_mutations=enable_mutations_val,
            
            # Use physics database
            use_physics_db=True,
        )
        
        logger.info("Configurations loaded successfully!")
        
        return sim_config, phase2_config
    
    def _count_atoms(self, formula: str) -> int:
        """Count atoms in formula (simple parser)"""
        count = 0
        i = 0
        while i < len(formula):
            if formula[i].isupper():
                # Found atom
                i += 1
                if i < len(formula) and formula[i].islower():
                    i += 1
                
                # Count
                num_str = ''
                while i < len(formula) and formula[i].isdigit():
                    num_str += formula[i]
                    i += 1
                
                count += int(num_str) if num_str else 1
            else:
                i += 1
        
        return count
    
    def initialize_taichi(self):
        """Initialize Taichi - Try GPU first, fallback to CPU with all available cores"""
        logger.info("Initializing Taichi...")
        
        # Get available CPU cores
        available_cores = multiprocessing.cpu_count()
        logger.info(f"Available CPU cores: {available_cores}")
        
        # Force CPU mode if requested (useful for AWS instances without GPU)
        if self.force_cpu:
            logger.info("Force CPU mode enabled - skipping GPU detection")
            ti.init(arch=ti.cpu, cpu_max_num_threads=available_cores)
            logger.info(f"✅ Using CPU with ALL {available_cores} threads")
            logger.info(f"   This will utilize all {available_cores} CPU cores for maximum performance")
            return
        
        try:
            # Try GPU first
            ti.init(arch=ti.cuda)
            
            # Quick test to verify GPU is actually being used
            @ti.kernel
            def quick_test(arr: ti.template()):
                for i in range(arr.shape[0]):
                    arr[i] = i * 2.0
            
            test_arr = ti.field(dtype=ti.f32, shape=1000)
            start = time.time()
            quick_test(test_arr)
            ti.sync()
            elapsed = time.time() - start
            
            if elapsed > 0.1:  # More than 100ms = likely not using GPU
                logger.warning(f"GPU initialization slow ({elapsed*1000:.0f}ms) - GPU may be busy!")
                logger.warning("Falling back to CPU with all available cores...")
                ti.reset()
                ti.init(arch=ti.cpu, cpu_max_num_threads=available_cores)
                logger.info(f"✅ Using CPU with ALL {available_cores} threads")
            else:
                logger.info("✅ Using GPU acceleration")
                
        except Exception as e:
            logger.warning(f"GPU not available: {e}")
            logger.info("Falling back to CPU with all available cores...")
            max_threads = available_cores
            ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)
            logger.info(f"✅ Using CPU with ALL {max_threads} threads")
            logger.info(f"   This will utilize all {max_threads} CPU cores for maximum performance")
    
    def run(self) -> dict:
        """
        Run complete Phase 2 simulation
        
        Returns:
            Results dictionary
        """
        start_time = time.time()
        
        # Load configurations
        sim_config, phase2_config = self.load_configs()
        
        # Initialize Taichi
        self.initialize_taichi()
        
        # Create simulation stepper
        logger.info("Creating simulation stepper...")
        stepper = SimulationStepper(sim_config)
        
        # Initialize Phase 2 molecules
        logger.info("Initializing Phase 2 molecules...")
        init_stats = initialize_phase2_simulation(
            sim_config,
            phase2_config,
            stepper.particles,  # Fixed: it's 'particles' not 'particle_system'
            stepper
        )
        
        # Run simulation
        logger.info("=" * 70)
        logger.info(f"STARTING SIMULATION - {self.max_steps:,} STEPS")
        logger.info("=" * 70)
        
        # CRITICAL: Start the simulation
        logger.info("Starting simulation stepper...")
        stepper.start()
        logger.info(f"Stepper started - is_running={stepper.is_running}")
        
        save_interval = max(self.max_steps // 20, 50000)  # Save 20 snapshots
        checkpoint_interval = max(self.max_steps // 10, 100000)  # 10 checkpoints
        
        for step in range(self.max_steps):
            # Step simulation
            stepper.step()
            
            # Progress logging
            if step % 10000 == 0 and step > 0:
                elapsed = time.time() - start_time
                progress = (step / self.max_steps) * 100
                steps_per_sec = step / elapsed if elapsed > 0 else 0
                eta = (self.max_steps - step) / steps_per_sec if steps_per_sec > 0 else 0
                
                logger.info(
                    f"Step {step:,}/{self.max_steps:,} ({progress:.1f}%) | "
                    f"Speed: {steps_per_sec:.1f} steps/s | "
                    f"Elapsed: {elapsed/60:.1f}min | "
                    f"ETA: {eta/60:.1f}min"
                )
            
            # Save snapshots
            if phase2_config.save_snapshots and step % save_interval == 0 and step > 0:
                self._save_snapshot(stepper, step, phase2_config)
            
            # Checkpoints
            if step % checkpoint_interval == 0 and step > 0:
                self._save_checkpoint(stepper, step, phase2_config)
        
        # Final snapshot
        self._save_snapshot(stepper, self.max_steps, phase2_config)
        
        # Extract results
        logger.info("Extracting final results...")
        results = self._extract_results(stepper, sim_config, phase2_config, init_stats)
        
        # Save results
        self._save_results(results)
        
        elapsed_time = time.time() - start_time
        logger.info("=" * 70)
        logger.info("SIMULATION COMPLETE!")
        logger.info(f"Total time: {elapsed_time/60:.1f} minutes ({elapsed_time/3600:.2f} hours)")
        logger.info(f"Average speed: {self.max_steps/elapsed_time:.1f} steps/second")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info("=" * 70)
        
        return results
    
    def _save_snapshot(self, stepper, step: int, phase2_config: Phase2Config):
        """Save simulation snapshot with full data for batch analysis"""
        snapshot_dir = self.output_dir / "snapshots"
        snapshot_dir.mkdir(exist_ok=True)
        
        snapshot_file = snapshot_dir / f"step_{step:08d}.json"
        
        # Extract full state for batch analysis
        n_particles = stepper.particles.particle_count[None]
        
        # Get particle data
        positions, velocities, attributes, active_mask, energies = stepper.particles.get_active_particles()
        
        # Get bonds
        bonds = stepper.binding.get_bonds()
        
        # Get clusters
        clusters = stepper.binding.get_clusters(min_size=2)
        
        snapshot = {
            'step': step,
            'time': stepper.current_time,
            'n_particles': n_particles,
            'scenario': phase2_config.scenario_name,
            'timestamp': datetime.now().isoformat(),
            # Full data for batch analysis
            'positions': positions.tolist() if hasattr(positions, 'tolist') else positions,
            'attributes': attributes.tolist() if hasattr(attributes, 'tolist') else attributes,
            'active_mask': active_mask.tolist() if hasattr(active_mask, 'tolist') else active_mask,
            'bonds': bonds,  # List of (i, j, strength) tuples
            'clusters': clusters  # List of particle index lists
        }
        
        with open(snapshot_file, 'w') as f:
            json.dump(snapshot, f, indent=2, default=self._json_serializer)
        
        logger.info(f"  [SNAPSHOT] Step {step:,} saved with {len(bonds)} bonds, {len(clusters)} clusters")
    
    def _save_checkpoint(self, stepper, step: int, phase2_config: Phase2Config):
        """Save checkpoint for restart capability"""
        checkpoint_dir = self.output_dir / "checkpoints"
        checkpoint_dir.mkdir(exist_ok=True)
        
        checkpoint_file = checkpoint_dir / f"checkpoint_{step:08d}.json"
        
        checkpoint = {
            'step': step,
            'time': stepper.current_time,
            'scenario': phase2_config.scenario_name,
            'timestamp': datetime.now().isoformat(),
            'note': 'Checkpoint for restart (full state would be saved in production)'
        }
        
        with open(checkpoint_file, 'w') as f:
            json.dump(checkpoint, f, indent=2, default=self._json_serializer)
        
        logger.info(f"  [CHECKPOINT] Step {step:,} saved")
    
    def _extract_results(self, stepper, sim_config: SimulationConfig,
                        phase2_config: Phase2Config, init_stats: dict) -> dict:
        """Extract final results"""
        n_particles = stepper.particles.particle_count[None]
        
        results = {
            'scenario': phase2_config.scenario_name,
            'description': phase2_config.description,
            
            'configuration': {
                'max_steps': self.max_steps,
                'dt': sim_config.dt,
                'temperature': phase2_config.temperature,
                'seed': self.seed,
                'initial_molecules': [
                    {'name': mol.name, 'formula': mol.formula, 'count': mol.count}
                    for mol in phase2_config.initial_molecules
                ]
            },
            
            'initialization': init_stats,
            
            'final_state': {
                'step': self.max_steps,
                'time': stepper.current_time,
                'n_particles': n_particles
            },
            
            # Extract molecules from stepper catalog
            'molecules_detected': self._extract_molecules_from_catalog(stepper),
            'reactions_observed': [],  # TODO: Extract from reaction detector
            'novel_molecules': self._extract_novel_molecules_from_catalog(stepper),
            
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'version': '2.0',
                'phase': 2
            }
        }
        
        return results
    
    def _extract_molecules_from_catalog(self, stepper) -> List[Dict]:
        """Extract all detected molecules from stepper catalog"""
        molecules = []
        try:
            # Get all substances from catalog (use most common as proxy for all)
            catalog_substances = stepper.catalog.get_most_common_substances(count=1000)  # Get up to 1000 most common
            
            for substance_record in catalog_substances:
                # Extract basic info from graph and properties
                formula = substance_record.graph.get_canonical_form()  # Use canonical form as formula
                mass = substance_record.properties.get('avg_mass', 0.0)  # Get mass from properties
                complexity = substance_record.graph.get_complexity()
                
                molecule_info = {
                    'id': substance_record.id,
                    'formula': formula,
                    'mass': mass,
                    'complexity': complexity,
                    'first_detected': substance_record.first_seen,
                    'last_seen': substance_record.last_seen,
                    'count': substance_record.occurrence_count,
                    'properties': substance_record.properties
                }
                molecules.append(molecule_info)
                
        except Exception as e:
            logger.warning(f"Failed to extract molecules from catalog: {e}")
            
        logger.info(f"Extracted {len(molecules)} molecules from catalog")
        return molecules
    
    def _extract_novel_molecules_from_catalog(self, stepper) -> List[Dict]:
        """Extract only novel (newly discovered) molecules"""
        novel_molecules = []
        try:
            # Get novel substances (those discovered during this run)
            catalog_substances = stepper.catalog.get_novel_substances(count=1000)  # Get up to 1000 novel
            
            for substance_record in catalog_substances:
                # Extract basic info from graph and properties
                formula = substance_record.graph.get_canonical_form()  # Use canonical form as formula
                mass = substance_record.properties.get('avg_mass', 0.0)  # Get mass from properties
                complexity = substance_record.graph.get_complexity()
                
                molecule_info = {
                    'id': substance_record.id,
                    'formula': formula,
                    'mass': mass,
                    'complexity': complexity,
                    'discovery_time': substance_record.first_seen,
                    'properties': substance_record.properties
                }
                novel_molecules.append(molecule_info)
                    
        except Exception as e:
            logger.warning(f"Failed to extract novel molecules from catalog: {e}")
            
        logger.info(f"Extracted {len(novel_molecules)} novel molecules from catalog")
        return novel_molecules
    
    def _save_results(self, results: dict):
        """Save final results"""
        # Main results file
        results_file = self.output_dir / "results.json"
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=self._json_serializer)
        
        logger.info(f"Results saved to {results_file}")
        
        # Molecules file (for MatcherV2 analysis)
        molecules_file = self.output_dir / "molecules.json"
        with open(molecules_file, 'w') as f:
            json.dump(results.get('molecules_detected', []), f, indent=2, default=self._json_serializer)
        
        # Summary file
        summary_file = self.output_dir / "summary.txt"
        with open(summary_file, 'w') as f:
            f.write(f"Phase 2 Simulation Summary\n")
            f.write(f"=" * 70 + "\n\n")
            f.write(f"Scenario: {results['scenario']}\n")
            f.write(f"Description: {results['description']}\n\n")
            f.write(f"Configuration:\n")
            f.write(f"  Steps: {results['configuration']['max_steps']:,}\n")
            f.write(f"  Temperature: {results['configuration']['temperature']}K\n")
            f.write(f"  Seed: {results['configuration']['seed']}\n\n")
            f.write(f"Initial Molecules:\n")
            for mol in results['configuration']['initial_molecules']:
                f.write(f"  - {mol['name']} ({mol['formula']}): {mol['count']}\n")
            f.write(f"\nFinal State:\n")
            f.write(f"  Particles: {results['final_state']['n_particles']}\n")
            f.write(f"  Time: {results['final_state']['time']:.2f}\n")


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(description="Run Phase 2 Full Simulation")
    parser.add_argument('--config', required=True,
                       help='Path to Phase 2 YAML configuration file')
    parser.add_argument('--output', required=True,
                       help='Output directory')
    parser.add_argument('--steps', type=int, default=10000000,
                       help='Maximum steps (default: 10M)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed (default: 42)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Dry run (load config only, don\'t run)')
    parser.add_argument('--force-cpu', action='store_true',
                       help='Force CPU mode (skip GPU detection, use all CPU cores)')
    
    args = parser.parse_args()
    
    # Create runner
    runner = Phase2FullRunner(
        config_path=args.config,
        output_dir=args.output,
        max_steps=args.steps,
        seed=args.seed,
        force_cpu=args.force_cpu
    )
    
    if args.dry_run:
        logger.info("DRY RUN MODE - Configuration loaded, no simulation will run")
        sim_config, phase2_config = runner.load_configs()
        logger.info(f"Would run {args.steps:,} steps for scenario: {phase2_config.scenario_name}")
        return
    
    # Run simulation
    try:
        results = runner.run()
        logger.info("[SUCCESS] Phase 2 simulation completed!")
        sys.exit(0)
    except Exception as e:
        logger.error(f"[FAILED] Simulation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

