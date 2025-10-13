"""
Standalone Simulation Runner for Live 2.0
==========================================

Runs simulations from YAML configuration files.
Designed for Phase 2 batch processing.

Usage:
    python backend/sim/run_simulation.py --config configs/phase2_miller_urey.yaml \\
                                         --output results/miller_urey/run01 \\
                                         --seed 42 \\
                                         --max-steps 10000000
"""

import sys
import argparse
import yaml
import json
import time
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import taichi as ti
import numpy as np

from backend.sim.config import SimulationConfig
from backend.sim.core.stepper import SimulationStepper
from backend.sim.core.particles import ParticleSystem
from backend.sim.io.snapshot import SnapshotManager


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SimulationRunner:
    """Standalone simulation runner"""
    
    def __init__(self,
                 config_path: str,
                 output_dir: str,
                 seed: int = 42,
                 max_steps: int = None):
        """
        Initialize simulation runner
        
        Args:
            config_path: Path to YAML configuration
            output_dir: Output directory for results
            seed: Random seed
            max_steps: Maximum steps (overrides config)
        """
        self.config_path = Path(config_path)
        self.output_dir = Path(output_dir)
        self.seed = seed
        self.max_steps_override = max_steps
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging to file
        log_file = self.output_dir / "simulation.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(file_handler)
        
        logger.info(f"Simulation Runner initialized")
        logger.info(f"Config: {self.config_path}")
        logger.info(f"Output: {self.output_dir}")
        logger.info(f"Seed: {self.seed}")
    
    def load_config(self) -> SimulationConfig:
        """Load configuration from YAML"""
        logger.info(f"Loading configuration from {self.config_path}")
        
        with open(self.config_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        
        # Extract simulation parameters
        sim_config = config_dict.get('simulation', {})
        
        # Create SimulationConfig
        # NOTE: This is simplified - you may need to map YAML structure to SimulationConfig
        config = SimulationConfig(
            n_particles=sim_config.get('n_particles', 2000),
            max_steps=self.max_steps_override or sim_config.get('max_steps', 1000000),
            dt=sim_config.get('dt', 0.001),
            box_size=sim_config.get('box_size', 100.0),
            temperature=sim_config.get('target_temperature', 298.0),
            mode='preset_prebiotic',  # Must be 'preset_prebiotic' or 'open_chemistry'
            seed=self.seed
        )
        
        logger.info(f"Configuration loaded:")
        logger.info(f"  Particles: {config.n_particles}")
        logger.info(f"  Max steps: {config.max_steps}")
        logger.info(f"  Temperature: {config.temperature}K")
        logger.info(f"  Box size: {config.box_size}")
        
        return config
    
    def initialize_taichi(self):
        """Initialize Taichi"""
        logger.info("Initializing Taichi...")
        
        # Try GPU first, fallback to CPU
        try:
            ti.init(arch=ti.cuda, device_memory_GB=4.0)
            logger.info("  Using CUDA (GPU)")
        except:
            try:
                ti.init(arch=ti.gpu)
                logger.info("  Using generic GPU")
            except:
                ti.init(arch=ti.cpu)
                logger.info("  Using CPU (fallback)")
    
    def run(self) -> Dict[str, Any]:
        """
        Run simulation
        
        Returns:
            Results dictionary
        """
        start_time = time.time()
        
        # Load configuration
        config = self.load_config()
        
        # Initialize Taichi
        self.initialize_taichi()
        
        # Create stepper
        logger.info("Creating simulation stepper...")
        stepper = SimulationStepper(config)
        
        # Initialize particles
        logger.info("Initializing particles...")
        # NOTE: For Phase 2, we need to load initial molecules from config
        # For now, using random initialization
        self._initialize_particles(stepper, config)
        
        # Run simulation
        logger.info(f"Starting simulation ({config.max_steps} steps)...")
        logger.info("=" * 70)
        
        save_interval = config.max_steps // 20  # Save 20 snapshots
        save_interval = max(save_interval, 50000)  # At least every 50k
        
        for step in range(config.max_steps):
            # Step simulation
            stepper.step()
            
            # Progress logging
            if step % 10000 == 0:
                progress = (step / config.max_steps) * 100
                elapsed = time.time() - start_time
                eta = (elapsed / (step + 1)) * (config.max_steps - step)
                
                logger.info(
                    f"Step {step:,}/{config.max_steps:,} ({progress:.1f}%) | "
                    f"Elapsed: {elapsed/60:.1f}min | ETA: {eta/60:.1f}min"
                )
            
            # Save snapshots
            if step % save_interval == 0:
                self._save_snapshot(stepper, step)
        
        # Final snapshot
        self._save_snapshot(stepper, config.max_steps)
        
        # Extract results
        results = self._extract_results(stepper, config)
        
        # Save results
        self._save_results(results)
        
        elapsed_time = time.time() - start_time
        logger.info("=" * 70)
        logger.info(f"Simulation complete!")
        logger.info(f"Total time: {elapsed_time/60:.1f} minutes")
        logger.info(f"Output: {self.output_dir}")
        
        return results
    
    def _initialize_particles(self, stepper: SimulationStepper, config: SimulationConfig):
        """Initialize particles (placeholder - needs proper molecule loading)"""
        # For Phase 2, we need to load initial molecules from YAML config
        # For now, just use random positions
        
        particles = stepper.particle_system
        
        for i in range(config.n_particles):
            # Random position in box
            x = np.random.uniform(0, config.box_size)
            y = np.random.uniform(0, config.box_size)
            
            # Random velocity (Maxwell-Boltzmann)
            vx = np.random.normal(0, 1.0)
            vy = np.random.normal(0, 1.0)
            
            # Random atom type (simplified)
            atom_type = np.random.choice([0, 1, 2, 3])  # H, C, N, O
            
            particles.positions[i] = [x, y]
            particles.velocities[i] = [vx, vy]
            particles.types[i] = atom_type
            particles.masses[i] = 1.0
            particles.active[i] = True
        
        logger.info(f"  Initialized {config.n_particles} particles")
    
    def _save_snapshot(self, stepper: SimulationStepper, step: int):
        """Save simulation snapshot"""
        snapshot_dir = self.output_dir / "snapshots"
        snapshot_dir.mkdir(exist_ok=True)
        
        snapshot_file = snapshot_dir / f"step_{step:08d}.json"
        
        # Extract state
        snapshot = {
            'step': step,
            'time': stepper.current_time,
            'positions': stepper.particle_system.positions.to_numpy().tolist(),
            'velocities': stepper.particle_system.velocities.to_numpy().tolist(),
            'types': stepper.particle_system.types.to_numpy().tolist(),
            'active': stepper.particle_system.active.to_numpy().tolist()
        }
        
        with open(snapshot_file, 'w') as f:
            json.dump(snapshot, f)
        
        logger.debug(f"  Snapshot saved: {snapshot_file.name}")
    
    def _extract_results(self, stepper: SimulationStepper, config: SimulationConfig) -> Dict:
        """Extract results from completed simulation"""
        logger.info("Extracting results...")
        
        # Get final state
        positions = stepper.particle_system.positions.to_numpy()
        velocities = stepper.particle_system.velocities.to_numpy()
        types = stepper.particle_system.types.to_numpy()
        active = stepper.particle_system.active.to_numpy()
        
        # Count active particles
        n_active = np.sum(active)
        
        # Get metrics
        metrics_history = {}
        if hasattr(stepper, 'metrics_collector'):
            # Extract metrics if available
            pass
        
        results = {
            'config': {
                'n_particles': config.n_particles,
                'max_steps': config.max_steps,
                'dt': config.dt,
                'temperature': config.temperature,
                'seed': self.seed
            },
            'final_state': {
                'n_active': int(n_active),
                'positions': positions[active].tolist(),
                'velocities': velocities[active].tolist(),
                'types': types[active].tolist()
            },
            'molecules': [],  # Placeholder - need cluster detection
            'reactions': [],  # Placeholder - need reaction tracking
            'metrics': metrics_history
        }
        
        logger.info(f"  Final active particles: {n_active}")
        
        return results
    
    def _save_results(self, results: Dict):
        """Save final results"""
        # Save main results
        results_file = self.output_dir / "results.json"
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Results saved to {results_file}")
        
        # Save molecules (placeholder)
        molecules_file = self.output_dir / "molecules.json"
        with open(molecules_file, 'w') as f:
            json.dump(results['molecules'], f, indent=2)
        
        # Save reactions (placeholder)
        reactions_file = self.output_dir / "reactions.json"
        with open(reactions_file, 'w') as f:
            json.dump(results['reactions'], f, indent=2)


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(description="Run Live 2.0 simulation")
    parser.add_argument('--config', required=True,
                       help='Path to YAML configuration file')
    parser.add_argument('--output', required=True,
                       help='Output directory')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed (default: 42)')
    parser.add_argument('--max-steps', type=int, default=None,
                       help='Maximum steps (overrides config)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Dry run (show config, don\'t run)')
    
    args = parser.parse_args()
    
    # Create runner
    runner = SimulationRunner(
        config_path=args.config,
        output_dir=args.output,
        seed=args.seed,
        max_steps=args.max_steps
    )
    
    if args.dry_run:
        logger.info("DRY RUN MODE - No simulation will run")
        config = runner.load_config()
        logger.info(f"Would run simulation with {config.max_steps:,} steps")
        return
    
    # Run simulation
    try:
        results = runner.run()
        logger.info("✅ Simulation completed successfully!")
        sys.exit(0)
    except Exception as e:
        logger.error(f"❌ Simulation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

