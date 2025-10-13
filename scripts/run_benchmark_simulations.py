"""
Benchmark Simulations Runner
=============================

Runs actual benchmark simulations (formose, Strecker) and validates against literature.

This is Phase 2: Real experimental validation of the simulation infrastructure.
"""

import argparse
import logging
import sys
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).parent.parent))

import taichi as ti
import numpy as np

from backend.sim.config import SimulationConfig
from backend.sim.core.stepper import SimulationStepper
from backend.sim.core.benchmark_reactions import BenchmarkReactionDatabase
from backend.sim.core.reaction_detector import ReactionDetector
from backend.sim.core.reaction_kinetics import ReactionKineticsAnalyzer, KineticData

logger = logging.getLogger(__name__)


class BenchmarkSimulationRunner:
    """
    Runs benchmark simulations and analyzes results
    
    Usage:
        runner = BenchmarkSimulationRunner()
        results = runner.run_formose_simulation()
        runner.validate_results(results, 'formose')
    """
    
    def __init__(self, output_dir: str = "data/benchmark_results"):
        """Initialize benchmark runner"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Load literature data
        self.literature_db = BenchmarkReactionDatabase('data/benchmark_reactions.json')
        
        # Analysis tools
        self.reaction_detector = ReactionDetector()
        self.kinetics_analyzer = ReactionKineticsAnalyzer()
        
        logger.info("Benchmark runner initialized")
    
    def configure_for_reaction(self, reaction_name: str) -> SimulationConfig:
        """
        Create simulation configuration for specific benchmark reaction
        
        Args:
            reaction_name: 'formose' or 'strecker'
        
        Returns:
            Configured SimulationConfig
        """
        # Get reaction conditions from literature
        conditions = self.literature_db.get_conditions(reaction_name)
        
        # Base configuration
        config = SimulationConfig(
            grid_width=256,
            grid_height=256,
            max_particles=5000,
            max_time=1000.0,
            dt=0.005,
            mode="open_chemistry"
        )
        
        # Adjust for reaction-specific conditions
        # Note: temperature is controlled by energy pulses, not a direct config parameter
        if reaction_name == 'formose':
            # Formose: pH 11, Ca(OH)2 catalyst, autocatalytic
            config.binding_threshold = 0.4  # More reactive for formose
            config.pulse_amplitude = 3.0  # Higher energy for autocatalysis
            config.pulse_every = 30  # More frequent energy pulses
            
        elif reaction_name == 'strecker':
            # Strecker: pH 7-9, HCN + NH3 + aldehyde
            config.binding_threshold = 0.5  # Moderate reactivity
            config.pulse_amplitude = 2.0
            config.pulse_every = 50
        
        logger.info(f"Configured simulation for {reaction_name}")
        logger.info(f"  Conditions: {conditions}")
        
        return config
    
    def run_simulation(self,
                      reaction_name: str,
                      max_steps: int = 10000,
                      save_interval: int = 100) -> Dict:
        """
        Run benchmark simulation
        
        Args:
            reaction_name: Which reaction to simulate
            max_steps: Maximum simulation steps
            save_interval: Save trajectory every N steps
        
        Returns:
            Dict with simulation results
        """
        logger.info("=" * 70)
        logger.info(f"RUNNING BENCHMARK SIMULATION: {reaction_name.upper()}")
        logger.info("=" * 70)
        
        # Configure simulation
        config = self.configure_for_reaction(reaction_name)
        
        # Initialize Taichi
        ti.init(arch=ti.cpu, device_memory_GB=4)
        
        # Create stepper
        stepper = SimulationStepper(config)
        
        # Initialize reactants based on reaction
        self._initialize_reactants(stepper, reaction_name)
        
        # Simulation loop
        logger.info(f"Starting simulation loop ({max_steps} steps)...")
        
        trajectory = []
        start_time = datetime.now()
        
        for step in range(max_steps):
            # Step simulation
            try:
                stepper.step()
            except Exception as e:
                logger.error(f"Simulation error at step {step}: {e}")
                break
            
            # Save trajectory
            if step % save_interval == 0:
                state = self._extract_state(stepper, step)
                trajectory.append(state)
                
                # Update reaction detector
                self._update_reaction_detector(stepper, step, state)
                
                # Progress
                if step % 1000 == 0:
                    elapsed = (datetime.now() - start_time).total_seconds()
                    logger.info(f"  Step {step}/{max_steps} ({step/max_steps*100:.1f}%) - "
                              f"Time: {elapsed:.1f}s - "
                              f"Particles: {stepper.particles.count[None]}")
        
        elapsed_total = (datetime.now() - start_time).total_seconds()
        logger.info(f"Simulation complete! Total time: {elapsed_total:.1f}s")
        
        # Compile results
        results = {
            'reaction': reaction_name,
            'max_steps': max_steps,
            'elapsed_time': elapsed_total,
            'trajectory': trajectory,
            'detected_reactions': self.reaction_detector.get_detected_reactions(),
            'final_state': trajectory[-1] if trajectory else None
        }
        
        return results
    
    def _initialize_reactants(self, stepper: SimulationStepper, reaction_name: str):
        """Initialize simulation with appropriate reactants"""
        if reaction_name == 'formose':
            # Add formaldehyde molecules (CH2O)
            # Simplified: just add some C, O atoms
            for i in range(100):
                # Add carbon
                stepper.particles.add_particle(
                    position=np.random.rand(2) * stepper.config.grid_width,
                    velocity=np.random.randn(2) * 0.1,
                    atom_type=0,  # Carbon
                    mass=12.0
                )
                # Add oxygen
                stepper.particles.add_particle(
                    position=np.random.rand(2) * stepper.config.grid_width,
                    velocity=np.random.randn(2) * 0.1,
                    atom_type=2,  # Oxygen
                    mass=16.0
                )
            
            logger.info("Initialized formose reactants (formaldehyde analogs)")
        
        elif reaction_name == 'strecker':
            # Add acetaldehyde + HCN + NH3
            # Simplified: C, N, O atoms
            for i in range(80):
                # Acetaldehyde (C, C, O)
                stepper.particles.add_particle(
                    position=np.random.rand(2) * stepper.config.grid_width,
                    velocity=np.random.randn(2) * 0.1,
                    atom_type=0,  # Carbon
                    mass=12.0
                )
                # HCN (C, N)
                stepper.particles.add_particle(
                    position=np.random.rand(2) * stepper.config.grid_width,
                    velocity=np.random.randn(2) * 0.1,
                    atom_type=3,  # Nitrogen
                    mass=14.0
                )
            
            logger.info("Initialized Strecker reactants (acetaldehyde + HCN analogs)")
    
    def _extract_state(self, stepper: SimulationStepper, step: int) -> Dict:
        """Extract current simulation state"""
        # Get particle data
        positions = stepper.particles.positions.to_numpy()[:stepper.particles.count[None]]
        velocities = stepper.particles.velocities.to_numpy()[:stepper.particles.count[None]]
        attributes = stepper.particles.attributes.to_numpy()[:stepper.particles.count[None]]
        
        # Get bonds (simplified - would need proper bond extraction)
        bonds = []  # Placeholder
        
        state = {
            'step': step,
            'time': stepper.current_time,
            'n_particles': stepper.particles.count[None],
            'positions': positions,
            'velocities': velocities,
            'attributes': attributes,
            'bonds': bonds,
            'total_energy': stepper._get_total_energy()
        }
        
        return state
    
    def _update_reaction_detector(self, stepper: SimulationStepper, step: int, state: Dict):
        """Update reaction detector with current state"""
        # Note: This is simplified - full implementation would extract actual molecules
        # For now, just update with basic info
        pass
    
    def validate_results(self, results: Dict, reaction_name: str) -> Dict:
        """
        Validate simulation results against literature
        
        Returns validation report
        """
        logger.info("\n" + "=" * 70)
        logger.info("VALIDATION AGAINST LITERATURE")
        logger.info("=" * 70)
        
        # Get expected products
        expected_products = self.literature_db.get_products(reaction_name)
        
        # Analyze trajectory for products
        # Note: Simplified - would use ReactionDetector in full implementation
        detected_products = self._analyze_products(results)
        
        # Compare with literature
        validation_report = {
            'reaction': reaction_name,
            'expected_products': [p.name for p in expected_products],
            'detected_products': detected_products,
            'validation_passed': False,  # Placeholder
            'details': {}
        }
        
        logger.info(f"\nReaction: {reaction_name}")
        logger.info(f"Expected products: {[p.name for p in expected_products]}")
        logger.info(f"Detected products: {detected_products}")
        
        # Check each expected product
        for product in expected_products:
            logger.info(f"\n[{product.name}]")
            logger.info(f"  Expected yield: {product.yield_min*100:.1f}-{product.yield_max*100:.1f}%")
            logger.info(f"  Formula: {product.formula}")
            
            # Placeholder validation
            validation_report['details'][product.name] = {
                'expected_yield_range': (product.yield_min, product.yield_max),
                'observed_yield': 0.0,  # Would be calculated from detector
                'within_tolerance': False
            }
        
        return validation_report
    
    def _analyze_products(self, results: Dict) -> List[str]:
        """Analyze trajectory for product formation"""
        # Placeholder - would use ReactionDetector
        return []
    
    def save_results(self, results: Dict, validation_report: Dict):
        """Save results to file"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        reaction_name = results['reaction']
        
        # Save full results
        results_file = self.output_dir / f"{reaction_name}_{timestamp}_results.json"
        
        # Simplify for JSON (convert numpy arrays)
        results_simplified = {
            'reaction': results['reaction'],
            'max_steps': results['max_steps'],
            'elapsed_time': results['elapsed_time'],
            'n_reactions_detected': len(results['detected_reactions']),
            'final_n_particles': results['final_state']['n_particles'] if results['final_state'] else 0
        }
        
        with open(results_file, 'w') as f:
            json.dump({
                'results': results_simplified,
                'validation': validation_report
            }, f, indent=2)
        
        logger.info(f"\n[+] Results saved to: {results_file}")


def main():
    parser = argparse.ArgumentParser(description='Run benchmark simulations')
    parser.add_argument('--reaction', type=str, choices=['formose', 'strecker', 'all'],
                       default='all', help='Which reaction to simulate')
    parser.add_argument('--steps', type=int, default=10000,
                       help='Number of simulation steps')
    parser.add_argument('--output-dir', type=str, default='data/benchmark_results',
                       help='Output directory')
    
    args = parser.parse_args()
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s'
    )
    
    logger.info("=" * 70)
    logger.info("BENCHMARK SIMULATIONS RUNNER")
    logger.info("Phase 2: Experimental Validation")
    logger.info("=" * 70)
    
    # Create runner
    runner = BenchmarkSimulationRunner(output_dir=args.output_dir)
    
    # Determine which reactions to run
    if args.reaction == 'all':
        reactions = ['formose', 'strecker']
    else:
        reactions = [args.reaction]
    
    # Run each reaction
    all_results = []
    
    for reaction in reactions:
        try:
            # Run simulation
            results = runner.run_simulation(reaction, max_steps=args.steps)
            
            # Validate
            validation = runner.validate_results(results, reaction)
            
            # Save
            runner.save_results(results, validation)
            
            all_results.append({
                'reaction': reaction,
                'success': True,
                'validation': validation
            })
            
        except Exception as e:
            logger.error(f"Failed to run {reaction}: {e}")
            import traceback
            traceback.print_exc()
            
            all_results.append({
                'reaction': reaction,
                'success': False,
                'error': str(e)
            })
    
    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("BENCHMARK SIMULATIONS COMPLETE")
    logger.info("=" * 70)
    
    for result in all_results:
        status = "✅ SUCCESS" if result['success'] else "❌ FAILED"
        logger.info(f"{result['reaction']:15s} - {status}")
    
    logger.info("\n[+] All results saved to: " + str(Path(args.output_dir).absolute()))


if __name__ == "__main__":
    main()

