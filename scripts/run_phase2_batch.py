"""
Phase 2 Batch Simulation Runner
================================

Runs multiple independent simulations for Miller-Urey, Hydrothermal, and Formamide scenarios.

Usage:
    python scripts/run_phase2_batch.py --scenario miller_urey --runs 10
    python scripts/run_phase2_batch.py --scenario hydrothermal --runs 10
    python scripts/run_phase2_batch.py --scenario formamide --runs 10
    python scripts/run_phase2_batch.py --all --runs 10  # All scenarios
"""

import sys
import json
import time
import argparse
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List, Dict

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


SCENARIOS = {
    'miller_urey': {
        'config': 'configs/phase2_miller_urey.yaml',
        'description': 'Miller-Urey reducing atmosphere',
        'expected_time_per_run': 120  # minutes (estimate)
    },
    'hydrothermal': {
        'config': 'configs/phase2_hydrothermal.yaml',
        'description': 'Hydrothermal vent conditions',
        'expected_time_per_run': 120
    },
    'formamide': {
        'config': 'configs/phase2_formamide.yaml',
        'description': 'Formamide-rich with UV radiation',
        'expected_time_per_run': 120
    }
}


class Phase2BatchRunner:
    """Manages batch simulation runs for Phase 2"""
    
    def __init__(self, output_dir: str = "results/phase2"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.log_file = self.output_dir / "batch_run_log.json"
        self.runs_completed = []
        
        # Load existing log if available
        if self.log_file.exists():
            with open(self.log_file, 'r') as f:
                self.runs_completed = json.load(f)
    
    def run_scenario_batch(self,
                          scenario: str,
                          num_runs: int,
                          seeds: List[int] = None) -> List[Dict]:
        """
        Run multiple independent simulations for a scenario
        
        Args:
            scenario: One of 'miller_urey', 'hydrothermal', 'formamide'
            num_runs: Number of independent runs
            seeds: Random seeds (if None, use sequential numbers)
        
        Returns:
            List of run results
        """
        if scenario not in SCENARIOS:
            raise ValueError(f"Unknown scenario: {scenario}")
        
        config = SCENARIOS[scenario]
        
        print("=" * 70)
        print(f"PHASE 2 BATCH RUN: {scenario.upper()}")
        print("=" * 70)
        print(f"Description: {config['description']}")
        print(f"Config: {config['config']}")
        print(f"Number of runs: {num_runs}")
        print(f"Estimated time: {num_runs * config['expected_time_per_run']:.0f} minutes")
        print("=" * 70)
        
        # Generate seeds
        if seeds is None:
            seeds = list(range(42, 42 + num_runs))
        
        results = []
        
        for i, seed in enumerate(seeds):
            run_id = f"{scenario}_run{i+1:02d}_seed{seed}"
            
            print(f"\n[{i+1}/{num_runs}] Running: {run_id}")
            print("-" * 70)
            
            # Check if already completed
            if self._is_run_completed(run_id):
                print(f"  ⏭️  Skipping (already completed)")
                continue
            
            # Run simulation
            start_time = time.time()
            success, output_path, error = self._run_single_simulation(
                scenario=scenario,
                config_path=config['config'],
                run_id=run_id,
                seed=seed
            )
            elapsed_time = time.time() - start_time
            
            # Record result
            result = {
                'run_id': run_id,
                'scenario': scenario,
                'seed': seed,
                'success': success,
                'output_path': str(output_path) if output_path else None,
                'elapsed_time': elapsed_time,
                'timestamp': datetime.now().isoformat(),
                'error': error
            }
            
            results.append(result)
            self.runs_completed.append(result)
            self._save_log()
            
            if success:
                print(f"  ✅ Completed in {elapsed_time/60:.1f} minutes")
                print(f"  📁 Output: {output_path}")
            else:
                print(f"  ❌ Failed: {error}")
        
        # Summary
        print("\n" + "=" * 70)
        print(f"BATCH SUMMARY: {scenario}")
        print("=" * 70)
        successful = sum(1 for r in results if r['success'])
        print(f"Completed: {successful}/{num_runs}")
        print(f"Total time: {sum(r['elapsed_time'] for r in results)/60:.1f} minutes")
        print("=" * 70)
        
        return results
    
    def _run_single_simulation(self,
                              scenario: str,
                              config_path: str,
                              run_id: str,
                              seed: int) -> tuple:
        """
        Run a single simulation
        
        Returns:
            (success, output_path, error_message)
        """
        output_path = self.output_dir / scenario / run_id
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Build command
        # NOTE: This assumes you have a simulation runner script
        # Adjust according to your actual simulation interface
        cmd = [
            sys.executable,
            "backend/sim/run_simulation.py",  # Adjust to your actual script
            "--config", config_path,
            "--seed", str(seed),
            "--output", str(output_path),
            "--max-steps", "10000000",  # 10M steps
            "--save-interval", "50000"
        ]
        
        try:
            # Run simulation
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2 hour timeout
            )
            
            if result.returncode == 0:
                return (True, output_path, None)
            else:
                return (False, None, result.stderr)
        
        except subprocess.TimeoutExpired:
            return (False, None, "Timeout (2 hours)")
        
        except Exception as e:
            return (False, None, str(e))
    
    def _is_run_completed(self, run_id: str) -> bool:
        """Check if run already completed"""
        return any(r['run_id'] == run_id and r['success'] for r in self.runs_completed)
    
    def _save_log(self):
        """Save run log to file"""
        with open(self.log_file, 'w') as f:
            json.dump(self.runs_completed, f, indent=2)
    
    def run_all_scenarios(self, runs_per_scenario: int = 10):
        """Run all three scenarios"""
        all_results = {}
        
        for scenario in SCENARIOS.keys():
            print(f"\n\n{'='*70}")
            print(f"STARTING SCENARIO: {scenario.upper()}")
            print(f"{'='*70}\n")
            
            results = self.run_scenario_batch(scenario, runs_per_scenario)
            all_results[scenario] = results
        
        # Final summary
        print("\n\n" + "=" * 70)
        print("PHASE 2 BATCH RUN - FINAL SUMMARY")
        print("=" * 70)
        
        for scenario, results in all_results.items():
            successful = sum(1 for r in results if r['success'])
            print(f"{scenario}: {successful}/{len(results)} successful")
        
        total_successful = sum(sum(1 for r in results if r['success']) 
                              for results in all_results.values())
        total_runs = sum(len(results) for results in all_results.values())
        
        print(f"\nTotal: {total_successful}/{total_runs} successful")
        print("=" * 70)
        
        return all_results


def main():
    parser = argparse.ArgumentParser(description="Phase 2 Batch Simulation Runner")
    parser.add_argument('--scenario', choices=['miller_urey', 'hydrothermal', 'formamide'],
                       help='Scenario to run (or use --all)')
    parser.add_argument('--all', action='store_true',
                       help='Run all scenarios')
    parser.add_argument('--runs', type=int, default=10,
                       help='Number of independent runs per scenario (default: 10)')
    parser.add_argument('--output', default='results/phase2',
                       help='Output directory (default: results/phase2)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be run without running')
    
    args = parser.parse_args()
    
    if not args.scenario and not args.all:
        parser.error("Must specify --scenario or --all")
    
    # Initialize runner
    runner = Phase2BatchRunner(output_dir=args.output)
    
    if args.dry_run:
        print("DRY RUN MODE - No simulations will be executed")
        print("")
    
    # Run simulations
    if args.all:
        print("Running all scenarios...")
        if not args.dry_run:
            runner.run_all_scenarios(runs_per_scenario=args.runs)
        else:
            for scenario in SCENARIOS.keys():
                print(f"\nWould run {args.runs} simulations of {scenario}")
    else:
        print(f"Running {args.scenario}...")
        if not args.dry_run:
            runner.run_scenario_batch(args.scenario, args.runs)
        else:
            print(f"\nWould run {args.runs} simulations of {args.scenario}")
    
    print("\n✅ Done!")


if __name__ == "__main__":
    main()

