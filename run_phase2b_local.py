#!/usr/bin/env python3
"""
Phase 2B Local Runner - dla GPU RTX 5070
==========================================

Runs Phase 2B additional simulations locally on GPU.
Based on AWS Phase 2B but optimized for local execution.

Usage:
    python run_phase2b_local.py --scenario miller_urey --runs 10
    python run_phase2b_local.py --all --runs 10
"""

import sys
import argparse
import time
import subprocess
from pathlib import Path
from datetime import datetime
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuration for Phase 2B extended runs
SCENARIOS = {
    'miller_urey': {
        'config': 'aws_test/configs/phase2_miller_urey_extended.yaml',
        'description': 'Miller-Urey Extended (500K steps)',
        'expected_time_per_run': 120  # minutes on GPU
    },
    'hydrothermal': {
        'config': 'aws_test/configs/phase2_hydrothermal_extended.yaml',
        'description': 'Hydrothermal Extended (500K steps)',
        'expected_time_per_run': 120
    },
    'formamide': {
        'config': 'aws_test/configs/phase2_formamide_extended.yaml',
        'description': 'Formamide Extended (500K steps)',
        'expected_time_per_run': 120
    }
}


def run_single_simulation(config_path: str, output_dir: Path, steps: int, seed: int) -> bool:
    """Run a single simulation using run_phase2_full.py"""
    
    cmd = [
        sys.executable,
        "scripts/run_phase2_full.py",
        "--config", config_path,
        "--output", str(output_dir),
        "--steps", str(steps),
        "--seed", str(seed)
    ]
    
    logger.info(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=14400  # 4 hour timeout per run
        )
        
        if result.returncode == 0:
            logger.info(f"✅ Simulation completed successfully")
            return True
        else:
            logger.error(f"❌ Simulation failed:")
            logger.error(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        logger.error(f"⏱️ Simulation timed out after 4 hours")
        return False
    except Exception as e:
        logger.error(f"💥 Simulation crashed: {e}")
        return False


def run_scenario_batch(scenario: str, num_runs: int, seeds: list = None, output_base: str = "results/phase2b_local"):
    """Run multiple simulations for a scenario"""
    
    if scenario not in SCENARIOS:
        logger.error(f"Unknown scenario: {scenario}")
        return []
    
    config = SCENARIOS[scenario]
    output_dir = Path(output_base) / scenario
    
    if seeds is None:
        seeds = list(range(100, 100 + num_runs))  # Start from seed 100
    
    logger.info("=" * 70)
    logger.info(f"PHASE 2B LOCAL BATCH: {scenario.upper()}")
    logger.info("=" * 70)
    logger.info(f"Config: {config['config']}")
    logger.info(f"Runs: {num_runs}")
    logger.info(f"Steps: 500,000 per run")
    logger.info(f"Estimated time: {num_runs * config['expected_time_per_run']:.0f} minutes ({num_runs * config['expected_time_per_run']/60:.1f} hours)")
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 70)
    
    results = []
    
    for i, seed in enumerate(seeds):
        run_id = f"run_{i+1:02d}"
        run_output_dir = output_dir / run_id
        
        logger.info(f"\n[{i+1}/{num_runs}] Running: {run_id} (seed={seed})")
        logger.info("-" * 70)
        
        # Check if already completed
        if run_output_dir.exists() and (run_output_dir / "results.json").exists():
            logger.info(f"⏭️  Skipping - already completed")
            results.append({'run_id': run_id, 'success': True, 'skipped': True})
            continue
        
        start_time = time.time()
        success = run_single_simulation(
            config_path=config['config'],
            output_dir=run_output_dir,
            steps=500000,  # 500K steps as in Phase 2B
            seed=seed
        )
        elapsed_time = time.time() - start_time
        
        if success:
            logger.info(f"✅ Completed in {elapsed_time/60:.1f} minutes")
        else:
            logger.info(f"❌ Failed after {elapsed_time/60:.1f} minutes")
        
        results.append({
            'run_id': run_id,
            'seed': seed,
            'success': success,
            'elapsed_time': elapsed_time
        })
    
    # Summary
    successful = sum(1 for r in results if r.get('success', False))
    logger.info("\n" + "=" * 70)
    logger.info(f"BATCH SUMMARY: {scenario}")
    logger.info(f"Successful: {successful}/{num_runs}")
    logger.info("=" * 70)
    
    return results


def run_all_scenarios(num_runs: int = 10):
    """Run all three scenarios"""
    
    logger.info("🎯 PHASE 2B LOCAL - COMPLETE RUN")
    logger.info("=" * 70)
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Runs per scenario: {num_runs}")
    logger.info(f"Total runs: {num_runs * 3}")
    logger.info("=" * 70)
    
    all_results = {}
    
    for scenario in SCENARIOS.keys():
        logger.info(f"\n\n{'='*70}")
        logger.info(f"STARTING: {scenario.upper()}")
        logger.info(f"{'='*70}\n")
        
        results = run_scenario_batch(scenario, num_runs)
        all_results[scenario] = results
    
    # Final summary
    logger.info("\n\n" + "=" * 70)
    logger.info("PHASE 2B LOCAL - FINAL SUMMARY")
    logger.info("=" * 70)
    
    total_successful = 0
    total_runs = 0
    
    for scenario, results in all_results.items():
        successful = sum(1 for r in results if r.get('success', False))
        total_successful += successful
        total_runs += len(results)
        logger.info(f"{scenario}: {successful}/{len(results)} successful")
    
    logger.info(f"\nTotal: {total_successful}/{total_runs} successful ({100*total_successful/total_runs:.1f}%)")
    logger.info("=" * 70)
    
    return all_results


def main():
    parser = argparse.ArgumentParser(description="Phase 2B Local Runner for RTX 5070")
    parser.add_argument('--scenario', choices=['miller_urey', 'hydrothermal', 'formamide'],
                       help='Scenario to run')
    parser.add_argument('--all', action='store_true',
                       help='Run all scenarios')
    parser.add_argument('--runs', type=int, default=10,
                       help='Number of runs per scenario (default: 10)')
    parser.add_argument('--output', default='results/phase2b_local',
                       help='Output directory (default: results/phase2b_local)')
    
    args = parser.parse_args()
    
    if not args.scenario and not args.all:
        parser.error("Must specify --scenario or --all")
    
    if args.all:
        run_all_scenarios(num_runs=args.runs)
    else:
        run_scenario_batch(args.scenario, args.runs, output_base=args.output)
    
    logger.info("\n✅ Phase 2B Local Runner completed!")


if __name__ == "__main__":
    main()

