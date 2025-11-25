#!/usr/bin/env python3
"""
Batch Script: Build Reaction Networks for All Runs
==================================================

Generates reaction_network.json for all runs in a scenario.
Supports parallel processing for faster execution.

Usage:
    # Sequential (local)
    python scripts/build_reaction_networks_batch.py \
        --scenario hydrothermal_extended
    
    # Parallel (AWS or multi-core)
    python scripts/build_reaction_networks_batch.py \
        --scenario hydrothermal_extended \
        --parallel 4
"""

import sys
import argparse
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from scripts.build_reaction_network_from_snapshots import build_reaction_network

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_single_run(run_dir: Path) -> tuple:
    """Process a single run and return (run_name, success, error)"""
    try:
        network = build_reaction_network(run_dir)
        
        if not network:
            return (run_dir.name, False, "Failed to build network")
        
        # Save network
        output_file = run_dir / "reaction_network.json"
        import json
        with open(output_file, 'w') as f:
            json.dump(network, f, indent=2)
        
        logger.info(f"  [OK] {run_dir.name}: {network['metadata']['n_molecules']} molecules, {network['metadata']['n_reactions']} reactions")
        return (run_dir.name, True, None)
        
    except Exception as e:
        logger.error(f"  [ERROR] {run_dir.name}: {e}")
        return (run_dir.name, False, str(e))


def main():
    parser = argparse.ArgumentParser(
        description="Build reaction networks for all runs in a scenario"
    )
    parser.add_argument(
        '--scenario',
        required=True,
        help='Scenario name (e.g., hydrothermal_extended)'
    )
    parser.add_argument(
        '--base-dir',
        type=str,
        default='results/phase2b_additional',
        help='Base directory (default: results/phase2b_additional)'
    )
    parser.add_argument(
        '--parallel',
        type=int,
        default=1,
        help='Number of parallel workers (default: 1 = sequential)'
    )
    
    args = parser.parse_args()
    
    base_dir = Path(args.base_dir)
    scenario_dir = base_dir / args.scenario
    
    if not scenario_dir.exists():
        logger.error(f"Scenario directory not found: {scenario_dir}")
        sys.exit(1)
    
    # Find all run directories
    run_dirs = sorted([
        d for d in scenario_dir.iterdir() 
        if d.is_dir() and d.name.startswith('run_')
    ], key=lambda x: int(x.name.split('_')[1]))
    
    if not run_dirs:
        logger.error(f"No run directories found in {scenario_dir}")
        sys.exit(1)
    
    logger.info("="*70)
    logger.info("BUILDING REACTION NETWORKS")
    logger.info("="*70)
    logger.info(f"Scenario: {args.scenario}")
    logger.info(f"Runs: {len(run_dirs)}")
    logger.info(f"Parallel workers: {args.parallel}")
    logger.info("="*70)
    
    results = []
    
    if args.parallel > 1:
        # Parallel processing
        logger.info(f"\nProcessing {len(run_dirs)} runs in parallel ({args.parallel} workers)...")
        
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            futures = {executor.submit(process_single_run, run_dir): run_dir 
                      for run_dir in run_dirs}
            
            for future in as_completed(futures):
                run_name, success, error = future.result()
                results.append((run_name, success, error))
    else:
        # Sequential processing
        logger.info(f"\nProcessing {len(run_dirs)} runs sequentially...")
        
        for run_dir in run_dirs:
            run_name, success, error = process_single_run(run_dir)
            results.append((run_name, success, error))
    
    # Summary
    logger.info("\n" + "="*70)
    logger.info("SUMMARY")
    logger.info("="*70)
    
    successful = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]
    
    logger.info(f"Successful: {len(successful)}/{len(results)}")
    logger.info(f"Failed: {len(failed)}/{len(results)}")
    
    if failed:
        logger.warning("\nFailed runs:")
        for run_name, _, error in failed:
            logger.warning(f"  - {run_name}: {error}")
    
    logger.info("="*70)
    
    if len(successful) == len(results):
        logger.info("[OK] All reaction networks generated successfully!")
        sys.exit(0)
    else:
        logger.warning(f"[WARNING] {len(failed)} runs failed")
        sys.exit(1)


if __name__ == "__main__":
    main()

