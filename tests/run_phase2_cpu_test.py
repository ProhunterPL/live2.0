#!/usr/bin/env python3
"""
Phase 2 CPU Mode Test Runner
=============================

Forces CPU mode instead of GPU for testing performance.
"""

import sys
import argparse
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from scripts.run_phase2_full import Phase2FullRunner

class CPUPhase2Runner(Phase2FullRunner):
    """Phase 2 runner that forces CPU mode"""
    
    def initialize_taichi(self):
        """Initialize Taichi - Force CPU mode"""
        import logging
        logger = logging.getLogger(__name__)
        
        logger.info("Initializing Taichi...")
        logger.info("FORCING CPU MODE for performance testing")
        
        import taichi as ti
        import multiprocessing
        
        max_threads = multiprocessing.cpu_count()
        ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)
        logger.info(f"Using CPU with {max_threads} threads")

def main():
    parser = argparse.ArgumentParser(description="Run Phase 2 Simulation on CPU")
    parser.add_argument('--config', required=True,
                       help='Path to Phase 2 YAML configuration file')
    parser.add_argument('--output', required=True,
                       help='Output directory')
    parser.add_argument('--steps', type=int, default=10000,
                       help='Maximum steps (default: 10K for test)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed (default: 42)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("PHASE 2 CPU MODE TEST")
    print("=" * 70)
    print(f"Config: {args.config}")
    print(f"Output: {args.output}")
    print(f"Steps: {args.steps:,}")
    print(f"Mode: CPU (forced)")
    print("=" * 70)
    
    # Create runner
    runner = CPUPhase2Runner(
        config_path=args.config,
        output_dir=args.output,
        max_steps=args.steps,
        seed=args.seed
    )
    
    # Run simulation
    try:
        results = runner.run()
        print("\n[SUCCESS] CPU mode test completed!")
        sys.exit(0)
    except Exception as e:
        print(f"\n[FAILED] CPU mode test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

