#!/usr/bin/env python3
"""
Phase 2 Master Orchestrator - 1M Steps Version
==============================================

Optimized version for 1M steps instead of 10M
10x faster execution, suitable for parallel batch processing

Usage:
    # Test mode (quick validation)
    python scripts/phase2_master_1M.py --mode test --scenarios miller_urey
    
    # Full mode (50 runs recommended for 1M steps)
    python scripts/phase2_master_1M.py --mode full --scenarios all --max-parallel 4
"""

import sys
import argparse
import logging
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List, Dict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Phase2Master1M:
    """Master orchestrator for 1M step simulations"""
    
    SCENARIOS = {
        'miller_urey': {
            'name': 'Miller-Urey',
            'config': 'configs/phase2_miller_urey_1M.yaml',
            'runs': 50  # More runs since each is 10x faster
        },
        'hydrothermal': {
            'name': 'Hydrothermal Vent',
            'config': 'configs/phase2_hydrothermal_1M.yaml',
            'runs': 50
        },
        'formamide': {
            'name': 'Formamide-rich',
            'config': 'configs/phase2_formamide_1M.yaml',
            'runs': 50
        }
    }
    
    def __init__(self, mode: str, scenarios: List[str], output_dir: str, max_parallel: int = 4):
        self.mode = mode
        self.scenarios = scenarios
        self.base_output_dir = Path(output_dir)
        self.max_parallel = max_parallel
        
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        
        self.results = {
            'scenarios': {},
            'total_runs': 0,
            'successful': 0,
            'failed': 0
        }
    
    def run_pipeline(self):
        """Run complete pipeline"""
        logger.info("=" * 70)
        logger.info("PHASE 2 MASTER ORCHESTRATOR - 1M STEPS VERSION")
        logger.info("=" * 70)
        logger.info(f"Mode: {self.mode}")
        logger.info(f"Output: {self.base_output_dir}")
        logger.info(f"Max parallel: {self.max_parallel}")
        logger.info("=" * 70)
        
        # Run simulations
        logger.info("\nRunning simulations...")
        self.run_simulations()
        
        # Summary
        logger.info("\n" + "=" * 70)
        logger.info("PIPELINE COMPLETE!")
        logger.info(f"Total runs: {self.results['total_runs']}")
        logger.info(f"Successful: {self.results['successful']}")
        logger.info(f"Failed: {self.results['failed']}")
        logger.info("=" * 70)
        
        # Save report
        self.save_report()
    
    def run_simulations(self):
        """Run all simulations"""
        for scenario_name in self.scenarios:
            if scenario_name not in self.SCENARIOS:
                logger.warning(f"Unknown scenario: {scenario_name}")
                continue
            
            scenario = self.SCENARIOS[scenario_name]
            n_runs = 1 if self.mode == 'test' else scenario['runs']
            
            logger.info(f"\n{scenario['name']}")
            logger.info("-" * 70)
            logger.info(f"  Config: {scenario['config']}")
            logger.info(f"  Runs: {n_runs}")
            
            for run_id in range(1, n_runs + 1):
                output_dir = self.base_output_dir / scenario_name / f"run_{run_id:03d}"
                seed = 42 + run_id
                
                logger.info(f"\n  [{run_id}/{n_runs}] Starting run {run_id}...")
                
                success = self.run_single_simulation(
                    scenario['config'],
                    output_dir,
                    seed
                )
                
                self.results['total_runs'] += 1
                if success:
                    self.results['successful'] += 1
                    logger.info(f"  [{run_id}/{n_runs}] ✅ SUCCESS")
                else:
                    self.results['failed'] += 1
                    logger.error(f"  [{run_id}/{n_runs}] ❌ FAILED")
    
    def run_single_simulation(self, config: str, output_dir: Path, seed: int) -> bool:
        """Run single simulation"""
        try:
            cmd = [
                'python',
                'scripts/run_phase2_full.py',
                '--config', config,
                '--output', str(output_dir),
                '--steps', '1000000',  # 1M steps
                '--seed', str(seed)
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=None
            )
            
            return result.returncode == 0
            
        except Exception as e:
            logger.error(f"    ERROR: {e}")
            return False
    
    def save_report(self):
        """Save final report"""
        report_file = self.base_output_dir / "PHASE2_1M_REPORT.txt"
        
        lines = []
        lines.append("=" * 70)
        lines.append("PHASE 2 - 1M STEPS EXPERIMENTS - FINAL REPORT")
        lines.append("=" * 70)
        lines.append(f"Generated: {datetime.now().isoformat()}")
        lines.append(f"Mode: {self.mode}")
        lines.append("")
        lines.append("SIMULATION SUMMARY:")
        lines.append(f"  Total runs: {self.results['total_runs']}")
        lines.append(f"  Successful: {self.results['successful']}")
        lines.append(f"  Failed: {self.results['failed']}")
        lines.append(f"  Success rate: {100*self.results['successful']/max(1,self.results['total_runs']):.1f}%")
        lines.append("")
        lines.append("NOTES:")
        lines.append("  - Each simulation: 1M steps (optimized)")
        lines.append("  - Expected duration: ~4 days per run")
        lines.append("  - Recommended: 50 runs per scenario for statistical power")
        lines.append("")
        lines.append("=" * 70)
        
        with open(report_file, 'w') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"\nReport saved: {report_file}")


def main():
    parser = argparse.ArgumentParser(description='Phase 2 Master - 1M Steps Version')
    parser.add_argument('--mode', choices=['test', 'full'], required=True)
    parser.add_argument('--scenarios', nargs='+', default=['all'])
    parser.add_argument('--output', default='results/phase2_1M')
    parser.add_argument('--max-parallel', type=int, default=4)
    
    args = parser.parse_args()
    
    # Handle 'all' scenarios
    if 'all' in args.scenarios:
        scenarios = ['miller_urey', 'hydrothermal', 'formamide']
    else:
        scenarios = args.scenarios
    
    # Create orchestrator
    orchestrator = Phase2Master1M(
        mode=args.mode,
        scenarios=scenarios,
        output_dir=args.output,
        max_parallel=args.max_parallel
    )
    
    # Run pipeline
    orchestrator.run_pipeline()


if __name__ == "__main__":
    main()

