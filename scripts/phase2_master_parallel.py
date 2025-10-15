#!/usr/bin/env python3
"""
Phase 2 Master Orchestrator - PARALLEL VERSION
==============================================

Optimized version that runs simulations in parallel for faster execution.

Usage:
    # Run full pipeline (30 simulations in parallel)
    python scripts/phase2_master_parallel.py --mode full --scenarios all --max-parallel 5
"""

import sys
import argparse
import logging
import json
import subprocess
import threading
import time
from pathlib import Path
from datetime import datetime
from typing import List, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ParallelPhase2MasterOrchestrator:
    """Parallel master orchestrator for Phase 2 pipeline"""
    
    SCENARIOS = {
        'miller_urey': {
            'name': 'Miller-Urey',
            'config': 'configs/phase2_miller_urey.yaml',
            'config_test': 'configs/phase2_miller_urey_test.yaml',
            'runs': 10
        },
        'hydrothermal': {
            'name': 'Hydrothermal Vent',
            'config': 'configs/phase2_hydrothermal.yaml',
            'config_test': 'configs/phase2_hydrothermal_test.yaml',
            'runs': 10
        },
        'formamide': {
            'name': 'Formamide-rich',
            'config': 'configs/phase2_formamide.yaml',
            'config_test': 'configs/phase2_formamide_test.yaml',
            'runs': 10
        }
    }
    
    def __init__(self, mode: str, scenarios: List[str], output_dir: str, max_parallel: int = 5):
        self.mode = mode
        self.scenarios = scenarios
        self.base_output_dir = Path(output_dir)
        self.max_parallel = max_parallel
        
        # Create output directory
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Thread-safe results storage
        self.results_lock = threading.Lock()
        self.results = {
            'scenarios': {},
            'total_runs': 0,
            'successful': 0,
            'failed': 0,
            'running': 0
        }
    
    def run_full_pipeline(self):
        """Run complete Phase 2 pipeline"""
        logger.info("=" * 70)
        logger.info("PHASE 2 MASTER ORCHESTRATOR - PARALLEL VERSION")
        logger.info("=" * 70)
        logger.info(f"Mode: {self.mode}")
        logger.info(f"Output: {self.base_output_dir}")
        logger.info(f"Max parallel: {self.max_parallel}")
        logger.info("=" * 70)
        
        # Phase 1: Run simulations in parallel
        logger.info("\nPhase 1: Running simulations in parallel...")
        simulation_results = self.run_simulations_parallel()
        
        # Phase 2: Analyze results
        logger.info("\nPhase 2: Analyzing results...")
        analysis_results = self.analyze_results(simulation_results)
        
        # Phase 3: Generate final report
        logger.info("\nPhase 3: Generating final report...")
        self.generate_final_report(analysis_results)
        
        logger.info("\n" + "=" * 70)
        logger.info("[SUCCESS] Full pipeline complete!")
        logger.info("=" * 70)
    
    def run_simulations_parallel(self) -> Dict:
        """Run all simulations in parallel"""
        logger.info("Running simulations in parallel...")
        
        # Prepare all simulation tasks
        tasks = []
        for scenario_name in self.scenarios:
            if scenario_name not in self.SCENARIOS:
                logger.warning(f"Unknown scenario: {scenario_name}, skipping")
                continue
            
            scenario = self.SCENARIOS[scenario_name]
            
            # Determine configuration
            if self.mode == 'test':
                config_file = scenario.get('config_test', scenario['config'])
                n_runs = 1
                steps = 10000
            else:
                config_file = scenario['config']
                n_runs = scenario['runs']
                steps = 10000000
            
            logger.info(f"\n{scenario['name']}")
            logger.info("-" * 70)
            logger.info(f"  Config: {config_file}")
            logger.info(f"  Runs: {n_runs}")
            logger.info(f"  Steps: {steps:,}")
            
            # Create tasks for each run
            for run_id in range(1, n_runs + 1):
                output_dir = self.base_output_dir / scenario_name / f"run_{run_id:02d}"
                seed = 42 + run_id
                
                task = {
                    'scenario': scenario_name,
                    'run_id': run_id,
                    'config_file': config_file,
                    'output_dir': str(output_dir),
                    'steps': steps,
                    'seed': seed
                }
                tasks.append(task)
        
        logger.info(f"\nTotal tasks: {len(tasks)}")
        logger.info(f"Max parallel: {self.max_parallel}")
        
        # Run tasks in parallel
        with ThreadPoolExecutor(max_workers=self.max_parallel) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(self._run_single_simulation_task, task): task 
                for task in tasks
            }
            
            # Process completed tasks
            completed = 0
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                completed += 1
                
                try:
                    success = future.result()
                    
                    with self.results_lock:
                        self.results['total_runs'] += 1
                        if success:
                            self.results['successful'] += 1
                        else:
                            self.results['failed'] += 1
                        
                        # Update scenario results
                        if task['scenario'] not in self.results['scenarios']:
                            self.results['scenarios'][task['scenario']] = {
                                'runs': [],
                                'total': 0,
                                'successful': 0
                            }
                        
                        self.results['scenarios'][task['scenario']]['runs'].append({
                            'run_id': task['run_id'],
                            'output_dir': task['output_dir'],
                            'success': success
                        })
                        
                        if success:
                            self.results['scenarios'][task['scenario']]['successful'] += 1
                        self.results['scenarios'][task['scenario']]['total'] += 1
                    
                    status = "SUCCESS" if success else "FAILED"
                    logger.info(f"  [{completed}/{len(tasks)}] {task['scenario']} run {task['run_id']}: {status}")
                    
                except Exception as e:
                    logger.error(f"  [{completed}/{len(tasks)}] {task['scenario']} run {task['run_id']}: ERROR - {e}")
                    
                    with self.results_lock:
                        self.results['total_runs'] += 1
                        self.results['failed'] += 1
        
        # Print summary
        logger.info("\n" + "-" * 70)
        logger.info(f"Simulation summary:")
        logger.info(f"  Total runs: {self.results['total_runs']}")
        logger.info(f"  Successful: {self.results['successful']}")
        logger.info(f"  Failed: {self.results['failed']}")
        
        return self.results
    
    def _run_single_simulation_task(self, task: Dict) -> bool:
        """Run single simulation task"""
        try:
            cmd = [
                'python',
                'scripts/run_phase2_full.py',
                '--config', task['config_file'],
                '--output', task['output_dir'],
                '--steps', str(task['steps']),
                '--seed', str(task['seed'])
            ]
            
            # Run simulation
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=None  # No timeout for long runs
            )
            
            return result.returncode == 0
        
        except subprocess.TimeoutExpired:
            logger.error(f"    [TIMEOUT] {task['scenario']} run {task['run_id']}")
            return False
        except Exception as e:
            logger.error(f"    [ERROR] {task['scenario']} run {task['run_id']}: {e}")
            return False
    
    def analyze_results(self, simulation_results: Dict) -> Dict:
        """Analyze all simulation results"""
        logger.info("Analyzing simulation results...")
        
        analysis_dir = self.base_output_dir.parent / "analysis"
        
        try:
            cmd = [
                'python',
                'scripts/analyze_phase2_batch.py',
                '--input', str(self.base_output_dir),
                '--output', str(analysis_dir),
                '--recursive',
                '--use-matcher'
            ]
            
            logger.info(f"  Command: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )
            
            if result.returncode == 0:
                logger.info(f"  [SUCCESS] Analysis complete")
                logger.info(f"  Output: {analysis_dir}")
                
                # Load analysis results
                analysis_file = analysis_dir / "batch_analysis.json"
                if analysis_file.exists():
                    with open(analysis_file, 'r') as f:
                        return json.load(f)
                
                return {'success': True}
            else:
                logger.error(f"  [FAILED] Analysis failed")
                logger.error(f"  Error: {result.stderr[:500]}")
                return {'success': False, 'error': result.stderr}
        
        except Exception as e:
            logger.error(f"  [ERROR] {e}")
            return {'success': False, 'error': str(e)}
    
    def generate_final_report(self, analysis_results: Dict):
        """Generate final comprehensive report"""
        logger.info("Generating final report...")
        
        report_file = self.base_output_dir.parent / "PHASE2_FINAL_REPORT.txt"
        
        lines = []
        lines.append("=" * 70)
        lines.append("PHASE 2 EXPERIMENTS - FINAL REPORT")
        lines.append("=" * 70)
        lines.append(f"Generated: {datetime.now().isoformat()}")
        lines.append(f"Mode: {self.mode}")
        lines.append(f"Max parallel: {self.max_parallel}")
        lines.append("")
        
        # Simulation summary
        lines.append("SIMULATION SUMMARY:")
        lines.append(f"  Total runs: {self.results['total_runs']}")
        lines.append(f"  Successful: {self.results['successful']}")
        lines.append(f"  Failed: {self.results['failed']}")
        lines.append("")
        
        # Per-scenario summary
        for scenario_name, scenario_data in self.results['scenarios'].items():
            lines.append(f"{scenario_name.upper()}:")
            lines.append(f"  Runs: {scenario_data['successful']}/{scenario_data['total']}")
            lines.append("")
        
        # Analysis summary
        if analysis_results.get('success'):
            lines.append("ANALYSIS: SUCCESS")
        else:
            lines.append("ANALYSIS: FAILED")
            if 'error' in analysis_results:
                lines.append(f"Error: {analysis_results['error']}")
        
        lines.append("")
        lines.append("=" * 70)
        
        # Write report
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Final report saved: {report_file}")


def main():
    parser = argparse.ArgumentParser(description='Phase 2 Master Orchestrator - Parallel Version')
    parser.add_argument('--mode', choices=['test', 'full'], required=True, help='Execution mode')
    parser.add_argument('--scenarios', nargs='+', default=['all'], help='Scenarios to run')
    parser.add_argument('--output', default='results/phase2_parallel', help='Output directory')
    parser.add_argument('--max-parallel', type=int, default=5, help='Maximum parallel simulations')
    
    args = parser.parse_args()
    
    # Handle 'all' scenarios
    if 'all' in args.scenarios:
        scenarios = ['miller_urey', 'hydrothermal', 'formamide']
    else:
        scenarios = args.scenarios
    
    # Create orchestrator
    orchestrator = ParallelPhase2MasterOrchestrator(
        mode=args.mode,
        scenarios=scenarios,
        output_dir=args.output,
        max_parallel=args.max_parallel
    )
    
    # Run pipeline
    orchestrator.run_full_pipeline()


if __name__ == "__main__":
    main()
