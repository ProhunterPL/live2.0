"""
Phase 2 Master Orchestrator
============================

Complete pipeline management for Phase 2 experiments.
Coordinates simulation runs, analysis, and reporting.

Usage:
    # Run full pipeline (30 simulations)
    python scripts/phase2_master.py --mode full --scenarios all
    
    # Test mode (quick)
    python scripts/phase2_master.py --mode test --scenarios miller_urey
    
    # Only analyze existing results
    python scripts/phase2_master.py --mode analyze --input results/phase2
"""

import sys
import argparse
import logging
import json
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List, Dict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Phase2MasterOrchestrator:
    """Master orchestrator for Phase 2 pipeline"""
    
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
            'config_test': 'configs/phase2_hydrothermal.yaml',  # Would create test version
            'runs': 10
        },
        'formamide': {
            'name': 'Formamide-rich',
            'config': 'configs/phase2_formamide.yaml',
            'config_test': 'configs/phase2_formamide.yaml',  # Would create test version
            'runs': 10
        }
    }
    
    def __init__(self, mode: str, base_output_dir: str = "results/phase2"):
        """
        Initialize orchestrator
        
        Args:
            mode: 'test', 'full', or 'analyze'
            base_output_dir: Base directory for all results
        """
        self.mode = mode
        self.base_output_dir = Path(base_output_dir)
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("=" * 70)
        logger.info("PHASE 2 MASTER ORCHESTRATOR")
        logger.info("=" * 70)
        logger.info(f"Mode: {mode}")
        logger.info(f"Output: {self.base_output_dir}")
        logger.info("=" * 70)
    
    def run_full_pipeline(self, scenarios: List[str], test_mode: bool = False):
        """
        Run complete Phase 2 pipeline
        
        Args:
            scenarios: List of scenario names to run
            test_mode: Use test configurations (shorter runs)
        """
        logger.info("\n" + "=" * 70)
        logger.info("STARTING FULL PIPELINE")
        logger.info("=" * 70)
        logger.info(f"Scenarios: {', '.join(scenarios)}")
        logger.info(f"Test mode: {test_mode}")
        
        # Phase 1: Run simulations
        logger.info("\nPhase 1: Running simulations...")
        simulation_results = self.run_simulations(scenarios, test_mode)
        
        # Phase 2: Analyze results
        logger.info("\nPhase 2: Analyzing results...")
        analysis_results = self.analyze_results(simulation_results)
        
        # Phase 3: Generate final report
        logger.info("\nPhase 3: Generating final report...")
        self.generate_final_report(analysis_results)
        
        logger.info("\n" + "=" * 70)
        logger.info("[SUCCESS] Full pipeline complete!")
        logger.info("=" * 70)
    
    def run_simulations(self, scenarios: List[str], test_mode: bool = False) -> Dict:
        """
        Run all simulations for specified scenarios
        
        Args:
            scenarios: List of scenario names
            test_mode: Use test configurations
        
        Returns:
            Dictionary of simulation results
        """
        logger.info("Running simulations...")
        
        results = {
            'scenarios': {},
            'total_runs': 0,
            'successful': 0,
            'failed': 0
        }
        
        for scenario_name in scenarios:
            if scenario_name not in self.SCENARIOS:
                logger.warning(f"Unknown scenario: {scenario_name}, skipping")
                continue
            
            scenario = self.SCENARIOS[scenario_name]
            logger.info(f"\n{scenario['name']}")
            logger.info("-" * 70)
            
            # Determine configuration
            if test_mode:
                config_file = scenario.get('config_test', scenario['config'])
                n_runs = 1  # Only 1 run in test mode
                steps = 10000  # Short run
            else:
                config_file = scenario['config']
                n_runs = scenario['runs']
                steps = 10000000  # Full 10M steps
            
            logger.info(f"  Config: {config_file}")
            logger.info(f"  Runs: {n_runs}")
            logger.info(f"  Steps: {steps:,}")
            
            # Run simulations
            scenario_results = []
            for run_id in range(1, n_runs + 1):
                output_dir = self.base_output_dir / scenario_name / f"run_{run_id:02d}"
                seed = 42 + run_id
                
                logger.info(f"\n  [{run_id}/{n_runs}] Running {scenario_name} run {run_id}")
                
                success = self._run_single_simulation(
                    config_file=config_file,
                    output_dir=str(output_dir),
                    steps=steps,
                    seed=seed
                )
                
                scenario_results.append({
                    'run_id': run_id,
                    'output_dir': str(output_dir),
                    'success': success
                })
                
                results['total_runs'] += 1
                if success:
                    results['successful'] += 1
                else:
                    results['failed'] += 1
            
            results['scenarios'][scenario_name] = {
                'runs': scenario_results,
                'total': n_runs,
                'successful': sum(1 for r in scenario_results if r['success'])
            }
        
        logger.info("\n" + "-" * 70)
        logger.info(f"Simulation summary:")
        logger.info(f"  Total runs: {results['total_runs']}")
        logger.info(f"  Successful: {results['successful']}")
        logger.info(f"  Failed: {results['failed']}")
        
        return results
    
    def _run_single_simulation(self, config_file: str, output_dir: str,
                               steps: int, seed: int) -> bool:
        """
        Run single simulation
        
        Returns:
            True if successful
        """
        try:
            cmd = [
                'python',
                'scripts/run_phase2_full.py',
                '--config', config_file,
                '--output', output_dir,
                '--steps', str(steps),
                '--seed', str(seed)
            ]
            
            logger.info(f"    Command: {' '.join(cmd)}")
            
            # Run simulation
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=None  # No timeout for long runs
            )
            
            if result.returncode == 0:
                logger.info(f"    [SUCCESS]")
                return True
            else:
                logger.error(f"    [FAILED] Exit code: {result.returncode}")
                logger.error(f"    Error: {result.stderr[:500]}")
                return False
        
        except subprocess.TimeoutExpired:
            logger.error(f"    [TIMEOUT] Simulation took too long")
            return False
        except Exception as e:
            logger.error(f"    [ERROR] {e}")
            return False
    
    def analyze_results(self, simulation_results: Dict) -> Dict:
        """
        Analyze all simulation results
        
        Args:
            simulation_results: Results from run_simulations
        
        Returns:
            Analysis results
        """
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
        lines.append("")
        
        # Add summary from analysis
        if analysis_results.get('success'):
            summary = analysis_results.get('summary', {})
            lines.append("SUMMARY")
            lines.append("-" * 70)
            lines.append(f"Total runs: {summary.get('total_runs', 0)}")
            lines.append(f"Successful: {summary.get('successful_runs', 0)}")
            lines.append(f"Total molecules: {summary.get('overall', {}).get('total_unique_molecules', 0)}")
            lines.append(f"Unique formulas: {summary.get('overall', {}).get('total_unique_formulas', 0)}")
        else:
            lines.append("ANALYSIS INCOMPLETE OR FAILED")
        
        lines.append("")
        lines.append("=" * 70)
        lines.append("See analysis/ directory for detailed results")
        lines.append("=" * 70)
        
        report_text = "\n".join(lines)
        
        with open(report_file, 'w') as f:
            f.write(report_text)
        
        logger.info(f"Final report: {report_file}")
        logger.info("\n" + report_text)


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(description="Phase 2 Master Orchestrator")
    parser.add_argument('--mode', required=True,
                       choices=['test', 'full', 'analyze'],
                       help='Pipeline mode')
    parser.add_argument('--scenarios', nargs='+',
                       choices=list(Phase2MasterOrchestrator.SCENARIOS.keys()) + ['all'],
                       default=['all'],
                       help='Scenarios to run')
    parser.add_argument('--output', default='results/phase2',
                       help='Base output directory')
    parser.add_argument('--input', 
                       help='Input directory for analyze mode')
    
    args = parser.parse_args()
    
    # Resolve scenarios
    if 'all' in args.scenarios:
        scenarios = list(Phase2MasterOrchestrator.SCENARIOS.keys())
    else:
        scenarios = args.scenarios
    
    # Create orchestrator
    orchestrator = Phase2MasterOrchestrator(
        mode=args.mode,
        base_output_dir=args.output
    )
    
    try:
        if args.mode in ['test', 'full']:
            test_mode = (args.mode == 'test')
            orchestrator.run_full_pipeline(scenarios, test_mode=test_mode)
        
        elif args.mode == 'analyze':
            if not args.input:
                logger.error("--input required for analyze mode")
                sys.exit(1)
            
            # Just run analysis
            results = orchestrator.analyze_results({})
            orchestrator.generate_final_report(results)
        
        logger.info("\n[SUCCESS] Pipeline complete!")
        sys.exit(0)
    
    except KeyboardInterrupt:
        logger.warning("\n[INTERRUPTED] Pipeline stopped by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"\n[FAILED] Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

