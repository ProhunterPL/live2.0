#!/usr/bin/env python3
"""
Phase 2B Additional Runs - Master Runner
========================================

Runs 30 additional simulations for Phase 2B:
- 10 Miller-Urey extended (500K steps)
- 10 Hydrothermal extended (500K steps)  
- 10 Formamide extended (500K steps)

Includes:
- Progress tracking
- Error handling
- Performance monitoring
- Automatic analysis
"""

import os
import sys
import json
import time
import subprocess
import argparse
from pathlib import Path
from datetime import datetime
import logging

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

class Phase2BRunner:
    def __init__(self, base_output_dir="results/phase2b_additional"):
        self.base_output_dir = Path(base_output_dir)
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self.setup_logging()
        
        # Simulation configurations (SUPER FAST MODE - optimized)
        self.configs = {
            "miller_urey_extended": {
                "config_file": "configs/phase2_miller_urey_extended_SUPER_FAST.yaml",
                "runs": 10,
                "seeds": list(range(100, 110)),
                "description": "Extended Miller-Urey (500K steps) - SUPER FAST MODE"
            },
            "hydrothermal_extended": {
                "config_file": "configs/phase2_hydrothermal_extended_SUPER_FAST.yaml", 
                "runs": 10,
                "seeds": list(range(110, 120)),
                "description": "Extended Hydrothermal (500K steps) - SUPER FAST MODE"
            },
            "formamide_extended": {
                "config_file": "configs/phase2_formamide_extended_SUPER_FAST.yaml",
                "runs": 10, 
                "seeds": list(range(120, 130)),
                "description": "Extended Formamide (500K steps) - SUPER FAST MODE"
            }
        }
        
        # Results tracking
        self.results = {
            "start_time": datetime.now().isoformat(),
            "total_runs": 30,
            "completed_runs": 0,
            "failed_runs": 0,
            "scenarios": {}
        }
        
    def setup_logging(self):
        """Setup logging for the runner"""
        log_dir = self.base_output_dir / "logs"
        log_dir.mkdir(exist_ok=True)
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_dir / "phase2b_runner.log"),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def run_single_simulation(self, scenario_name, run_id, seed, config_file):
        """Run a single simulation"""
        self.logger.info(f"Starting {scenario_name} run {run_id} (seed {seed})")
        
        # Create output directory
        output_dir = self.base_output_dir / scenario_name / f"run_{run_id}"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare command - use python3 on Linux/AWS
        import sys
        python_cmd = sys.executable  # Use current Python interpreter
        
        cmd = [
            python_cmd, "scripts/run_phase2_full.py",
            "--config", config_file,
            "--output", str(output_dir),
            "--seed", str(seed),
            "--steps", "500000"  # 500K steps
        ]
        
        # Run simulation
        start_time = time.time()
        try:
            result = subprocess.run(
                cmd,
                cwd=project_root,
                capture_output=True,
                text=True,
                timeout=21600  # 6 hours timeout (for 500K steps in SUPER FAST MODE)
            )
            
            duration = time.time() - start_time
            
            if result.returncode == 0:
                self.logger.info(f"✅ {scenario_name} run {run_id} completed in {duration:.1f}s")
                return {
                    "status": "success",
                    "duration": duration,
                    "output_dir": str(output_dir)
                }
            else:
                self.logger.error(f"❌ {scenario_name} run {run_id} failed: {result.stderr}")
                return {
                    "status": "failed",
                    "duration": duration,
                    "error": result.stderr
                }
                
        except subprocess.TimeoutExpired:
            self.logger.error(f"⏰ {scenario_name} run {run_id} timed out after 6 hours")
            return {
                "status": "timeout",
                "duration": 21600,
                "error": "Simulation timed out after 6 hours"
            }
        except Exception as e:
            self.logger.error(f"💥 {scenario_name} run {run_id} crashed: {str(e)}")
            return {
                "status": "crashed",
                "duration": time.time() - start_time,
                "error": str(e)
            }
    
    def run_scenario(self, scenario_name, config_info):
        """Run all simulations for a scenario"""
        self.logger.info(f"🚀 Starting {scenario_name}: {config_info['description']}")
        
        scenario_results = {
            "scenario": scenario_name,
            "description": config_info["description"],
            "total_runs": config_info["runs"],
            "completed_runs": 0,
            "failed_runs": 0,
            "runs": []
        }
        
        for i in range(config_info["runs"]):
            run_id = i + 1
            seed = config_info["seeds"][i]
            
            result = self.run_single_simulation(
                scenario_name, run_id, seed, config_info["config_file"]
            )
            
            scenario_results["runs"].append({
                "run_id": run_id,
                "seed": seed,
                **result
            })
            
            if result["status"] == "success":
                scenario_results["completed_runs"] += 1
                self.results["completed_runs"] += 1
            else:
                scenario_results["failed_runs"] += 1
                self.results["failed_runs"] += 1
            
            # Save intermediate results
            self.save_results()
            
            # Progress update
            total_progress = self.results["completed_runs"] + self.results["failed_runs"]
            self.logger.info(f"📊 Progress: {total_progress}/30 runs completed")
        
        self.results["scenarios"][scenario_name] = scenario_results
        self.logger.info(f"✅ {scenario_name} completed: {scenario_results['completed_runs']}/{scenario_results['total_runs']} successful")
        
    def run_all_simulations(self):
        """Run all Phase 2B simulations"""
        self.logger.info("🎯 Starting Phase 2B Additional Runs")
        self.logger.info(f"📁 Output directory: {self.base_output_dir}")
        
        start_time = time.time()
        
        # Run each scenario
        for scenario_name, config_info in self.configs.items():
            self.run_scenario(scenario_name, config_info)
        
        total_duration = time.time() - start_time
        
        # Final results
        self.results["end_time"] = datetime.now().isoformat()
        self.results["total_duration"] = total_duration
        
        self.logger.info(f"🏁 Phase 2B completed in {total_duration/3600:.1f} hours")
        self.logger.info(f"📊 Final results: {self.results['completed_runs']}/{self.results['total_runs']} successful")
        
        # Save final results
        self.save_results()
        
        # Generate summary report
        self.generate_summary_report()
        
    def save_results(self):
        """Save results to JSON file"""
        results_file = self.base_output_dir / "phase2b_results.json"
        with open(results_file, 'w') as f:
            json.dump(self.results, f, indent=2)
    
    def generate_summary_report(self):
        """Generate summary report"""
        report_file = self.base_output_dir / "phase2b_summary_report.md"
        
        with open(report_file, 'w') as f:
            f.write("# Phase 2B Additional Runs - Summary Report\n\n")
            f.write(f"**Generated**: {datetime.now().isoformat()}\n")
            f.write(f"**Total Duration**: {self.results['total_duration']/3600:.1f} hours\n\n")
            
            f.write("## Overall Results\n\n")
            f.write(f"- **Total Runs**: {self.results['total_runs']}\n")
            f.write(f"- **Successful**: {self.results['completed_runs']}\n")
            f.write(f"- **Failed**: {self.results['failed_runs']}\n")
            f.write(f"- **Success Rate**: {self.results['completed_runs']/self.results['total_runs']*100:.1f}%\n\n")
            
            f.write("## Scenario Results\n\n")
            for scenario_name, scenario_data in self.results["scenarios"].items():
                f.write(f"### {scenario_name}\n")
                f.write(f"- **Description**: {scenario_data['description']}\n")
                f.write(f"- **Successful**: {scenario_data['completed_runs']}/{scenario_data['total_runs']}\n")
                f.write(f"- **Success Rate**: {scenario_data['completed_runs']/scenario_data['total_runs']*100:.1f}%\n\n")
            
            f.write("## Next Steps\n\n")
            f.write("1. Run analysis: `python scripts/analyze_additional_results.py`\n")
            f.write("2. Compare with Phase 2A results\n")
            f.write("3. Generate publication figures\n")
            f.write("4. Proceed to Phase 3 (Paper Writing)\n")
        
        self.logger.info(f"📄 Summary report saved: {report_file}")

def main():
    parser = argparse.ArgumentParser(description="Phase 2B Additional Runs")
    parser.add_argument("--output-dir", default="results/phase2b_additional",
                       help="Base output directory")
    parser.add_argument("--scenario", choices=["miller_urey_extended", "hydrothermal_extended", "formamide_extended"],
                       help="Run only specific scenario")
    parser.add_argument("--dry-run", action="store_true",
                       help="Show what would be run without executing")
    
    args = parser.parse_args()
    
    runner = Phase2BRunner(args.output_dir)
    
    if args.dry_run:
        print("🔍 DRY RUN - What would be executed:")
        for scenario_name, config_info in runner.configs.items():
            print(f"  {scenario_name}: {config_info['runs']} runs, seeds {config_info['seeds']}")
        return
    
    if args.scenario:
        # Run only specific scenario
        if args.scenario in runner.configs:
            runner.run_scenario(args.scenario, runner.configs[args.scenario])
        else:
            print(f"❌ Unknown scenario: {args.scenario}")
            return
    else:
        # Run all scenarios
        runner.run_all_simulations()

if __name__ == "__main__":
    main()
