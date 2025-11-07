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
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add project root to path
# __file__ is at aws_test/scripts/run_phase2b_additional.py
# project_root should be live2.0/ (two levels up from scripts/)
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

class Phase2BRunner:
    def __init__(self, base_output_dir="results/phase2b_additional", max_parallel=2):
        # Make base_output_dir relative to project root
        self.base_output_dir = Path(base_output_dir)
        if not self.base_output_dir.is_absolute():
            self.base_output_dir = project_root / self.base_output_dir
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Maximum parallel simulations (default: 2 for 64 CPU cores, ~32 cores each)
        self.max_parallel = max_parallel
        
        # Setup logging
        self.setup_logging()
        
        # Log project root for debugging
        self.logger.info(f"📁 Project root: {project_root}")
        self.logger.info(f"📁 Script path: {project_root / 'scripts' / 'run_phase2_full.py'}")
        self.logger.info(f"⚡ Max parallel simulations: {self.max_parallel}")
        
        # Simulation configurations (SUPER FAST MODE - optimized)
        # Config files are in aws_test/configs/ relative to project root
        self.configs = {
            "miller_urey_extended": {
                "config_file": str(project_root / "aws_test" / "configs" / "phase2_miller_urey_extended_SUPER_FAST.yaml"),
                "runs": 10,
                "seeds": list(range(100, 110)),
                "description": "Extended Miller-Urey (500K steps) - SUPER FAST MODE"
            },
            "hydrothermal_extended": {
                "config_file": str(project_root / "aws_test" / "configs" / "phase2_hydrothermal_extended_SUPER_FAST.yaml"),
                "runs": 10,
                "seeds": list(range(110, 120)),
                "description": "Extended Hydrothermal (500K steps) - SUPER FAST MODE"
            },
            "formamide_extended": {
                "config_file": str(project_root / "aws_test" / "configs" / "phase2_formamide_extended_SUPER_FAST.yaml"),
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
        
        # Use absolute path to script
        script_path = project_root / "scripts" / "run_phase2_full.py"
        if not script_path.exists():
            self.logger.error(f"❌ Script not found: {script_path}")
            return {
                "status": "failed",
                "duration": 0,
                "error": f"Script not found: {script_path}"
            }
        
        cmd = [
            python_cmd, str(script_path),
            "--config", config_file,
            "--output", str(output_dir),
            "--seed", str(seed),
            "--steps", "500000",  # 500K steps
            "--force-cpu"  # Force CPU mode to use all available cores on AWS
        ]
        
        # Run simulation
        start_time = time.time()
        try:
            result = subprocess.run(
                cmd,
                cwd=project_root,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hours timeout (500K steps can take 10-14 hours)
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
            self.logger.error(f"⏰ {scenario_name} run {run_id} timed out after 24 hours")
            return {
                "status": "timeout",
                "duration": 86400,
                "error": "Simulation timed out after 24 hours"
            }
        except Exception as e:
            self.logger.error(f"💥 {scenario_name} run {run_id} crashed: {str(e)}")
            return {
                "status": "crashed",
                "duration": time.time() - start_time,
                "error": str(e)
            }
    
    def run_scenario(self, scenario_name, config_info):
        """Run all simulations for a scenario (in parallel)"""
        self.logger.info(f"🚀 Starting {scenario_name}: {config_info['description']}")
        
        scenario_results = {
            "scenario": scenario_name,
            "description": config_info["description"],
            "total_runs": config_info["runs"],
            "completed_runs": 0,
            "failed_runs": 0,
            "runs": []
        }
        
        # Prepare all simulation tasks
        tasks = []
        for i in range(config_info["runs"]):
            run_id = i + 1
            seed = config_info["seeds"][i]
            tasks.append({
                "scenario_name": scenario_name,
                "run_id": run_id,
                "seed": seed,
                "config_file": config_info["config_file"]
            })
        
        # Run simulations in parallel
        self.logger.info(f"⚡ Running {len(tasks)} simulations with max {self.max_parallel} parallel")
        
        with ThreadPoolExecutor(max_workers=self.max_parallel) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(self.run_single_simulation, 
                              task["scenario_name"], 
                              task["run_id"], 
                              task["seed"], 
                              task["config_file"]): task
                for task in tasks
            }
            
            # Process completed tasks as they finish
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    result = future.result()
                    
                    scenario_results["runs"].append({
                        "run_id": task["run_id"],
                        "seed": task["seed"],
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
                    self.logger.info(f"📊 Progress: {total_progress}/30 runs completed ({scenario_results['completed_runs']}/{scenario_results['total_runs']} for {scenario_name})")
                    
                except Exception as e:
                    self.logger.error(f"💥 Exception in {scenario_name} run {task['run_id']}: {e}")
                    scenario_results["failed_runs"] += 1
                    self.results["failed_runs"] += 1
        
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
    parser.add_argument("--max-parallel", type=int, default=2,
                       help="Maximum parallel simulations (default: 2 for 64 CPU cores)")
    
    args = parser.parse_args()
    
    runner = Phase2BRunner(args.output_dir, max_parallel=args.max_parallel)
    
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
