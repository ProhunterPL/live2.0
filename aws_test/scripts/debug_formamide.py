#!/usr/bin/env python3
"""
Formamide Debug Tool
====================

Debug tool for formamide scenario to identify why it produces 0 molecules.
Runs short simulations with detailed logging and analysis.
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
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

class FormamideDebugger:
    def __init__(self, output_dir="results/phase2b_additional/formamide_debug"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self.setup_logging()
        
        # Debug configurations
        self.debug_configs = {
            "short_test": {
                "steps": 10000,
                "description": "10K steps - quick test",
                "seeds": [200, 201, 202]
            },
            "medium_test": {
                "steps": 50000,
                "description": "50K steps - medium test", 
                "seeds": [203, 204, 205]
            },
            "long_test": {
                "steps": 100000,
                "description": "100K steps - extended test",
                "seeds": [206, 207, 208]
            }
        }
        
        self.results = {
            "start_time": datetime.now().isoformat(),
            "debug_tests": {}
        }
        
    def setup_logging(self):
        """Setup logging for debugging"""
        log_dir = self.output_dir / "logs"
        log_dir.mkdir(exist_ok=True)
        
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_dir / "formamide_debug.log"),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def run_debug_test(self, test_name, config):
        """Run a debug test"""
        self.logger.info(f"🔍 Running {test_name}: {config['description']}")
        
        test_results = {
            "test_name": test_name,
            "description": config["description"],
            "steps": config["steps"],
            "runs": []
        }
        
        for i, seed in enumerate(config["seeds"]):
            run_id = i + 1
            self.logger.info(f"  Run {run_id}/3 (seed {seed})")
            
            # Create output directory
            run_dir = self.output_dir / test_name / f"run_{run_id}"
            run_dir.mkdir(parents=True, exist_ok=True)
            
            # Prepare command
            cmd = [
                "python", "scripts/run_phase2_full.py",
                "--config", "configs/phase2_formamide_debug.yaml",
                "--output", str(run_dir),
                "--seed", str(seed),
                "--steps", str(config["steps"])
            ]
            
            # Run simulation
            start_time = time.time()
            try:
                result = subprocess.run(
                    cmd,
                    cwd=project_root,
                    capture_output=True,
                    text=True,
                    timeout=1800  # 30 minute timeout
                )
                
                duration = time.time() - start_time
                
                # Analyze results
                analysis = self.analyze_run_results(run_dir)
                
                run_result = {
                    "run_id": run_id,
                    "seed": seed,
                    "duration": duration,
                    "return_code": result.returncode,
                    "stdout": result.stdout,
                    "stderr": result.stderr,
                    "analysis": analysis
                }
                
                test_results["runs"].append(run_result)
                
                if result.returncode == 0:
                    self.logger.info(f"    ✅ Run {run_id} completed in {duration:.1f}s")
                    self.logger.info(f"    📊 Molecules detected: {analysis['molecules_count']}")
                else:
                    self.logger.error(f"    ❌ Run {run_id} failed: {result.stderr}")
                    
            except subprocess.TimeoutExpired:
                self.logger.error(f"    ⏰ Run {run_id} timed out")
                test_results["runs"].append({
                    "run_id": run_id,
                    "seed": seed,
                    "status": "timeout",
                    "error": "Simulation timed out"
                })
            except Exception as e:
                self.logger.error(f"    💥 Run {run_id} crashed: {str(e)}")
                test_results["runs"].append({
                    "run_id": run_id,
                    "seed": seed,
                    "status": "crashed",
                    "error": str(e)
                })
        
        self.results["debug_tests"][test_name] = test_results
        return test_results
    
    def analyze_run_results(self, run_dir):
        """Analyze results from a single run"""
        analysis = {
            "molecules_count": 0,
            "reactions_count": 0,
            "final_particles": 0,
            "energy_drift": 0,
            "temperature": 0,
            "files_present": [],
            "issues": []
        }
        
        # Check if files exist
        expected_files = ["summary.txt", "results.json", "molecules.json", "simulation.log"]
        for file_name in expected_files:
            file_path = run_dir / file_name
            if file_path.exists():
                analysis["files_present"].append(file_name)
            else:
                analysis["issues"].append(f"Missing file: {file_name}")
        
        # Analyze results.json
        results_file = run_dir / "results.json"
        if results_file.exists():
            try:
                with open(results_file, 'r') as f:
                    results_data = json.load(f)
                
                analysis["final_particles"] = results_data.get("final_state", {}).get("n_particles", 0)
                analysis["molecules_count"] = len(results_data.get("molecules_detected", []))
                analysis["reactions_count"] = len(results_data.get("reactions_observed", []))
                
            except Exception as e:
                analysis["issues"].append(f"Error reading results.json: {str(e)}")
        
        # Analyze molecules.json
        molecules_file = run_dir / "molecules.json"
        if molecules_file.exists():
            try:
                with open(molecules_file, 'r') as f:
                    molecules_data = json.load(f)
                
                analysis["molecules_count"] = len(molecules_data)
                
            except Exception as e:
                analysis["issues"].append(f"Error reading molecules.json: {str(e)}")
        
        # Analyze simulation.log
        log_file = run_dir / "simulation.log"
        if log_file.exists():
            try:
                with open(log_file, 'r') as f:
                    log_content = f.read()
                
                # Look for key metrics in log
                if "Energy drift" in log_content:
                    # Extract energy drift value
                    lines = log_content.split('\n')
                    for line in lines:
                        if "Energy drift" in line:
                            try:
                                drift_value = float(line.split(":")[-1].strip())
                                analysis["energy_drift"] = drift_value
                                break
                            except:
                                pass
                
                if "Temperature" in log_content:
                    # Extract temperature value
                    lines = log_content.split('\n')
                    for line in lines:
                        if "Temperature" in line and ":" in line:
                            try:
                                temp_value = float(line.split(":")[-1].strip())
                                analysis["temperature"] = temp_value
                                break
                            except:
                                pass
                
            except Exception as e:
                analysis["issues"].append(f"Error reading simulation.log: {str(e)}")
        
        return analysis
    
    def run_all_debug_tests(self):
        """Run all debug tests"""
        self.logger.info("🔍 Starting Formamide Debug Analysis")
        
        for test_name, config in self.debug_configs.items():
            self.run_debug_test(test_name, config)
        
        # Generate debug report
        self.generate_debug_report()
        
    def generate_debug_report(self):
        """Generate debug report"""
        report_file = self.output_dir / "formamide_debug_report.md"
        
        with open(report_file, 'w') as f:
            f.write("# Formamide Debug Report\n\n")
            f.write(f"**Generated**: {datetime.now().isoformat()}\n\n")
            
            f.write("## Summary\n\n")
            
            total_molecules = 0
            total_runs = 0
            
            for test_name, test_data in self.results["debug_tests"].items():
                test_molecules = sum(run.get("analysis", {}).get("molecules_count", 0) for run in test_data["runs"])
                test_runs = len(test_data["runs"])
                
                total_molecules += test_molecules
                total_runs += test_runs
                
                f.write(f"### {test_name}\n")
                f.write(f"- **Steps**: {test_data['steps']}\n")
                f.write(f"- **Runs**: {test_runs}\n")
                f.write(f"- **Molecules detected**: {test_molecules}\n")
                f.write(f"- **Average per run**: {test_molecules/test_runs if test_runs > 0 else 0:.1f}\n\n")
            
            f.write(f"## Overall Results\n\n")
            f.write(f"- **Total runs**: {total_runs}\n")
            f.write(f"- **Total molecules**: {total_molecules}\n")
            f.write(f"- **Average per run**: {total_molecules/total_runs if total_runs > 0 else 0:.1f}\n\n")
            
            f.write("## Detailed Analysis\n\n")
            
            for test_name, test_data in self.results["debug_tests"].items():
                f.write(f"### {test_name}\n\n")
                
                for run in test_data["runs"]:
                    f.write(f"#### Run {run['run_id']} (seed {run['seed']})\n")
                    
                    if "analysis" in run:
                        analysis = run["analysis"]
                        f.write(f"- **Molecules**: {analysis['molecules_count']}\n")
                        f.write(f"- **Reactions**: {analysis['reactions_count']}\n")
                        f.write(f"- **Final particles**: {analysis['final_particles']}\n")
                        f.write(f"- **Energy drift**: {analysis['energy_drift']}\n")
                        f.write(f"- **Temperature**: {analysis['temperature']}\n")
                        
                        if analysis["issues"]:
                            f.write(f"- **Issues**: {', '.join(analysis['issues'])}\n")
                    
                    f.write("\n")
            
            f.write("## Recommendations\n\n")
            
            if total_molecules == 0:
                f.write("❌ **CRITICAL ISSUE**: No molecules detected in any run\n\n")
                f.write("**Possible causes**:\n")
                f.write("1. **Configuration error**: Check formamide_debug.yaml\n")
                f.write("2. **Physics parameters**: Bond formation thresholds too strict\n")
                f.write("3. **Energy input**: UV radiation not effective\n")
                f.write("4. **Catalysis**: Mineral catalysts not working\n")
                f.write("5. **Molecule detection**: MoleculeExtractor not working\n\n")
                
                f.write("**Next steps**:\n")
                f.write("1. Check configuration files\n")
                f.write("2. Verify physics parameters\n")
                f.write("3. Test with simpler molecules first\n")
                f.write("4. Check MoleculeExtractor functionality\n")
            elif total_molecules < 5:
                f.write("⚠️ **LOW ACTIVITY**: Very few molecules detected\n\n")
                f.write("**Recommendations**:\n")
                f.write("1. Increase simulation time\n")
                f.write("2. Adjust physics parameters\n")
                f.write("3. Enhance energy input\n")
                f.write("4. Add more catalysts\n")
            else:
                f.write("✅ **GOOD ACTIVITY**: Molecules detected\n\n")
                f.write("**Recommendations**:\n")
                f.write("1. Proceed with extended simulations\n")
                f.write("2. Optimize parameters for more molecules\n")
                f.write("3. Focus on autocatalytic cycles\n")
        
        self.logger.info(f"📄 Debug report saved: {report_file}")
        
        # Save JSON results
        json_file = self.output_dir / "formamide_debug_results.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        self.logger.info(f"📊 Debug results saved: {json_file}")

def main():
    parser = argparse.ArgumentParser(description="Formamide Debug Tool")
    parser.add_argument("--output-dir", default="results/phase2b_additional/formamide_debug",
                       help="Output directory for debug results")
    parser.add_argument("--test", choices=["short_test", "medium_test", "long_test"],
                       help="Run only specific test")
    
    args = parser.parse_args()
    
    debugger = FormamideDebugger(args.output_dir)
    
    if args.test:
        # Run only specific test
        if args.test in debugger.debug_configs:
            debugger.run_debug_test(args.test, debugger.debug_configs[args.test])
            debugger.generate_debug_report()
        else:
            print(f"❌ Unknown test: {args.test}")
            return
    else:
        # Run all debug tests
        debugger.run_all_debug_tests()

if __name__ == "__main__":
    main()
