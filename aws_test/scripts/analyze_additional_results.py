#!/usr/bin/env python3
"""
Phase 2B Results Analyzer
========================

Analyzes results from Phase 2B additional runs and compares with Phase 2A.
Generates comprehensive reports and figures for publication.
"""

import os
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import statistics

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

class Phase2BAnalyzer:
    def __init__(self, phase2b_dir="results/phase2b_additional", phase2a_dir="aws_test"):
        self.phase2b_dir = Path(phase2b_dir)
        self.phase2a_dir = Path(phase2a_dir)
        
        self.analysis_results = {
            "timestamp": datetime.now().isoformat(),
            "phase2a_results": {},
            "phase2b_results": {},
            "comparison": {},
            "recommendations": []
        }
        
    def load_phase2a_results(self):
        """Load Phase 2A results from AWS test"""
        print("📊 Loading Phase 2A results...")
        
        phase2a_data = {
            "scenarios": defaultdict(list),
            "total_runs": 0,
            "successful_runs": 0,
            "failed_runs": 0
        }
        
        # Load from AWS test directories
        aws_dirs = ["results_16_completed", "results_28_completed", "results_all_completed"]
        
        for aws_dir in aws_dirs:
            results_dir = self.phase2a_dir / aws_dir / "results"
            if not results_dir.exists():
                continue
                
            for scenario_dir in results_dir.iterdir():
                if not scenario_dir.is_dir():
                    continue
                    
                scenario_name = scenario_dir.name
                
                for run_dir in scenario_dir.iterdir():
                    if not run_dir.is_dir():
                        continue
                        
                    phase2a_data["total_runs"] += 1
                    
                    summary_file = run_dir / "summary.txt"
                    results_file = run_dir / "results.json"
                    molecules_file = run_dir / "molecules.json"
                    
                    if summary_file.exists() and results_file.exists():
                        phase2a_data["successful_runs"] += 1
                        
                        with open(results_file, 'r') as f:
                            results_json = json.load(f)
                        
                        molecules_data = []
                        if molecules_file.exists():
                            with open(molecules_file, 'r') as f:
                                molecules_data = json.load(f)
                        
                        run_data = {
                            "run_id": run_dir.name,
                            "scenario": scenario_name,
                            "results": results_json,
                            "molecules": molecules_data,
                            "molecules_count": len(molecules_data),
                            "final_particles": results_json.get("final_state", {}).get("n_particles", 0),
                            "steps": results_json.get("final_state", {}).get("step", 0),
                            "time": results_json.get("final_state", {}).get("time", 0)
                        }
                        
                        phase2a_data["scenarios"][scenario_name].append(run_data)
                    else:
                        phase2a_data["failed_runs"] += 1
        
        self.analysis_results["phase2a_results"] = phase2a_data
        print(f"  ✅ Loaded {phase2a_data['successful_runs']}/{phase2a_data['total_runs']} Phase 2A runs")
        
    def load_phase2b_results(self):
        """Load Phase 2B results"""
        print("📊 Loading Phase 2B results...")
        
        phase2b_data = {
            "scenarios": defaultdict(list),
            "total_runs": 0,
            "successful_runs": 0,
            "failed_runs": 0
        }
        
        scenarios = ["miller_urey_extended", "hydrothermal_extended", "formamide_extended"]
        
        for scenario in scenarios:
            scenario_dir = self.phase2b_dir / scenario
            if not scenario_dir.exists():
                continue
                
            for run_dir in scenario_dir.iterdir():
                if not run_dir.is_dir():
                    continue
                    
                phase2b_data["total_runs"] += 1
                
                summary_file = run_dir / "summary.txt"
                results_file = run_dir / "results.json"
                molecules_file = run_dir / "molecules.json"
                
                if summary_file.exists() and results_file.exists():
                    phase2b_data["successful_runs"] += 1
                    
                    with open(results_file, 'r') as f:
                        results_json = json.load(f)
                    
                    molecules_data = []
                    if molecules_file.exists():
                        with open(molecules_file, 'r') as f:
                            molecules_data = json.load(f)
                    
                    run_data = {
                        "run_id": run_dir.name,
                        "scenario": scenario,
                        "results": results_json,
                        "molecules": molecules_data,
                        "molecules_count": len(molecules_data),
                        "final_particles": results_json.get("final_state", {}).get("n_particles", 0),
                        "steps": results_json.get("final_state", {}).get("step", 0),
                        "time": results_json.get("final_state", {}).get("time", 0)
                    }
                    
                    phase2b_data["scenarios"][scenario].append(run_data)
                else:
                    phase2b_data["failed_runs"] += 1
        
        self.analysis_results["phase2b_results"] = phase2b_data
        print(f"  ✅ Loaded {phase2b_data['successful_runs']}/{phase2b_data['total_runs']} Phase 2B runs")
        
    def analyze_molecular_diversity(self):
        """Analyze molecular diversity across phases"""
        print("🧬 Analyzing molecular diversity...")
        
        diversity_analysis = {
            "phase2a": {},
            "phase2b": {},
            "comparison": {}
        }
        
        # Phase 2A analysis
        phase2a_data = self.analysis_results["phase2a_results"]
        for scenario_name, runs in phase2a_data["scenarios"].items():
            all_molecules = set()
            molecules_per_run = []
            
            for run in runs:
                run_molecules = set()
                for mol in run["molecules"]:
                    mol_id = mol.get("id", "unknown")
                    all_molecules.add(mol_id)
                    run_molecules.add(mol_id)
                
                molecules_per_run.append(len(run_molecules))
            
            diversity_analysis["phase2a"][scenario_name] = {
                "total_unique_molecules": len(all_molecules),
                "average_per_run": statistics.mean(molecules_per_run) if molecules_per_run else 0,
                "std_per_run": statistics.stdev(molecules_per_run) if len(molecules_per_run) > 1 else 0,
                "runs": len(runs)
            }
        
        # Phase 2B analysis
        phase2b_data = self.analysis_results["phase2b_results"]
        for scenario_name, runs in phase2b_data["scenarios"].items():
            all_molecules = set()
            molecules_per_run = []
            
            for run in runs:
                run_molecules = set()
                for mol in run["molecules"]:
                    mol_id = mol.get("id", "unknown")
                    all_molecules.add(mol_id)
                    run_molecules.add(mol_id)
                
                molecules_per_run.append(len(run_molecules))
            
            diversity_analysis["phase2b"][scenario_name] = {
                "total_unique_molecules": len(all_molecules),
                "average_per_run": statistics.mean(molecules_per_run) if molecules_per_run else 0,
                "std_per_run": statistics.stdev(molecules_per_run) if len(molecules_per_run) > 1 else 0,
                "runs": len(runs)
            }
        
        # Comparison
        for scenario_name in diversity_analysis["phase2a"].keys():
            if scenario_name in diversity_analysis["phase2b"]:
                phase2a_mol = diversity_analysis["phase2a"][scenario_name]["total_unique_molecules"]
                phase2b_mol = diversity_analysis["phase2b"][scenario_name]["total_unique_molecules"]
                
                diversity_analysis["comparison"][scenario_name] = {
                    "phase2a_molecules": phase2a_mol,
                    "phase2b_molecules": phase2b_mol,
                    "improvement": phase2b_mol - phase2a_mol,
                    "improvement_percent": ((phase2b_mol - phase2a_mol) / phase2a_mol * 100) if phase2a_mol > 0 else 0
                }
        
        self.analysis_results["molecular_diversity"] = diversity_analysis
        
    def analyze_autocatalytic_cycles(self):
        """Analyze autocatalytic cycles (placeholder for now)"""
        print("🔄 Analyzing autocatalytic cycles...")
        
        # This would integrate with autocatalytic_detector.py
        # For now, placeholder analysis
        
        cycle_analysis = {
            "phase2a": {"total_cycles": 0, "cycles_per_scenario": {}},
            "phase2b": {"total_cycles": 0, "cycles_per_scenario": {}},
            "comparison": {}
        }
        
        # Placeholder - would need to run autocatalytic detector
        self.analysis_results["autocatalytic_cycles"] = cycle_analysis
        
    def generate_recommendations(self):
        """Generate recommendations based on analysis"""
        print("💡 Generating recommendations...")
        
        recommendations = []
        
        # Check molecular diversity
        diversity = self.analysis_results.get("molecular_diversity", {})
        phase2b_total = sum(
            scenario["total_unique_molecules"] 
            for scenario in diversity.get("phase2b", {}).values()
        )
        
        if phase2b_total >= 100:
            recommendations.append("✅ SUCCESS: Molecular diversity target (100+) achieved")
        elif phase2b_total >= 50:
            recommendations.append("⚠️ PARTIAL: Molecular diversity improved but below target")
        else:
            recommendations.append("❌ INSUFFICIENT: Molecular diversity still too low")
        
        # Check formamide activity
        formamide_2b = diversity.get("phase2b", {}).get("formamide_extended", {})
        if formamide_2b.get("total_unique_molecules", 0) > 0:
            recommendations.append("✅ SUCCESS: Formamide scenario now active")
        else:
            recommendations.append("❌ ISSUE: Formamide scenario still inactive")
        
        # Check completion rates
        phase2b_data = self.analysis_results["phase2b_results"]
        completion_rate = phase2b_data["successful_runs"] / phase2b_data["total_runs"] * 100 if phase2b_data["total_runs"] > 0 else 0
        
        if completion_rate >= 90:
            recommendations.append("✅ SUCCESS: High completion rate achieved")
        elif completion_rate >= 80:
            recommendations.append("⚠️ ACCEPTABLE: Completion rate acceptable")
        else:
            recommendations.append("❌ ISSUE: Low completion rate")
        
        # Overall recommendation
        if phase2b_total >= 100 and completion_rate >= 90:
            recommendations.append("🎉 PHASE 2 COMPLETE: Ready for Phase 3 (Paper Writing)")
        elif phase2b_total >= 50 and completion_rate >= 80:
            recommendations.append("⚠️ PHASE 2 MOSTLY COMPLETE: Consider additional runs")
        else:
            recommendations.append("❌ PHASE 2 INCOMPLETE: Need more work")
        
        self.analysis_results["recommendations"] = recommendations
        
    def generate_report(self):
        """Generate comprehensive analysis report"""
        print("📄 Generating analysis report...")
        
        report_file = self.phase2b_dir / "phase2b_analysis_report.md"
        
        with open(report_file, 'w') as f:
            f.write("# Phase 2B Analysis Report\n\n")
            f.write(f"**Generated**: {datetime.now().isoformat()}\n\n")
            
            # Executive Summary
            f.write("## Executive Summary\n\n")
            
            phase2a_data = self.analysis_results["phase2a_results"]
            phase2b_data = self.analysis_results["phase2b_results"]
            
            f.write(f"- **Phase 2A**: {phase2a_data['successful_runs']}/{phase2a_data['total_runs']} runs successful\n")
            f.write(f"- **Phase 2B**: {phase2b_data['successful_runs']}/{phase2b_data['total_runs']} runs successful\n")
            
            diversity = self.analysis_results.get("molecular_diversity", {})
            phase2a_total = sum(
                scenario["total_unique_molecules"] 
                for scenario in diversity.get("phase2a", {}).values()
            )
            phase2b_total = sum(
                scenario["total_unique_molecules"] 
                for scenario in diversity.get("phase2b", {}).values()
            )
            
            f.write(f"- **Phase 2A Molecules**: {phase2a_total}\n")
            f.write(f"- **Phase 2B Molecules**: {phase2b_total}\n")
            f.write(f"- **Improvement**: {phase2b_total - phase2a_total} molecules\n\n")
            
            # Detailed Analysis
            f.write("## Detailed Analysis\n\n")
            
            f.write("### Molecular Diversity\n\n")
            for scenario_name in diversity.get("phase2a", {}).keys():
                f.write(f"#### {scenario_name.replace('_', ' ').title()}\n\n")
                
                phase2a_mol = diversity["phase2a"][scenario_name]["total_unique_molecules"]
                phase2b_mol = diversity["phase2b"].get(scenario_name.replace("_extended", ""), {}).get("total_unique_molecules", 0)
                
                f.write(f"- **Phase 2A**: {phase2a_mol} molecules\n")
                f.write(f"- **Phase 2B**: {phase2b_mol} molecules\n")
                f.write(f"- **Improvement**: {phase2b_mol - phase2a_mol} molecules\n\n")
            
            # Recommendations
            f.write("## Recommendations\n\n")
            for rec in self.analysis_results["recommendations"]:
                f.write(f"- {rec}\n")
            
            f.write("\n## Next Steps\n\n")
            f.write("1. Review analysis results\n")
            f.write("2. Generate publication figures\n")
            f.write("3. Proceed to Phase 3 (Paper Writing)\n")
        
        print(f"📄 Analysis report saved: {report_file}")
        
        # Save JSON results
        json_file = self.phase2b_dir / "phase2b_analysis_results.json"
        with open(json_file, 'w') as f:
            json.dump(self.analysis_results, f, indent=2)
        
        print(f"📊 Analysis results saved: {json_file}")
        
    def run_analysis(self):
        """Run complete analysis"""
        print("🔍 Starting Phase 2B Analysis")
        print("=" * 50)
        
        # Load data
        self.load_phase2a_results()
        self.load_phase2b_results()
        
        # Run analyses
        self.analyze_molecular_diversity()
        self.analyze_autocatalytic_cycles()
        self.generate_recommendations()
        
        # Generate report
        self.generate_report()
        
        print("=" * 50)
        print("✅ Analysis complete!")
        
        # Print summary
        print("\n📊 SUMMARY:")
        for rec in self.analysis_results["recommendations"]:
            print(f"  {rec}")

def main():
    parser = argparse.ArgumentParser(description="Phase 2B Results Analyzer")
    parser.add_argument("--phase2b-dir", default="results/phase2b_additional",
                       help="Phase 2B results directory")
    parser.add_argument("--phase2a-dir", default="aws_test",
                       help="Phase 2A results directory")
    
    args = parser.parse_args()
    
    analyzer = Phase2BAnalyzer(args.phase2b_dir, args.phase2a_dir)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
