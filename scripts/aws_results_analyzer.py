#!/usr/bin/env python3
"""
AWS Results Batch Analyzer
===========================

Automatically analyze all downloaded AWS simulation results.
Generates comprehensive report with statistics, figures, and data ready for paper.

Usage:
    python scripts/aws_results_analyzer.py --input results/aws_batch
    
Author: Live 2.0 Team
Date: October 2025
"""

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import subprocess

class AWSBatchAnalyzer:
    """Analyze batch of AWS simulation results"""
    
    def __init__(self, results_dir):
        self.results_dir = Path(results_dir)
        self.analysis_dir = self.results_dir / 'analysis'
        self.analysis_dir.mkdir(exist_ok=True)
        
        self.report = {
            'timestamp': datetime.now().isoformat(),
            'scenarios': {},
            'summary': {},
            'figures': [],
            'errors': []
        }
    
    def discover_simulations(self):
        """Find all completed simulations"""
        print("🔍 Discovering simulations...")
        
        scenarios = defaultdict(list)
        
        for scenario_dir in self.results_dir.iterdir():
            if not scenario_dir.is_dir() or scenario_dir.name == 'analysis':
                continue
            
            scenario = scenario_dir.name
            
            for run_dir in sorted(scenario_dir.iterdir()):
                if not run_dir.is_dir():
                    continue
                
                # Check if complete
                summary_file = run_dir / 'summary.txt'
                if summary_file.exists():
                    scenarios[scenario].append({
                        'path': run_dir,
                        'run_id': run_dir.name
                    })
        
        for scenario, runs in scenarios.items():
            print(f"  {scenario}: {len(runs)} runs")
        
        return scenarios
    
    def extract_molecules(self, scenarios):
        """Extract molecules from all simulations"""
        print("\n📊 Extracting molecules...")
        
        for scenario, runs in scenarios.items():
            print(f"\n{scenario}:")
            scenario_analysis = self.analysis_dir / scenario
            scenario_analysis.mkdir(exist_ok=True)
            
            for i, run in enumerate(runs, 1):
                print(f"  [{i}/{len(runs)}] {run['run_id']}...", end=' ')
                
                # Run quick_analyze.py
                cmd = [
                    'python', 'scripts/quick_analyze.py',
                    str(run['path']),
                    '--output', str(scenario_analysis / f"{run['run_id']}_molecules.json")
                ]
                
                try:
                    result = subprocess.run(cmd, capture_output=True, 
                                          text=True, timeout=300)
                    if result.returncode == 0:
                        print("✅")
                        run['molecules_extracted'] = True
                    else:
                        print("❌")
                        run['molecules_extracted'] = False
                        self.report['errors'].append({
                            'scenario': scenario,
                            'run': run['run_id'],
                            'step': 'extraction',
                            'error': result.stderr[:200]
                        })
                except Exception as e:
                    print(f"❌ {e}")
                    run['molecules_extracted'] = False
    
    def build_reaction_networks(self, scenarios):
        """Build reaction networks for each scenario"""
        print("\n🕸️  Building reaction networks...")
        
        for scenario, runs in scenarios.items():
            print(f"\n{scenario}:")
            
            # Collect all molecule files for this scenario
            scenario_analysis = self.analysis_dir / scenario
            molecule_files = list(scenario_analysis.glob("*_molecules.json"))
            
            if not molecule_files:
                print("  ⚠️  No molecule data found, skipping")
                continue
            
            print(f"  Merging {len(molecule_files)} runs...")
            
            # Run reaction_network_analyzer.py with merge
            network_output = scenario_analysis / 'network'
            network_output.mkdir(exist_ok=True)
            
            # Use first run's data as base, then merge others
            # TODO: This should actually merge multiple snapshots
            # For now, just analyze the merged molecules
            cmd = [
                'python', 'scripts/reaction_network_analyzer.py',
                str(scenario_analysis),
                '--merge',
                '--output', str(network_output)
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, 
                                      text=True, timeout=600)
                if result.returncode == 0:
                    print("  ✅ Network built")
                    
                    # Save network info
                    network_file = network_output / 'reaction_network.json'
                    if network_file.exists():
                        self.report['scenarios'][scenario] = {
                            'network_file': str(network_file)
                        }
                else:
                    print(f"  ❌ Error: {result.stderr[:200]}")
            except Exception as e:
                print(f"  ❌ Error: {e}")
    
    def detect_autocatalytic_cycles(self, scenarios):
        """Detect autocatalytic cycles in each scenario"""
        print("\n🔄 Detecting autocatalytic cycles...")
        
        for scenario, runs in scenarios.items():
            print(f"\n{scenario}:")
            
            network_file = self.analysis_dir / scenario / 'network' / 'reaction_network.json'
            
            if not network_file.exists():
                print("  ⚠️  No network file, skipping")
                continue
            
            # Run autocatalytic_detector.py
            cycles_output = self.analysis_dir / scenario / 'cycles'
            cycles_output.mkdir(exist_ok=True)
            
            cmd = [
                'python', 'scripts/autocatalytic_detector.py',
                str(network_file),
                '--output', str(cycles_output)
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, 
                                      text=True, timeout=300)
                if result.returncode == 0:
                    # Parse output for cycle count
                    if 'cycles detected' in result.stdout:
                        # Extract number
                        import re
                        match = re.search(r'(\d+) cycles detected', result.stdout)
                        if match:
                            cycle_count = int(match.group(1))
                            print(f"  ✅ {cycle_count} cycles detected")
                            
                            if scenario in self.report['scenarios']:
                                self.report['scenarios'][scenario]['cycles'] = cycle_count
                    else:
                        print("  ✅ Analysis complete")
                else:
                    print(f"  ❌ Error: {result.stderr[:200]}")
            except Exception as e:
                print(f"  ❌ Error: {e}")
    
    def compare_scenarios(self):
        """Compare all scenarios"""
        print("\n📈 Comparing scenarios...")
        
        cmd = [
            'python', 'scripts/compare_scenarios.py',
            str(self.analysis_dir),
            '--output', str(self.analysis_dir / 'comparison')
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, 
                                  text=True, timeout=300)
            if result.returncode == 0:
                print("  ✅ Comparison complete")
                self.report['comparison'] = str(self.analysis_dir / 'comparison')
            else:
                print(f"  ❌ Error: {result.stderr[:200]}")
        except Exception as e:
            print(f"  ❌ Error: {e}")
    
    def generate_figures(self):
        """Generate all publication figures"""
        print("\n🎨 Generating figures...")
        
        # This would call plot_molecular_diversity.py, plot_reaction_networks.py, etc.
        figure_scripts = [
            'plot_molecular_diversity.py',
            'plot_reaction_networks.py',
            'plot_autocatalytic_cycles.py',
            'plot_emergence_timeline.py'
        ]
        
        figures_dir = self.analysis_dir / 'figures'
        figures_dir.mkdir(exist_ok=True)
        
        for script in figure_scripts:
            print(f"  Running {script}...", end=' ')
            
            cmd = [
                'python', f'scripts/{script}',
                str(self.analysis_dir),
                '--output', str(figures_dir)
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, 
                                      text=True, timeout=300)
                if result.returncode == 0:
                    print("✅")
                    self.report['figures'].append(script)
                else:
                    print("⚠️")
            except Exception as e:
                print(f"❌ {e}")
    
    def generate_report(self):
        """Generate comprehensive analysis report"""
        print("\n📝 Generating report...")
        
        report_file = self.analysis_dir / 'batch_analysis_report.json'
        
        # Add summary statistics
        self.report['summary'] = {
            'total_scenarios': len(self.report['scenarios']),
            'total_cycles': sum(s.get('cycles', 0) 
                              for s in self.report['scenarios'].values()),
            'figures_generated': len(self.report['figures']),
            'errors_count': len(self.report['errors'])
        }
        
        with open(report_file, 'w') as f:
            json.dump(self.report, f, indent=2)
        
        print(f"✅ Report saved: {report_file}")
        
        # Print summary
        print("\n" + "="*60)
        print("📊 ANALYSIS SUMMARY")
        print("="*60)
        print(f"Scenarios analyzed: {self.report['summary']['total_scenarios']}")
        print(f"Total autocatalytic cycles: {self.report['summary']['total_cycles']}")
        print(f"Figures generated: {self.report['summary']['figures_generated']}")
        print(f"Errors: {self.report['summary']['errors_count']}")
        
        if self.report['errors']:
            print("\n⚠️  Errors occurred:")
            for err in self.report['errors'][:5]:
                print(f"  - {err['scenario']}/{err['run']}: {err['step']}")
    
    def run_full_analysis(self):
        """Run complete analysis pipeline"""
        print("="*60)
        print("AWS BATCH ANALYSIS PIPELINE")
        print("="*60)
        
        # 1. Discover
        scenarios = self.discover_simulations()
        
        if not scenarios:
            print("\n❌ No simulations found!")
            return
        
        # 2. Extract molecules
        self.extract_molecules(scenarios)
        
        # 3. Build networks
        self.build_reaction_networks(scenarios)
        
        # 4. Detect cycles
        self.detect_autocatalytic_cycles(scenarios)
        
        # 5. Compare scenarios
        self.compare_scenarios()
        
        # 6. Generate figures
        self.generate_figures()
        
        # 7. Report
        self.generate_report()
        
        print("\n✅ Analysis complete!")
        print(f"📁 Results: {self.analysis_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze batch of AWS simulation results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--input', default='./results/aws_batch',
                       help='Directory with downloaded AWS results')
    
    args = parser.parse_args()
    
    results_dir = Path(args.input)
    
    if not results_dir.exists():
        print(f"❌ Directory not found: {results_dir}")
        sys.exit(1)
    
    analyzer = AWSBatchAnalyzer(results_dir)
    analyzer.run_full_analysis()


if __name__ == '__main__':
    main()

