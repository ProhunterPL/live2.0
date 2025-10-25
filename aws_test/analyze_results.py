#!/usr/bin/env python3
"""
AWS Test Results Analysis Script
Analyzes simulation results from all three test directories
"""

import json
import os
from pathlib import Path
from collections import defaultdict
import statistics

def analyze_results_directory(results_dir):
    """Analyze results from a specific directory"""
    results_data = {
        'scenarios': defaultdict(list),
        'total_runs': 0,
        'successful_runs': 0,
        'failed_runs': 0
    }
    
    for scenario_dir in results_dir.iterdir():
        if not scenario_dir.is_dir():
            continue
            
        scenario_name = scenario_dir.name
        print(f"\nAnalyzing {scenario_name} scenario...")
        
        for run_dir in scenario_dir.iterdir():
            if not run_dir.is_dir():
                continue
                
            results_data['total_runs'] += 1
            
            # Check if run completed successfully
            summary_file = run_dir / 'summary.txt'
            results_file = run_dir / 'results.json'
            molecules_file = run_dir / 'molecules.json'
            
            if summary_file.exists() and results_file.exists():
                results_data['successful_runs'] += 1
                
                # Parse summary
                with open(summary_file, 'r') as f:
                    summary_content = f.read()
                
                # Parse results.json
                with open(results_file, 'r') as f:
                    results_json = json.load(f)
                
                # Parse molecules.json
                molecules_data = []
                if molecules_file.exists():
                    with open(molecules_file, 'r') as f:
                        molecules_data = json.load(f)
                
                run_data = {
                    'run_id': run_dir.name,
                    'summary': summary_content,
                    'results': results_json,
                    'molecules': molecules_data,
                    'molecules_count': len(molecules_data),
                    'final_particles': results_json.get('final_state', {}).get('n_particles', 0),
                    'steps': results_json.get('final_state', {}).get('step', 0),
                    'time': results_json.get('final_state', {}).get('time', 0)
                }
                
                results_data['scenarios'][scenario_name].append(run_data)
            else:
                results_data['failed_runs'] += 1
                print(f"  Failed run: {run_dir.name}")
    
    return results_data

def print_scenario_summary(scenario_name, runs):
    """Print summary for a specific scenario"""
    print(f"\n{'='*60}")
    print(f"SCENARIO: {scenario_name.upper()}")
    print(f"{'='*60}")
    
    if not runs:
        print("No successful runs found.")
        return
    
    # Basic statistics
    total_runs = len(runs)
    final_particles = [run['final_particles'] for run in runs]
    molecules_counts = [run['molecules_count'] for run in runs]
    steps = [run['steps'] for run in runs]
    
    print(f"Total runs: {total_runs}")
    print(f"Average final particles: {statistics.mean(final_particles):.0f}")
    print(f"Particle range: {min(final_particles)} - {max(final_particles)}")
    print(f"Average molecules detected: {statistics.mean(molecules_counts):.1f}")
    print(f"Molecules range: {min(molecules_counts)} - {max(molecules_counts)}")
    print(f"Steps: {min(steps)} - {max(steps)}")
    
    # Show configuration from first run
    if runs:
        config = runs[0]['results'].get('configuration', {})
        print(f"\nConfiguration:")
        print(f"  Temperature: {config.get('temperature', 'N/A')}K")
        print(f"  Initial molecules: {len(config.get('initial_molecules', []))}")
        for mol in config.get('initial_molecules', []):
            print(f"    - {mol.get('name', 'unknown')} ({mol.get('formula', 'unknown')}): {mol.get('count', 0)}")
    
    # Show molecules detected across all runs
    all_molecules = defaultdict(int)
    for run in runs:
        for mol in run['molecules']:
            all_molecules[mol['id']] += 1
    
    if all_molecules:
        print(f"\nMolecules detected across runs:")
        for mol_id, count in sorted(all_molecules.items(), key=lambda x: x[1], reverse=True):
            print(f"  {mol_id}: {count}/{total_runs} runs")
    
    # Show individual run details
    print(f"\nIndividual run details:")
    for run in runs:
        print(f"  {run['run_id']}: {run['final_particles']} particles, {run['molecules_count']} molecules")

def main():
    """Main analysis function"""
    print("AWS Test Results Analysis")
    print("=" * 60)
    
    # Analyze all three result directories
    directories = [
        ('results_16_completed', '16-hour runs'),
        ('results_28_completed', '28-hour runs'), 
        ('results_all_completed', 'All completed runs')
    ]
    
    all_results = {}
    
    for dir_name, description in directories:
        results_dir = Path(f"{dir_name}/results")
        if results_dir.exists():
            print(f"\n{description.upper()}")
            print("-" * 40)
            results_data = analyze_results_directory(results_dir)
            all_results[dir_name] = results_data
            
            print(f"Total runs: {results_data['total_runs']}")
            print(f"Successful: {results_data['successful_runs']}")
            print(f"Failed: {results_data['failed_runs']}")
            
            # Print scenario summaries
            for scenario_name, runs in results_data['scenarios'].items():
                print_scenario_summary(scenario_name, runs)
        else:
            print(f"Directory not found: {results_dir}")
    
    # Overall summary
    print(f"\n{'='*60}")
    print("OVERALL SUMMARY")
    print(f"{'='*60}")
    
    total_runs_all = sum(data['total_runs'] for data in all_results.values())
    total_successful_all = sum(data['successful_runs'] for data in all_results.values())
    total_failed_all = sum(data['failed_runs'] for data in all_results.values())
    
    print(f"Total runs across all directories: {total_runs_all}")
    print(f"Total successful: {total_successful_all}")
    print(f"Total failed: {total_failed_all}")
    print(f"Success rate: {(total_successful_all/total_runs_all*100):.1f}%" if total_runs_all > 0 else "N/A")
    
    # Scenario comparison
    print(f"\nScenario comparison:")
    scenario_stats = defaultdict(lambda: {'runs': 0, 'particles': [], 'molecules': []})
    
    for dir_name, data in all_results.items():
        for scenario_name, runs in data['scenarios'].items():
            scenario_stats[scenario_name]['runs'] += len(runs)
            scenario_stats[scenario_name]['particles'].extend([run['final_particles'] for run in runs])
            scenario_stats[scenario_name]['molecules'].extend([run['molecules_count'] for run in runs])
    
    for scenario_name, stats in scenario_stats.items():
        if stats['particles']:
            avg_particles = statistics.mean(stats['particles'])
            avg_molecules = statistics.mean(stats['molecules'])
            print(f"  {scenario_name}: {stats['runs']} runs, avg {avg_particles:.0f} particles, avg {avg_molecules:.1f} molecules")

if __name__ == "__main__":
    main()
