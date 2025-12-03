#!/usr/bin/env python3
"""Calculate statistical tests (Kruskal-Wallis, Mann-Whitney) for Phase 2B results"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import json
import numpy as np
from scipy import stats
from collections import defaultdict

def load_scenario_data():
    """Load data for each scenario"""
    base_dir = Path('results/phase2b_additional')
    scenarios = {}
    
    for scenario_dir in base_dir.iterdir():
        if not scenario_dir.is_dir():
            continue
        
        scenario = scenario_dir.name
        species_counts = []
        cycle_counts = []
        
        for run_dir in scenario_dir.iterdir():
            if not run_dir.is_dir() or not run_dir.name.startswith('run_'):
                continue
            
            results_file = run_dir / "results.json"
            if not results_file.exists():
                continue
            
            try:
                with open(results_file) as f:
                    data = json.load(f)
                
                # Count species
                molecules = data.get('molecules_detected', [])
                species_counts.append(len(molecules))
                
                # Count cycles
                cycles_file = run_dir / "autocatalytic_cycles.json"
                if cycles_file.exists():
                    with open(cycles_file) as cf:
                        cycles_data = json.load(cf)
                        cycle_counts.append(len(cycles_data))
                else:
                    cycle_counts.append(0)
            
            except Exception as e:
                print(f"Error reading {results_file}: {e}")
                continue
        
        if species_counts:
            scenarios[scenario] = {
                'species_counts': species_counts,
                'cycle_counts': cycle_counts
            }
    
    return scenarios

def calculate_kruskal_wallis(data_dict):
    """Calculate Kruskal-Wallis test for multiple groups"""
    groups = list(data_dict.values())
    if len(groups) < 2:
        return None, None
    
    # Flatten groups for scipy
    all_data = []
    group_labels = []
    for i, (name, values) in enumerate(data_dict.items()):
        all_data.extend(values)
        group_labels.extend([i] * len(values))
    
    if len(set(group_labels)) < 2:
        return None, None
    
    # Perform test
    h_stat, p_value = stats.kruskal(*groups)
    return h_stat, p_value

def calculate_mann_whitney(group1, group2):
    """Calculate Mann-Whitney U test between two groups"""
    if len(group1) == 0 or len(group2) == 0:
        return None, None
    
    u_stat, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')
    return u_stat, p_value

def main():
    scenarios = load_scenario_data()
    
    print("="*70)
    print("STATISTICAL TESTS FOR PHASE 2B RESULTS")
    print("="*70)
    
    # Extract data
    species_data = {k: v['species_counts'] for k, v in scenarios.items()}
    cycle_data = {k: v['cycle_counts'] for k, v in scenarios.items()}
    
    # Kruskal-Wallis for species diversity
    print("\n1. Species Diversity (Kruskal-Wallis H-test):")
    print(f"   Groups: {list(species_data.keys())}")
    for name, values in species_data.items():
        print(f"   {name}: n={len(values)}, mean={np.mean(values):.1f} ± {np.std(values):.1f}")
    
    h_stat, p_value = calculate_kruskal_wallis(species_data)
    if h_stat is not None:
        print(f"   H-statistic: {h_stat:.4f}")
        print(f"   p-value: {p_value:.6f}")
        print(f"   Significant: {'Yes' if p_value < 0.05 else 'No'} (p < 0.05)")
    else:
        print("   Cannot calculate (insufficient data)")
    
    # Kruskal-Wallis for cycle counts
    print("\n2. Autocatalytic Cycles (Kruskal-Wallis H-test):")
    print(f"   Groups: {list(cycle_data.keys())}")
    for name, values in cycle_data.items():
        print(f"   {name}: n={len(values)}, mean={np.mean(values):.1f} ± {np.std(values):.1f}")
    
    h_stat_cycles, p_value_cycles = calculate_kruskal_wallis(cycle_data)
    if h_stat_cycles is not None:
        print(f"   H-statistic: {h_stat_cycles:.4f}")
        print(f"   p-value: {p_value_cycles:.6f}")
        print(f"   Significant: {'Yes' if p_value_cycles < 0.05 else 'No'} (p < 0.05)")
    else:
        print("   Cannot calculate (insufficient data)")
    
    # Pairwise Mann-Whitney tests
    print("\n3. Pairwise Comparisons (Mann-Whitney U-test):")
    scenario_names = list(species_data.keys())
    for i in range(len(scenario_names)):
        for j in range(i+1, len(scenario_names)):
            name1, name2 = scenario_names[i], scenario_names[j]
            u_stat, p_val = calculate_mann_whitney(
                species_data[name1], 
                species_data[name2]
            )
            if u_stat is not None:
                print(f"   {name1} vs {name2}: U={u_stat:.1f}, p={p_val:.6f}")
    
    # Save results
    results = {
        'species_diversity': {
            'test': 'Kruskal-Wallis',
            'h_statistic': float(h_stat) if h_stat is not None else None,
            'p_value': float(p_value) if p_value is not None else None,
            'significant': bool(p_value < 0.05) if p_value is not None else None
        },
        'autocatalytic_cycles': {
            'test': 'Kruskal-Wallis',
            'h_statistic': float(h_stat_cycles) if h_stat_cycles is not None else None,
            'p_value': float(p_value_cycles) if p_value_cycles is not None else None,
            'significant': bool(p_value_cycles < 0.05) if p_value_cycles is not None else None
        },
        'scenario_data': {
            k: {
                'n_runs': len(v['species_counts']),
                'species_mean': float(np.mean(v['species_counts'])),
                'species_std': float(np.std(v['species_counts'])),
                'cycles_mean': float(np.mean(v['cycle_counts'])),
                'cycles_std': float(np.std(v['cycle_counts']))
            }
            for k, v in scenarios.items()
        }
    }
    
    output_file = Path('paper/results_data/statistical_tests.json')
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to {output_file}")
    print("="*70)

if __name__ == "__main__":
    main()

