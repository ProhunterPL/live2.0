#!/usr/bin/env python3
"""Extract novel molecules data from Phase 2B results"""

import json
from pathlib import Path
from typing import Dict, List
from collections import defaultdict
import statistics


def extract_novel_molecules_stats() -> Dict:
    """Extract statistics about novel molecules across all runs"""
    base_dir = Path('results/phase2b_additional')
    
    all_novel = []
    all_known = []
    novel_by_scenario = defaultdict(list)
    known_by_scenario = defaultdict(list)
    
    total_molecules = 0
    total_novel = 0
    
    for scenario_dir in base_dir.iterdir():
        if not scenario_dir.is_dir():
            continue
        
        scenario = scenario_dir.name
        
        for run_dir in scenario_dir.iterdir():
            if not run_dir.is_dir() or not run_dir.name.startswith('run_'):
                continue
            
            results_file = run_dir / "results.json"
            if not results_file.exists():
                continue
            
            try:
                with open(results_file) as f:
                    data = json.load(f)
                
                molecules = data.get('molecules_detected', [])
                novel_molecules = data.get('novel_molecules', [])
                
                total_molecules += len(molecules)
                total_novel += len(novel_molecules)
                
                # Extract detection steps
                for mol in molecules:
                    step = mol.get('first_seen', mol.get('step', 0))
                    if mol.get('id') in [n.get('id') for n in novel_molecules]:
                        all_novel.append(step)
                        novel_by_scenario[scenario].append(step)
                    else:
                        all_known.append(step)
                        known_by_scenario[scenario].append(step)
            
            except Exception as e:
                print(f"Error reading {results_file}: {e}")
                continue
    
    # Calculate statistics
    stats = {
        'total_novel': total_novel,
        'total_molecules': total_molecules,
        'novel_percent': (total_novel / total_molecules * 100) if total_molecules > 0 else 0,
        'novel_median_step': int(statistics.median(all_novel)) if all_novel else 0,
        'known_median_step': int(statistics.median(all_known)) if all_known else 0,
        'scenario_specific': {}
    }
    
    # Scenario-specific percentages (simplified - would need overlap analysis)
    for scenario in novel_by_scenario:
        scenario_novel = len(novel_by_scenario[scenario])
        scenario_total = scenario_novel + len(known_by_scenario[scenario])
        stats['scenario_specific'][scenario] = {
            'novel_count': scenario_novel,
            'novel_percent': (scenario_novel / scenario_total * 100) if scenario_total > 0 else 0
        }
    
    return stats


if __name__ == "__main__":
    stats = extract_novel_molecules_stats()
    
    print("="*70)
    print("NOVEL MOLECULES STATISTICS")
    print("="*70)
    print(f"Total novel molecules: {stats['total_novel']:,}")
    print(f"Total molecules: {stats['total_molecules']:,}")
    print(f"Novel percentage: {stats['novel_percent']:.1f}%")
    print(f"Novel median detection step: {stats['novel_median_step']:,}")
    print(f"Known median detection step: {stats['known_median_step']:,}")
    print("\nScenario-specific:")
    for scenario, data in stats['scenario_specific'].items():
        print(f"  {scenario}: {data['novel_percent']:.1f}% novel")
    print("="*70)
    
    # Save to JSON
    with open('paper/results_data/novel_molecules_stats.json', 'w') as f:
        json.dump(stats, f, indent=2)
    
    print("\n✅ Stats saved to paper/results_data/novel_molecules_stats.json")

