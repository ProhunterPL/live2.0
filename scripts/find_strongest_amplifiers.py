#!/usr/bin/env python3
"""Find strongest amplifying molecules and formose-like cycles"""

import json
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple


def analyze_all_cycles() -> Dict:
    """Analyze all autocatalytic cycles to find strongest amplifiers"""
    results = {
        'strongest_amplifiers': [],
        'formose_runs': 0,
        'glycolaldehyde_max_amp': 0.0,
        'glycolaldehyde_time': 0
    }
    
    cycles_files = list(Path('results/phase2b_additional').rglob('autocatalytic_cycles.json'))
    print(f"Analyzing {len(cycles_files)} cycles files...")
    
    max_amp = 0.0
    max_amp_molecules = []
    max_amp_scenario = None
    
    formose_runs = set()
    glycolaldehyde_amps = []
    
    for cycles_file in cycles_files:
        scenario = cycles_file.parent.parent.name
        run_id = cycles_file.parent.name
        
        try:
            with open(cycles_file) as f:
                cycles = json.load(f)
            
            for cycle in cycles:
                amp = cycle.get('amplification', 0.0)
                molecules = cycle.get('molecules', [])
                cycle_type = cycle.get('type', '')
                
                # Track strongest amplifiers
                if amp > max_amp:
                    max_amp = amp
                    max_amp_molecules = molecules[:3]  # Top 3 molecules
                    max_amp_scenario = scenario
                
                # Check for formose-like cycles (formaldehyde or glycolaldehyde)
                mol_str = ' '.join(molecules).upper()
                if 'CH2O' in mol_str or 'FORMALDEHYDE' in mol_str or 'GLYCOLALDEHYDE' in mol_str:
                    formose_runs.add((scenario, run_id))
                    if 'GLYCOLALDEHYDE' in mol_str or 'C2H4O2' in mol_str:
                        glycolaldehyde_amps.append({
                            'amp': amp,
                            'step': cycle.get('first_step', 0),
                            'scenario': scenario,
                            'run': run_id
                        })
        
        except Exception as e:
            print(f"Error reading {cycles_file}: {e}")
            continue
    
    results['strongest_amplifiers'] = max_amp_molecules
    results['strongest_scenario'] = max_amp_scenario
    results['max_amplification'] = max_amp
    results['formose_runs'] = len(formose_runs)
    
    if glycolaldehyde_amps:
        max_glyc = max(glycolaldehyde_amps, key=lambda x: x['amp'])
        results['glycolaldehyde_max_amp'] = max_glyc['amp']
        results['glycolaldehyde_time'] = max_glyc['step']
    
    return results


if __name__ == "__main__":
    results = analyze_all_cycles()
    
    print("\n" + "="*70)
    print("STRONGEST AMPLIFIERS ANALYSIS")
    print("="*70)
    print(f"Strongest amplifiers: {', '.join(results['strongest_amplifiers'])}")
    print(f"Scenario: {results['strongest_scenario']}")
    print(f"Max amplification: {results['max_amplification']:.2f}-fold")
    print(f"\nFormose-like cycles found in {results['formose_runs']} runs")
    print(f"Max glycolaldehyde amplification: {results['glycolaldehyde_max_amp']:.2f}-fold")
    print(f"Time to max: {results['glycolaldehyde_time']:,} steps")
    print("="*70)
    
    # Save to JSON
    with open('paper/results_data/strongest_amplifiers.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n✅ Results saved to paper/results_data/strongest_amplifiers.json")

