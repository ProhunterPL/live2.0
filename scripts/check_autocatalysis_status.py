#!/usr/bin/env python3
"""Check status of autocatalysis analysis for Phase 2B runs"""

import json
from pathlib import Path
from collections import defaultdict

base = Path('results/phase2b_additional')
scenarios = ['miller_urey_extended', 'hydrothermal_extended']

results = defaultdict(lambda: {
    'with_history': [],
    'without_history': [],
    'with_cycles_file': [],
    'cycles_detected': [],
    'cycles_empty': [],
    'timeout_issues': []
})

for scenario in scenarios:
    scenario_dir = base / scenario
    if not scenario_dir.exists():
        continue
    
    for run_dir in sorted(scenario_dir.iterdir()):
        if not run_dir.is_dir() or not run_dir.name.startswith('run_'):
            continue
        
        network_file = run_dir / 'reaction_network.json'
        cycles_file = run_dir / 'autocatalytic_cycles.json'
        
        # Check reaction_network.json
        if network_file.exists():
            try:
                with open(network_file) as f:
                    network = json.load(f)
                
                has_history = 'abundance_history' in network
                if has_history:
                    results[scenario]['with_history'].append(run_dir.name)
                else:
                    results[scenario]['without_history'].append(run_dir.name)
            except Exception as e:
                print(f"Error reading {network_file}: {e}")
        
        # Check autocatalytic_cycles.json
        if cycles_file.exists():
            results[scenario]['with_cycles_file'].append(run_dir.name)
            try:
                with open(cycles_file) as f:
                    cycles = json.load(f)
                
                if isinstance(cycles, list):
                    if len(cycles) > 0:
                        results[scenario]['cycles_detected'].append(run_dir.name)
                    else:
                        results[scenario]['cycles_empty'].append(run_dir.name)
                elif isinstance(cycles, dict):
                    if cycles.get('cycles'):
                        results[scenario]['cycles_detected'].append(run_dir.name)
                    else:
                        results[scenario]['cycles_empty'].append(run_dir.name)
                    # Check for timeout
                    if 'timeout' in str(cycles.get('metadata', {})).lower():
                        results[scenario]['timeout_issues'].append(run_dir.name)
            except Exception as e:
                print(f"Error reading {cycles_file}: {e}")

# Print summary
print("="*70)
print("AUTOCATALYSIS ANALYSIS STATUS")
print("="*70)

for scenario in scenarios:
    print(f"\n{scenario.upper()}:")
    print(f"  Total runs: {len(results[scenario]['with_history']) + len(results[scenario]['without_history'])}")
    print(f"  With abundance_history: {len(results[scenario]['with_history'])}")
    print(f"  Without abundance_history: {len(results[scenario]['without_history'])}")
    if results[scenario]['without_history']:
        print(f"    Missing: {', '.join(results[scenario]['without_history'])}")
    print(f"  With cycles file: {len(results[scenario]['with_cycles_file'])}")
    print(f"  Cycles detected: {len(results[scenario]['cycles_detected'])}")
    print(f"  Cycles empty: {len(results[scenario]['cycles_empty'])}")
    if results[scenario]['timeout_issues']:
        print(f"  ⚠️  Timeout issues: {', '.join(results[scenario]['timeout_issues'])}")

print("\n" + "="*70)

