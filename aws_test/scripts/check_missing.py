#!/usr/bin/env python3
"""Check what's missing"""
import os
from pathlib import Path

def check_run(run_path):
    """Check if run is complete"""
    has_results = os.path.exists(os.path.join(run_path, 'results.json'))
    has_molecules = os.path.exists(os.path.join(run_path, 'molecules.json'))
    has_snapshots = os.path.exists(os.path.join(run_path, 'snapshots'))
    return has_results and has_molecules and has_snapshots

scenarios = ['miller_urey_extended', 'hydrothermal_extended', 'formamide_extended']

print("📊 Status danych lokalnych:\n")

for scenario in scenarios:
    path = Path(f'results/phase2b_additional/{scenario}')
    if not path.exists():
        print(f"{scenario}: BRAK KATALOGU")
        continue
    
    runs = sorted([d for d in os.listdir(path) if os.path.isdir(path / d) and d.startswith('run_')], 
                  key=lambda x: int(x.split('_')[1]))
    
    complete = [r for r in runs if check_run(path / r)]
    incomplete = [r for r in runs if not check_run(path / r)]
    
    print(f"{scenario}:")
    print(f"  Total: {len(runs)} runs")
    print(f"  Complete: {len(complete)}")
    print(f"  Incomplete: {len(incomplete)}")
    if incomplete:
        print(f"  Missing: {', '.join(incomplete)}")
    print()

