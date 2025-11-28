#!/usr/bin/env python3
"""Check completeness of Phase 2B data"""
import os
from pathlib import Path

def check_scenario(scenario_path, scenario_name):
    """Check completeness of a scenario"""
    if not os.path.exists(scenario_path):
        return {'name': scenario_name, 'runs': 0, 'complete': 0, 'incomplete': []}
    
    runs = []
    for d in os.listdir(scenario_path):
        run_path = os.path.join(scenario_path, d)
        if os.path.isdir(run_path) and d.startswith('run_'):
            has_results = os.path.exists(os.path.join(run_path, 'results.json'))
            has_molecules = os.path.exists(os.path.join(run_path, 'molecules.json'))
            has_snapshots = os.path.exists(os.path.join(run_path, 'snapshots'))
            has_summary = os.path.exists(os.path.join(run_path, 'summary.txt'))
            
            complete = has_results and has_molecules and has_snapshots
            runs.append({
                'run': d,
                'complete': complete,
                'results': has_results,
                'molecules': has_molecules,
                'snapshots': has_snapshots,
                'summary': has_summary
            })
    
    complete_runs = [r for r in runs if r['complete']]
    incomplete_runs = [r for r in runs if not r['complete']]
    
    return {
        'name': scenario_name,
        'runs': len(runs),
        'complete': len(complete_runs),
        'incomplete': incomplete_runs
    }

def main():
    base_path = Path('results/phase2b_additional')
    
    scenarios = [
        ('miller_urey_extended', 'Miller-Urey'),
        ('hydrothermal_extended', 'Hydrothermal'),
        ('formamide_extended', 'Formamide')
    ]
    
    print("📊 Phase 2B Data Completeness Check\n")
    print("=" * 60)
    
    total_runs = 0
    total_complete = 0
    
    for scenario_dir, scenario_name in scenarios:
        scenario_path = base_path / scenario_dir
        result = check_scenario(str(scenario_path), scenario_name)
        
        total_runs += result['runs']
        total_complete += result['complete']
        
        print(f"\n{scenario_name}:")
        print(f"  Total runs: {result['runs']}")
        print(f"  Complete: {result['complete']}")
        print(f"  Incomplete: {len(result['incomplete'])}")
        
        if result['incomplete']:
            print(f"  Incomplete runs:")
            for r in result['incomplete'][:5]:
                missing = []
                if not r['results']: missing.append('results.json')
                if not r['molecules']: missing.append('molecules.json')
                if not r['snapshots']: missing.append('snapshots')
                print(f"    - {r['run']}: missing {', '.join(missing)}")
            if len(result['incomplete']) > 5:
                print(f"    ... and {len(result['incomplete']) - 5} more")
    
    print("\n" + "=" * 60)
    print(f"\n📈 Summary:")
    print(f"  Total runs: {total_runs}")
    print(f"  Complete runs: {total_complete}")
    print(f"  Incomplete runs: {total_runs - total_complete}")
    print(f"  Completion rate: {total_complete/total_runs*100:.1f}%" if total_runs > 0 else "  Completion rate: N/A")
    
    if total_complete == total_runs and total_runs > 0:
        print("\n✅ All runs are complete!")
    elif total_complete >= total_runs * 0.9:
        print("\n⚠️ Most runs complete, but some are missing data")
    else:
        print("\n❌ Many runs are incomplete")

if __name__ == "__main__":
    main()

