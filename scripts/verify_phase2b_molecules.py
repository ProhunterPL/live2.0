#!/usr/bin/env python3
"""
Verify Phase 2B molecules.json files availability and content
"""

import json
from pathlib import Path
from collections import defaultdict

def verify_molecules_files():
    """Check all Phase 2B molecules.json files"""
    
    base_dir = Path("results/phase2b_additional")
    
    scenarios = {
        "miller_urey_extended": 18,
        "hydrothermal_extended": 17,
        "formamide_extended": 8
    }
    
    results = defaultdict(dict)
    total_runs = 0
    total_with_molecules = 0
    total_with_content = 0
    total_empty = 0
    
    print("=" * 70)
    print("Phase 2B molecules.json Verification")
    print("=" * 70)
    print()
    
    for scenario, expected_runs in scenarios.items():
        scenario_dir = base_dir / scenario
        
        if not scenario_dir.exists():
            print(f"❌ {scenario}: Directory not found")
            continue
        
        runs = sorted([d for d in scenario_dir.iterdir() if d.is_dir() and d.name.startswith("run_")])
        actual_runs = len(runs)
        
        results[scenario]["expected"] = expected_runs
        results[scenario]["actual_runs"] = actual_runs
        results[scenario]["with_molecules_file"] = 0
        results[scenario]["with_content"] = 0
        results[scenario]["empty"] = 0
        results[scenario]["missing"] = []
        
        for run_dir in runs:
            total_runs += 1
            molecules_file = run_dir / "molecules.json"
            results_file = run_dir / "results.json"
            
            if molecules_file.exists():
                total_with_molecules += 1
                results[scenario]["with_molecules_file"] += 1
                
                # Check content
                try:
                    with open(molecules_file, 'r') as f:
                        content = json.load(f)
                    
                    if isinstance(content, list) and len(content) > 0:
                        total_with_content += 1
                        results[scenario]["with_content"] += 1
                    elif isinstance(content, dict) and len(content) > 0:
                        total_with_content += 1
                        results[scenario]["with_content"] += 1
                    else:
                        total_empty += 1
                        results[scenario]["empty"] += 1
                except (json.JSONDecodeError, Exception) as e:
                    print(f"⚠️  Error reading {molecules_file}: {e}")
                    total_empty += 1
                    results[scenario]["empty"] += 1
            else:
                results[scenario]["missing"].append(run_dir.name)
        
        # Print scenario summary
        print(f"📊 {scenario}:")
        print(f"   Expected runs: {expected_runs}")
        print(f"   Actual runs: {actual_runs}")
        print(f"   With molecules.json: {results[scenario]['with_molecules_file']}/{actual_runs}")
        print(f"   With content: {results[scenario]['with_content']}")
        print(f"   Empty: {results[scenario]['empty']}")
        if results[scenario]["missing"]:
            print(f"   ⚠️  Missing molecules.json: {', '.join(results[scenario]['missing'])}")
        print()
    
    # Overall summary
    print("=" * 70)
    print("Overall Summary")
    print("=" * 70)
    print(f"Total runs: {total_runs}")
    print(f"Runs with molecules.json file: {total_with_molecules}/{total_runs} ({100*total_with_molecules/total_runs:.1f}%)")
    print(f"Runs with content in molecules.json: {total_with_content}/{total_runs} ({100*total_with_content/total_runs:.1f}%)")
    print(f"Runs with empty molecules.json: {total_empty}/{total_runs} ({100*total_empty/total_runs:.1f}%)")
    print()
    
    # Status
    if total_with_molecules == total_runs and total_with_content == total_runs:
        print("✅ STATUS: All runs have molecules.json files with content!")
    elif total_with_molecules == total_runs:
        print("⚠️  STATUS: All runs have molecules.json files, but some are empty")
        print("   → May need to extract molecules from snapshots")
    else:
        print("❌ STATUS: Some runs are missing molecules.json files")
    
    return results

if __name__ == "__main__":
    verify_molecules_files()

