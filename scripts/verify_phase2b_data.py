#!/usr/bin/env python3
"""
Verify Phase 2B data completeness for Paper 2 analysis.
"""

import json
import os
from pathlib import Path
from collections import defaultdict

def verify_phase2b_data():
    """Verify Phase 2B data completeness"""
    
    base_dir = Path("results/phase2b_additional")
    
    scenarios = {
        "miller_urey_extended": 18,
        "hydrothermal_extended": 17,
        "formamide_extended": 8
    }
    
    results = {
        "total_runs_expected": 43,
        "total_runs_found": 0,
        "scenarios": {},
        "autocatalytic_cycles": {
            "total": 0,
            "by_scenario": defaultdict(int),
            "runs_with_cycles": 0,
            "runs_without_cycles": 0
        },
        "missing_files": [],
        "incomplete_runs": []
    }
    
    print("=" * 80)
    print("PHASE 2B DATA VERIFICATION")
    print("=" * 80)
    print()
    
    # Check each scenario
    for scenario, expected_count in scenarios.items():
        scenario_dir = base_dir / scenario
        scenario_results = {
            "expected": expected_count,
            "found": 0,
            "runs": [],
            "cycles_count": 0,
            "runs_with_cycles": 0
        }
        
        if not scenario_dir.exists():
            print(f"❌ Scenario directory not found: {scenario_dir}")
            results["missing_files"].append(str(scenario_dir))
            continue
        
        # Find all run directories
        run_dirs = sorted([d for d in scenario_dir.iterdir() 
                          if d.is_dir() and d.name.startswith("run_")])
        
        scenario_results["found"] = len(run_dirs)
        results["total_runs_found"] += len(run_dirs)
        
        print(f"📊 {scenario}:")
        print(f"   Expected: {expected_count} runs")
        print(f"   Found: {len(run_dirs)} runs")
        
        if len(run_dirs) != expected_count:
            print(f"   ⚠️  Mismatch: expected {expected_count}, found {len(run_dirs)}")
        
        # Check each run
        for run_dir in run_dirs:
            run_info = {
                "run_id": run_dir.name,
                "has_results_json": False,
                "has_molecules_json": False,
                "has_cycles_json": False,
                "has_snapshots": False,
                "cycles_count": 0
            }
            
            # Check results.json
            results_file = run_dir / "results.json"
            if results_file.exists():
                run_info["has_results_json"] = True
                try:
                    with open(results_file) as f:
                        results_data = json.load(f)
                except Exception as e:
                    print(f"   ⚠️  Error reading {results_file}: {e}")
            
            # Check molecules.json
            molecules_file = run_dir / "molecules.json"
            if molecules_file.exists():
                run_info["has_molecules_json"] = True
            
            # Check autocatalytic_cycles.json
            cycles_file = run_dir / "autocatalytic_cycles.json"
            if cycles_file.exists():
                run_info["has_cycles_json"] = True
                try:
                    with open(cycles_file) as f:
                        cycles_data = json.load(f)
                        if isinstance(cycles_data, list):
                            run_info["cycles_count"] = len(cycles_data)
                        elif isinstance(cycles_data, dict):
                            # Could be structured differently
                            run_info["cycles_count"] = cycles_data.get("total_cycles", 
                                                                      len(cycles_data.get("cycles", [])))
                        scenario_results["cycles_count"] += run_info["cycles_count"]
                        results["autocatalytic_cycles"]["total"] += run_info["cycles_count"]
                        results["autocatalytic_cycles"]["by_scenario"][scenario] += run_info["cycles_count"]
                        scenario_results["runs_with_cycles"] += 1
                        results["autocatalytic_cycles"]["runs_with_cycles"] += 1
                except Exception as e:
                    print(f"   ⚠️  Error reading {cycles_file}: {e}")
            else:
                results["autocatalytic_cycles"]["runs_without_cycles"] += 1
            
            # Check snapshots directory
            snapshots_dir = run_dir / "snapshots"
            if snapshots_dir.exists() and any(snapshots_dir.iterdir()):
                run_info["has_snapshots"] = True
            
            # Check if run is complete
            if not (run_info["has_results_json"] and run_info["has_molecules_json"]):
                results["incomplete_runs"].append(f"{scenario}/{run_dir.name}")
            
            scenario_results["runs"].append(run_info)
        
        results["scenarios"][scenario] = scenario_results
        
        print(f"   Cycles detected: {scenario_results['cycles_count']:,}")
        print(f"   Runs with cycles: {scenario_results['runs_with_cycles']}/{scenario_results['found']}")
        print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total runs expected: {results['total_runs_expected']}")
    print(f"Total runs found: {results['total_runs_found']}")
    print(f"Total autocatalytic cycles: {results['autocatalytic_cycles']['total']:,}")
    print(f"Runs with cycles: {results['autocatalytic_cycles']['runs_with_cycles']}")
    print(f"Runs without cycles: {results['autocatalytic_cycles']['runs_without_cycles']}")
    print()
    
    if results["incomplete_runs"]:
        print(f"⚠️  Incomplete runs ({len(results['incomplete_runs'])}):")
        for run in results["incomplete_runs"][:10]:  # Show first 10
            print(f"   - {run}")
        if len(results["incomplete_runs"]) > 10:
            print(f"   ... and {len(results['incomplete_runs']) - 10} more")
        print()
    
    if results["missing_files"]:
        print(f"❌ Missing files/directories ({len(results['missing_files'])}):")
        for f in results["missing_files"]:
            print(f"   - {f}")
        print()
    
    # Check if we have 769,315 cycles
    expected_cycles = 769315
    found_cycles = results["autocatalytic_cycles"]["total"]
    
    print("=" * 80)
    print("CYCLE COUNT VERIFICATION")
    print("=" * 80)
    print(f"Expected cycles (from Paper 1): {expected_cycles:,}")
    print(f"Found cycles (from autocatalytic_cycles.json files): {found_cycles:,}")
    
    if found_cycles == 0:
        print("⚠️  WARNING: No cycles found in autocatalytic_cycles.json files")
        print("   This might mean:")
        print("   - Cycles are stored in a different format")
        print("   - Cycles need to be extracted from snapshots")
        print("   - Analysis needs to be re-run")
    elif abs(found_cycles - expected_cycles) > 1000:
        print(f"⚠️  Mismatch: difference of {abs(found_cycles - expected_cycles):,} cycles")
    else:
        print("✅ Cycle count matches (within tolerance)")
    
    print()
    
    # Save results
    output_file = Path("docs/plans/PHASE2B_DATA_VERIFICATION.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"✅ Verification results saved to: {output_file}")
    
    return results

if __name__ == "__main__":
    verify_phase2b_data()

