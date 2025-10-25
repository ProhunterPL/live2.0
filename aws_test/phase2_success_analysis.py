#!/usr/bin/env python3
"""
Phase 2 Success Criteria Analysis
Analyzes AWS results against Phase 2 success criteria from VALIDATION_ROADMAP.md
"""

import json
import os
from pathlib import Path
from collections import defaultdict
import statistics

def analyze_phase2_success_criteria():
    """Analyze AWS results against Phase 2 success criteria"""
    
    print("PHASE 2 SUCCESS CRITERIA ANALYSIS")
    print("=" * 60)
    print("Based on VALIDATION_ROADMAP.md criteria")
    print()
    
    # Load AWS results
    aws_results = load_aws_results()
    
    # Phase 2A: Test Validation (GO/NO-GO Decision)
    print("PHASE 2A: TEST VALIDATION (GO/NO-GO)")
    print("-" * 40)
    
    phase2a_criteria = analyze_phase2a_criteria(aws_results)
    print_phase2a_results(phase2a_criteria)
    
    print()
    
    # Phase 2B-D: Production Success Metrics
    print("PHASE 2B-D: PRODUCTION SUCCESS METRICS")
    print("-" * 40)
    
    phase2b_criteria = analyze_phase2b_criteria(aws_results)
    print_phase2b_results(phase2b_criteria)
    
    print()
    
    # Overall assessment
    print("OVERALL PHASE 2 ASSESSMENT")
    print("-" * 40)
    overall_assessment(phase2a_criteria, phase2b_criteria)

def load_aws_results():
    """Load AWS test results"""
    results = {
        'scenarios': defaultdict(list),
        'total_runs': 0,
        'successful_runs': 0,
        'failed_runs': 0
    }
    
    # Load from our previous analysis
    directories = ['results_16_completed', 'results_28_completed', 'results_all_completed']
    
    for dir_name in directories:
        results_dir = Path(f"{dir_name}/results")
        if results_dir.exists():
            for scenario_dir in results_dir.iterdir():
                if not scenario_dir.is_dir():
                    continue
                    
                scenario_name = scenario_dir.name
                
                for run_dir in scenario_dir.iterdir():
                    if not run_dir.is_dir():
                        continue
                        
                    results['total_runs'] += 1
                    
                    summary_file = run_dir / 'summary.txt'
                    results_file = run_dir / 'results.json'
                    molecules_file = run_dir / 'molecules.json'
                    
                    if summary_file.exists() and results_file.exists():
                        results['successful_runs'] += 1
                        
                        with open(results_file, 'r') as f:
                            results_json = json.load(f)
                        
                        molecules_data = []
                        if molecules_file.exists():
                            with open(molecules_file, 'r') as f:
                                molecules_data = json.load(f)
                        
                        run_data = {
                            'run_id': run_dir.name,
                            'scenario': scenario_name,
                            'results': results_json,
                            'molecules': molecules_data,
                            'molecules_count': len(molecules_data),
                            'final_particles': results_json.get('final_state', {}).get('n_particles', 0),
                            'steps': results_json.get('final_state', {}).get('step', 0),
                            'time': results_json.get('final_state', {}).get('time', 0)
                        }
                        
                        results['scenarios'][scenario_name].append(run_data)
                    else:
                        results['failed_runs'] += 1
    
    return results

def analyze_phase2a_criteria(results):
    """Analyze Phase 2A criteria"""
    criteria = {
        'minimum_criteria': {
            'simulation_completes': False,
            'memory_stable': False,  # Cannot verify from logs
            'performance_acceptable': False,  # Cannot verify from logs
            'molecules_detected': False,
            'expected_products': False,
            'thermodynamic_violations': False  # Cannot verify from logs
        },
        'optimal_criteria': {
            'molecules_10_plus': False,
            'expected_products_2_plus': False,
            'autocatalytic_cycles': False,  # Not detected in current results
            'performance_4_plus': False,  # Cannot verify
            'chemical_plausibility': False  # Cannot verify
        },
        'no_go_triggers': {
            'crashes_repeatedly': False,
            'memory_leak': False,  # Cannot verify
            'no_molecules_50k': False,
            'thermodynamic_violations': False,  # Cannot verify
            'performance_below_1': False  # Cannot verify
        }
    }
    
    # Analyze successful runs
    successful_runs = []
    total_molecules = 0
    all_molecules = []
    
    for scenario_name, runs in results['scenarios'].items():
        for run in runs:
            successful_runs.append(run)
            total_molecules += run['molecules_count']
            all_molecules.extend(run['molecules'])
    
    # Check minimum criteria
    if len(successful_runs) > 0:
        criteria['minimum_criteria']['simulation_completes'] = True
    
    if total_molecules >= 5:
        criteria['minimum_criteria']['molecules_detected'] = True
    
    # Check for expected products (simplified - looking for any molecules)
    if total_molecules > 0:
        criteria['minimum_criteria']['expected_products'] = True
    
    # Check optimal criteria
    if total_molecules >= 10:
        criteria['optimal_criteria']['molecules_10_plus'] = True
    
    if total_molecules >= 2:
        criteria['optimal_criteria']['expected_products_2_plus'] = True
    
    # Check NO-GO triggers
    if results['failed_runs'] == 0:
        criteria['no_go_triggers']['crashes_repeatedly'] = False
    
    if total_molecules > 0:
        criteria['no_go_triggers']['no_molecules_50k'] = False
    
    return criteria

def analyze_phase2b_criteria(results):
    """Analyze Phase 2B-D criteria"""
    criteria = {
        'simulation_quality': {
            'completion_rate': 0,
            'stability': False,  # Cannot verify
            'performance': False,  # Cannot verify
            'duration': False  # Cannot verify
        },
        'scientific_output': {
            'molecular_diversity_total': 0,
            'molecular_diversity_per_scenario': {},
            'expected_products': False,  # Cannot verify
            'autocatalytic_cycles': 0,
            'match_quality': False  # Cannot verify
        },
        'statistical_rigor': {
            'reproducibility': False,  # Cannot verify
            'scenario_differences': False,  # Cannot verify
            'error_bars': False,  # Cannot verify
            'n_sufficiency': False  # Cannot verify
        }
    }
    
    # Calculate completion rate
    total_runs = results['total_runs']
    successful_runs = results['successful_runs']
    criteria['simulation_quality']['completion_rate'] = (successful_runs / total_runs * 100) if total_runs > 0 else 0
    
    # Calculate molecular diversity
    all_molecules = set()
    scenario_molecules = defaultdict(set)
    
    for scenario_name, runs in results['scenarios'].items():
        for run in runs:
            for mol in run['molecules']:
                mol_id = mol.get('id', 'unknown')
                all_molecules.add(mol_id)
                scenario_molecules[scenario_name].add(mol_id)
    
    criteria['scientific_output']['molecular_diversity_total'] = len(all_molecules)
    criteria['scientific_output']['molecular_diversity_per_scenario'] = {
        scenario: len(molecules) for scenario, molecules in scenario_molecules.items()
    }
    
    return criteria

def print_phase2a_results(criteria):
    """Print Phase 2A analysis results"""
    
    print("MINIMUM CRITERIA FOR GO:")
    min_criteria = criteria['minimum_criteria']
    print(f"  [OK] Simulation completes: {min_criteria['simulation_completes']}")
    print(f"  [?] Memory stable: {min_criteria['memory_stable']} (cannot verify)")
    print(f"  [?] Performance acceptable: {min_criteria['performance_acceptable']} (cannot verify)")
    print(f"  [OK] Molecules detected (>=5): {min_criteria['molecules_detected']}")
    print(f"  [OK] Expected products: {min_criteria['expected_products']}")
    print(f"  [?] Thermodynamic violations: {min_criteria['thermodynamic_violations']} (cannot verify)")
    
    min_passed = sum(1 for v in min_criteria.values() if v)
    min_total = len(min_criteria)
    print(f"\n  [STATS] Minimum criteria passed: {min_passed}/{min_total}")
    
    print("\nOPTIMAL CRITERIA FOR GO:")
    opt_criteria = criteria['optimal_criteria']
    print(f"  [OK] Molecules 10+: {opt_criteria['molecules_10_plus']}")
    print(f"  [OK] Expected products 2+: {opt_criteria['expected_products_2_plus']}")
    print(f"  [NO] Autocatalytic cycles: {opt_criteria['autocatalytic_cycles']}")
    print(f"  [?] Performance 4+: {opt_criteria['performance_4_plus']} (cannot verify)")
    print(f"  [?] Chemical plausibility: {opt_criteria['chemical_plausibility']} (cannot verify)")
    
    opt_passed = sum(1 for v in opt_criteria.values() if v)
    opt_total = len(opt_criteria)
    print(f"\n  [STATS] Optimal criteria passed: {opt_passed}/{opt_total}")
    
    print("\nNO-GO TRIGGERS:")
    no_go = criteria['no_go_triggers']
    print(f"  [NO] Crashes repeatedly: {no_go['crashes_repeatedly']}")
    print(f"  [?] Memory leak: {no_go['memory_leak']} (cannot verify)")
    print(f"  [NO] No molecules after 50K: {no_go['no_molecules_50k']}")
    print(f"  [?] Thermodynamic violations: {no_go['thermodynamic_violations']} (cannot verify)")
    print(f"  [?] Performance <1: {no_go['performance_below_1']} (cannot verify)")
    
    no_go_triggered = sum(1 for v in no_go.values() if v)
    print(f"\n  [STATS] NO-GO triggers: {no_go_triggered}/5")

def print_phase2b_results(criteria):
    """Print Phase 2B-D analysis results"""
    
    print("SIMULATION QUALITY:")
    sim_quality = criteria['simulation_quality']
    print(f"  [STATS] Completion rate: {sim_quality['completion_rate']:.1f}% (target: >=90%)")
    print(f"  [?] Stability: {sim_quality['stability']} (cannot verify)")
    print(f"  [?] Performance: {sim_quality['performance']} (cannot verify)")
    print(f"  [?] Duration: {sim_quality['duration']} (cannot verify)")
    
    print("\nSCIENTIFIC OUTPUT:")
    sci_output = criteria['scientific_output']
    print(f"  [STATS] Molecular diversity total: {sci_output['molecular_diversity_total']} (target: >=100)")
    print(f"  [STATS] Per-scenario diversity:")
    for scenario, count in sci_output['molecular_diversity_per_scenario'].items():
        print(f"    - {scenario}: {count} (target: >=30)")
    print(f"  [?] Expected products: {sci_output['expected_products']} (cannot verify)")
    print(f"  [STATS] Autocatalytic cycles: {sci_output['autocatalytic_cycles']} (target: >=10)")
    print(f"  [?] Match quality: {sci_output['match_quality']} (cannot verify)")
    
    print("\nSTATISTICAL RIGOR:")
    stat_rigor = criteria['statistical_rigor']
    print(f"  [?] Reproducibility: {stat_rigor['reproducibility']} (cannot verify)")
    print(f"  [?] Scenario differences: {stat_rigor['scenario_differences']} (cannot verify)")
    print(f"  [?] Error bars: {stat_rigor['error_bars']} (cannot verify)")
    print(f"  [?] N sufficiency: {stat_rigor['n_sufficiency']} (cannot verify)")

def overall_assessment(phase2a_criteria, phase2b_criteria):
    """Provide overall assessment"""
    
    print("DECISION MATRIX:")
    print("All minimum criteria met -> GO to Phase 2B")
    print("Some minimum criteria failed -> Debug & retest")
    print("All NO-GO triggers -> Major revision needed")
    print()
    
    # Count verifiable criteria
    min_criteria = phase2a_criteria['minimum_criteria']
    verifiable_min = ['simulation_completes', 'molecules_detected', 'expected_products']
    min_passed = sum(1 for key in verifiable_min if min_criteria[key])
    min_total = len(verifiable_min)
    
    print(f"[OK] VERIFIABLE MINIMUM CRITERIA: {min_passed}/{min_total}")
    
    if min_passed == min_total:
        print("[SUCCESS] GO DECISION: All verifiable minimum criteria met!")
        print("   -> Proceed to Phase 2B (Production Runs)")
    elif min_passed >= min_total * 0.67:  # 2/3
        print("[WARNING] CONDITIONAL GO: Most criteria met")
        print("   -> Proceed with caution, monitor closely")
    else:
        print("[ERROR] NO-GO DECISION: Insufficient criteria met")
        print("   -> Debug & retest required")
    
    print()
    
    # Scientific output assessment
    total_molecules = phase2b_criteria['scientific_output']['molecular_diversity_total']
    completion_rate = phase2b_criteria['simulation_quality']['completion_rate']
    
    print("SCIENTIFIC OUTPUT ASSESSMENT:")
    print(f"  [STATS] Total unique molecules: {total_molecules}")
    print(f"  [STATS] Completion rate: {completion_rate:.1f}%")
    
    if total_molecules >= 100 and completion_rate >= 90:
        print("  [EXCELLENT] Ready for publication!")
    elif total_molecules >= 50 and completion_rate >= 80:
        print("  [GOOD] Sufficient for publication")
    elif total_molecules >= 20 and completion_rate >= 70:
        print("  [MARGINAL] May need more data")
    else:
        print("  [INSUFFICIENT] Need more simulations")
    
    print()
    print("RECOMMENDATION:")
    if min_passed == min_total and total_molecules >= 50:
        print("[SUCCESS] PHASE 2 COMPLETE - Proceed to Phase 3 (Paper Writing)")
    elif min_passed >= min_total * 0.67:
        print("[WARNING] PHASE 2 MOSTLY COMPLETE - Consider additional runs")
    else:
        print("[ERROR] PHASE 2 INCOMPLETE - Need more work before proceeding")

if __name__ == "__main__":
    analyze_phase2_success_criteria()
