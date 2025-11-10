#!/usr/bin/env python3
"""
Fix Phase 2B molecules by extracting from snapshots
===================================================

This script applies the solution that worked locally:
- Extract molecules from snapshots for all runs
- Update results.json with extracted molecules
- This fixes the empty catalog problem

Usage:
    # Fix all runs in a scenario
    python scripts/fix_phase2b_molecules.py --scenario miller_urey_extended
    
    # Fix specific run
    python scripts/fix_phase2b_molecules.py --run results/phase2b_additional/miller_urey_extended/run_1
    
    # Fix all runs in all scenarios
    python scripts/fix_phase2b_molecules.py --all
"""

import sys
import json
import argparse
from pathlib import Path
from typing import List, Optional

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def extract_molecules_from_snapshots(run_dir: Path) -> Optional[List[dict]]:
    """Extract molecules from snapshots using post_detect_batch logic"""
    try:
        import numpy as np
        from collections import defaultdict, Counter
        from backend.sim.core.graphs import MolecularGraph
        from backend.sim.core.catalog import SubstanceCatalog
        
        snapshot_dir = run_dir / "snapshots"
        if not snapshot_dir.exists():
            return None
        
        snapshot_files = sorted(snapshot_dir.glob("step_*.json"))
        if not snapshot_files:
            return None
        
        # Process all snapshots
        all_molecules = []
        catalog = SubstanceCatalog()
        
        for snapshot_file in snapshot_files:
            with open(snapshot_file, 'r') as f:
                data = json.load(f)
            
            step = data.get('step', 0)
            positions = np.array(data.get('positions', []))
            attributes = np.array(data.get('attributes', []))
            bonds = data.get('bonds', [])
            clusters = data.get('clusters', [])
            
            # If no clusters, infer from bonds
            if not clusters and bonds:
                cluster_graph = defaultdict(set)
                for bond in bonds:
                    if len(bond) >= 2:
                        i, j = bond[0], bond[1]
                        cluster_graph[i].add(j)
                        cluster_graph[j].add(i)
                
                visited = set()
                clusters = []
                for node in cluster_graph:
                    if node not in visited:
                        cluster = []
                        stack = [node]
                        while stack:
                            n = stack.pop()
                            if n not in visited:
                                visited.add(n)
                                cluster.append(n)
                                stack.extend(cluster_graph[n])
                        if len(cluster) >= 2:
                            clusters.append(cluster)
            
            # Process each cluster
            for cluster in clusters:
                if len(cluster) >= 2:  # Minimum 2 atoms for a molecule
                    # Get cluster bonds
                    cluster_bonds = []
                    bond_dict = {}
                    for b in bonds:
                        if len(b) >= 2:
                            i, j = b[0], b[1]
                            if i in cluster and j in cluster:
                                bond_dict[(min(i,j), max(i,j))] = b[2] if len(b) > 2 else 1.0
                    
                    for (i, j), strength in bond_dict.items():
                        cluster_bonds.append((i, j))
                    
                    # Get attributes for cluster particles
                    particle_attributes = {}
                    for idx in cluster:
                        if idx < len(attributes):
                            attr = attributes[idx]
                            if hasattr(attr, 'tolist'):
                                particle_attributes[idx] = attr.tolist()
                            else:
                                particle_attributes[idx] = list(attr) if isinstance(attr, (list, tuple)) else [attr]
                    
                    # Create graph and add to catalog
                    try:
                        graph = MolecularGraph(cluster, cluster_bonds, particle_attributes)
                        is_novel, substance_id = catalog.add_substance(graph, 0.0, {})
                        
                        # Get formula from graph
                        formula = graph.get_canonical_form()
                        
                        molecule_info = {
                            'id': substance_id,
                            'formula': formula,
                            'step': step,
                            'cluster': cluster,
                            'size': len(cluster),
                            'bonds': cluster_bonds,
                            'first_seen': step,
                            'last_seen': step,
                            'count': 1
                        }
                        all_molecules.append(molecule_info)
                    except Exception as e:
                        continue
        
        # Aggregate molecules by formula
        formula_counts = Counter(m['formula'] for m in all_molecules)
        
        molecules = []
        seen_formulas = set()
        for mol in all_molecules:
            formula = mol['formula']
            if formula not in seen_formulas:
                mol['count'] = formula_counts[formula]
                molecules.append(mol)
                seen_formulas.add(formula)
        
        return molecules if molecules else None
        
    except Exception as e:
        print(f"   [ERROR] Extraction failed: {e}")
        return None

def fix_single_run(run_dir: Path, verbose: bool = True) -> bool:
    """Extract molecules from snapshots and update results.json for a single run"""
    
    if not run_dir.exists():
        if verbose:
            print(f"[ERROR] Directory not found: {run_dir}")
        return False
    
    if verbose:
        print(f"\nProcessing: {run_dir.name}")
        print("  [1/3] Extracting molecules from snapshots...")
    
    # Extract molecules
    molecules = extract_molecules_from_snapshots(run_dir)
    
    if not molecules:
        if verbose:
            print("  [WARNING] No molecules extracted - snapshots may not contain bond/cluster data")
        return False
    
    if verbose:
        print(f"  [OK] Extracted {len(molecules)} unique molecules")
    
    # Update results.json
    if verbose:
        print("  [2/3] Updating results.json...")
    results_file = run_dir / "results.json"
    
    if not results_file.exists():
        if verbose:
            print(f"  [ERROR] results.json not found: {results_file}")
        return False
    
    try:
        # Load current results
        with open(results_file, 'r') as f:
            results_data = json.load(f)
        
        # Update with extracted molecules
        results_data['molecules_detected'] = molecules
        results_data['novel_molecules'] = molecules  # For now, treat all as novel
        
        # Add extraction metadata
        if 'metadata' not in results_data:
            results_data['metadata'] = {}
        results_data['metadata']['molecules_extracted_from'] = 'snapshots'
        
        # Save updated results
        with open(results_file, 'w') as f:
            json.dump(results_data, f, indent=2)
        
        if verbose:
            print(f"  [OK] Updated results.json with {len(molecules)} molecules")
        
    except Exception as e:
        if verbose:
            print(f"  [ERROR] Failed to update results.json: {e}")
        return False
    
    # Update molecules.json
    if verbose:
        print("  [3/3] Updating molecules.json...")
    molecules_file = run_dir / "molecules.json"
    
    try:
        with open(molecules_file, 'w') as f:
            json.dump(molecules, f, indent=2)
        
        if verbose:
            print(f"  [OK] Updated molecules.json")
        
    except Exception as e:
        if verbose:
            print(f"  [WARNING] Failed to update molecules.json: {e}")
    
    return True

def fix_all_runs(base_dir: Path = Path("results/phase2b_additional"), 
                  scenario: Optional[str] = None,
                  verbose: bool = True) -> dict:
    """Fix molecules for all runs in specified scenarios"""
    
    scenarios = []
    if scenario:
        scenarios = [scenario]
    else:
        # Find all scenarios
        for item in base_dir.iterdir():
            if item.is_dir() and not item.name.startswith('.'):
                scenarios.append(item.name)
    
    results = {
        'total_runs': 0,
        'successful': 0,
        'failed': 0,
        'runs': []
    }
    
    for scenario_name in scenarios:
        scenario_dir = base_dir / scenario_name
        if not scenario_dir.exists():
            continue
        
        if verbose:
            print(f"\n{'='*70}")
            print(f"SCENARIO: {scenario_name}")
            print(f"{'='*70}")
        
        # Find all run directories
        run_dirs = sorted([d for d in scenario_dir.iterdir() 
                          if d.is_dir() and d.name.startswith('run_')])
        
        for run_dir in run_dirs:
            results['total_runs'] += 1
            success = fix_single_run(run_dir, verbose=verbose)
            
            if success:
                results['successful'] += 1
                results['runs'].append({
                    'scenario': scenario_name,
                    'run': run_dir.name,
                    'status': 'success'
                })
            else:
                results['failed'] += 1
                results['runs'].append({
                    'scenario': scenario_name,
                    'run': run_dir.name,
                    'status': 'failed'
                })
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description="Fix Phase 2B molecules by extracting from snapshots"
    )
    parser.add_argument(
        '--run',
        type=str,
        help='Fix specific run directory (e.g., results/phase2b_additional/miller_urey_extended/run_1)'
    )
    parser.add_argument(
        '--scenario',
        type=str,
        help='Fix all runs in a scenario (e.g., miller_urey_extended)'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Fix all runs in all scenarios'
    )
    parser.add_argument(
        '--base-dir',
        type=str,
        default='results/phase2b_additional',
        help='Base directory for results (default: results/phase2b_additional)'
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress detailed output'
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("FIXING PHASE 2B MOLECULES")
    print("=" * 70)
    print(f"Base directory: {args.base_dir}\n")
    
    base_dir = Path(args.base_dir)
    
    if args.run:
        # Fix single run
        run_dir = Path(args.run)
        success = fix_single_run(run_dir, verbose=not args.quiet)
        sys.exit(0 if success else 1)
    
    elif args.scenario:
        # Fix all runs in scenario
        results = fix_all_runs(base_dir, scenario=args.scenario, verbose=not args.quiet)
    
    elif args.all:
        # Fix all runs in all scenarios
        results = fix_all_runs(base_dir, scenario=None, verbose=not args.quiet)
    
    else:
        parser.print_help()
        sys.exit(1)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total runs processed: {results['total_runs']}")
    print(f"Successful: {results['successful']}")
    print(f"Failed: {results['failed']}")
    print("=" * 70)
    
    if results['failed'] > 0:
        print("\nFailed runs:")
        for run_info in results['runs']:
            if run_info['status'] == 'failed':
                print(f"  - {run_info['scenario']}/{run_info['run']}")
    
    sys.exit(0 if results['failed'] == 0 else 1)

if __name__ == "__main__":
    main()

