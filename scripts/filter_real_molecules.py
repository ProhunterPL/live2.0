#!/usr/bin/env python3
"""
Filter real chemical molecules from clusters/aggregates
And prepare for PubChem matching
"""
import json
from pathlib import Path
from collections import Counter

def is_real_molecule(mol):
    """
    Determine if a molecule is a real chemical species or just a cluster
    
    Criteria for REAL molecule:
    1. Size >= 2 (at least a dimer)
    2. Bonds >= 0.4 * (size - 1)  [relaxed from 0.5 for some branching]
    3. Not too sparse: bonds/size >= 0.3
    """
    size = mol.get('size', 0)
    n_bonds = len(mol.get('bonds', []))
    
    if size < 2:
        return False
    
    # Minimum bonds for real molecule (linear would have size-1)
    min_bonds_linear = 0.4 * (size - 1)
    
    # Bond density
    bond_density = n_bonds / size if size > 0 else 0
    
    # Real molecule criteria
    is_real = (n_bonds >= min_bonds_linear) and (bond_density >= 0.3)
    
    return is_real

def filter_molecules(input_file, output_file):
    """Filter batch analysis to keep only real molecules"""
    
    print("Loading data...")
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    print(f"Original runs: {len(data['runs'])}")
    
    # Statistics
    total_original = 0
    total_filtered = 0
    total_clusters = 0
    
    filtered_runs = []
    real_molecules = Counter()
    
    for run in data['runs']:
        total_original += len(run.get('molecules', []))
        
        # Filter molecules
        real_mols = []
        clusters = []
        
        for mol in run.get('molecules', []):
            if is_real_molecule(mol):
                real_mols.append(mol)
                formula = mol.get('formula', 'Unknown')
                count = mol.get('count', 1)
                real_molecules[formula] += count
            else:
                clusters.append(mol)
        
        total_filtered += len(real_mols)
        total_clusters += len(clusters)
        
        # Create filtered run
        filtered_run = run.copy()
        filtered_run['molecules'] = real_mols
        filtered_run['molecules_original'] = len(run.get('molecules', []))
        filtered_run['molecules_filtered'] = len(real_mols)
        filtered_run['clusters_removed'] = len(clusters)
        filtered_runs.append(filtered_run)
    
    # Create filtered data
    filtered_data = {
        'runs': filtered_runs,
        'summary': {
            'total_runs': len(filtered_runs),
            'successful_runs': len(filtered_runs),
            'total_unique_molecules': len(real_molecules),
            'total_instances': sum(real_molecules.values()),
            'filtering_stats': {
                'original_molecules': total_original,
                'filtered_molecules': total_filtered,
                'clusters_removed': total_clusters,
                'retention_rate': f"{100 * total_filtered / total_original:.1f}%" if total_original > 0 else "0%"
            }
        },
        'most_common_real_molecules': [
            {'formula': formula, 'count': count}
            for formula, count in real_molecules.most_common(50)
        ]
    }
    
    # Save filtered data
    with open(output_file, 'w') as f:
        json.dump(filtered_data, f, indent=2)
    
    # Print report
    print()
    print("=" * 80)
    print("MOLECULE FILTERING REPORT")
    print("=" * 80)
    print()
    print(f"Original molecules (all):        {total_original}")
    print(f"Real molecules (filtered):       {total_filtered}")
    print(f"Clusters removed:                {total_clusters}")
    print(f"Retention rate:                  {100 * total_filtered / total_original:.1f}%")
    print()
    print(f"Unique real molecules:           {len(real_molecules)}")
    print(f"Total instances (real):          {sum(real_molecules.values())}")
    print()
    print("Top 20 Real Molecules:")
    for i, (formula, count) in enumerate(real_molecules.most_common(20), 1):
        print(f"  {i:2d}. {formula[:40]:<40} : {count:6d} instances")
    print()
    print(f"Filtered data saved to: {output_file}")
    print("=" * 80)
    
    return filtered_data

def main():
    input_file = Path("analysis/phase2b_miller_urey/batch_analysis.json")
    output_file = Path("analysis/phase2b_miller_urey/batch_analysis_filtered.json")
    
    filtered_data = filter_molecules(input_file, output_file)
    
    print()
    print("Next step: Use matcher to identify molecules")
    print("Command: python scripts/match_molecules_pubchem.py")

if __name__ == "__main__":
    main()

