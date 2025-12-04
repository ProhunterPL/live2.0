#!/usr/bin/env python3
"""
Validate Novel Molecules with TruthFilter 2.0
==============================================

Validates the 5 novel molecules from Figure 6B using TruthFilterV2.

Usage:
    python scripts/validate_novel_molecules_tf2.py
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from backend.validation.truth_filter_v2 import TruthFilterV2
import json

# Novel molecules from Figure 6B
NOVEL_MOLECULES = [
    {
        'formula': 'C8H12N2O3',
        'mass': 184,
        'name': 'Novel compound 1',
        'smiles': 'CC(=O)NC1CCCC1NC(=O)C',  # diketopiperazine-like
        'metadata': {
            'occurrence_count': 1,
            'steps_seen': [100000],
            'pubchem_match': False
        }
    },
    {
        'formula': 'C7H9NO4',
        'mass': 171,
        'name': 'Novel compound 2',
        'smiles': 'CC(=O)OC1=CC=CC=C1N',  # N-acetyl derivative (aromatic)
        'metadata': {
            'occurrence_count': 1,
            'steps_seen': [150000],
            'pubchem_match': False
        }
    },
    {
        'formula': 'C9H11N3O2',
        'mass': 193,
        'name': 'Novel compound 3',
        'smiles': 'CC1=CC=C(C=C1)N(C)C(=O)N',  # N-methyl derivative (aromatic)
        'metadata': {
            'occurrence_count': 1,
            'steps_seen': [200000],
            'pubchem_match': False
        }
    },
    {
        'formula': 'C6H8N2O3',
        'mass': 156,
        'name': 'Novel compound 4',
        'smiles': 'CC(=O)NC1CCCC1N',  # cyclic amide
        'metadata': {
            'occurrence_count': 1,
            'steps_seen': [120000],
            'pubchem_match': False
        }
    },
    {
        'formula': 'C10H14NO2',
        'mass': 180,
        'name': 'Novel compound 5',
        'smiles': 'CC1=CC=C(C=C1)OC(=O)NC',  # aromatic ester
        'metadata': {
            'occurrence_count': 1,
            'steps_seen': [180000],
            'pubchem_match': False
        }
    },
]


def main():
    print("=" * 70)
    print("TRUTHFILTER 2.0 - Novel Molecules Validation")
    print("=" * 70)
    print()
    
    filter_v2 = TruthFilterV2()
    
    results = []
    for mol_data in NOVEL_MOLECULES:
        print(f"\nValidating: {mol_data['formula']} ({mol_data['name']})")
        print(f"  SMILES: {mol_data['smiles']}")
        
        result = filter_v2.validate_molecule(
            smiles=mol_data['smiles'],
            formula=mol_data['formula'],
            mass=mol_data['mass'],
            metadata=mol_data['metadata']
        )
        
        results.append({
            **mol_data,
            'validation': result
        })
        
        print(f"  Validity: {result['validity']}")
        print(f"  Confidence: {result['confidence']:.2f}")
        print(f"  Model Compatibility: {result['model_compatibility']}")
        print(f"  Reasons: {', '.join(result['reasons'])}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    accept_count = sum(1 for r in results if r['validation']['validity'] == 'ACCEPT')
    flag_count = sum(1 for r in results if r['validation']['validity'] == 'FLAG')
    reject_count = sum(1 for r in results if r['validation']['validity'] == 'REJECT')
    
    print(f"ACCEPT: {accept_count}")
    print(f"FLAG: {flag_count}")
    print(f"REJECT: {reject_count}")
    
    print("\nDetailed Results:")
    for r in results:
        val = r['validation']
        print(f"  {r['formula']}: {val['validity']} (confidence: {val['confidence']:.2f}, "
              f"compatibility: {val['model_compatibility']})")
        if 'AROMATIC_UNSUPPORTED_BY_MODEL' in val['reasons']:
            print(f"    ⚠️  Aromatic - model incompatible")
        if 'HIGH_STRAIN_RING' in val['reasons']:
            print(f"    ⚠️  High strain ring")
    
    # Save results
    output_file = project_root / 'paper' / 'novel_molecules_tf2_validation.json'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to: {output_file}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

