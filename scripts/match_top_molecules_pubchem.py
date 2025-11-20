#!/usr/bin/env python3
"""
Match top real molecules to PubChem database
"""
import sys
import json
from pathlib import Path
from datetime import datetime

# Add parent directory
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from matcher.matcher_v2 import MatcherV2
    MATCHER_AVAILABLE = True
except ImportError as e:
    print(f"Warning: MatcherV2 not available: {e}")
    MATCHER_AVAILABLE = False

def load_filtered_molecules(filtered_file):
    """Load filtered molecules"""
    with open(filtered_file, 'r') as f:
        data = json.load(f)
    return data

def prepare_molecule_for_matching(formula, runs_data):
    """
    Prepare molecule data for matching
    Find first occurrence in runs to get structure
    """
    # Find molecule in runs
    for run in runs_data['runs']:
        for mol in run.get('molecules', []):
            if mol.get('formula') == formula:
                # Prepare minimal cluster data for matcher
                return {
                    'formula': formula,
                    'atoms': mol.get('cluster', []),  # Particle indices
                    'bonds': mol.get('bonds', []),
                    'size': mol.get('size', 0)
                }
    return None

def match_top_molecules(filtered_file, batch_file, output_file, top_n=20):
    """Match top N real molecules to PubChem"""
    
    if not MATCHER_AVAILABLE:
        print("ERROR: MatcherV2 not available. Cannot proceed.")
        return None
    
    print("=" * 80)
    print("PUBCHEM MATCHING - TOP REAL MOLECULES")
    print("=" * 80)
    print()
    
    # Load data
    print("[1] Loading filtered molecules...")
    filtered_data = load_filtered_molecules(filtered_file)
    
    with open(batch_file, 'r') as f:
        batch_data = json.load(f)
    
    top_molecules = filtered_data.get('most_common_real_molecules', [])[:top_n]
    print(f"    Found {len(top_molecules)} molecules to match")
    print()
    
    # Initialize matcher
    print("[2] Initializing MatcherV2...")
    print("    (This may take a moment...)")
    matcher = MatcherV2(use_ml_classifier=False)  # Disable ML for speed
    print("    Matcher ready!")
    print()
    
    # Match each molecule
    print("[3] Matching molecules to PubChem...")
    print("    (Requires internet connection)")
    print("-" * 80)
    
    results = []
    for i, mol_data in enumerate(top_molecules, 1):
        formula = mol_data['formula']
        count = mol_data['count']
        
        print(f"\n[{i}/{len(top_molecules)}] {formula[:40]}")
        print(f"    Instances: {count}")
        
        # Prepare molecule
        mol_cluster = prepare_molecule_for_matching(formula, batch_data)
        
        if not mol_cluster:
            print("    ERROR: Could not find molecule structure")
            results.append({
                'formula': formula,
                'count': count,
                'success': False,
                'error': 'Structure not found'
            })
            continue
        
        try:
            # Match to PubChem
            result = matcher.match_cluster(mol_cluster, top_n=5, min_similarity=0.3)
            
            if result.success:
                print(f"    MATCH: {result.pubchem_name}")
                print(f"      CID: {result.pubchem_cid}")
                print(f"      Similarity: {result.similarity_score.overall:.3f}")
                print(f"      Confidence: {result.confidence.confidence_score:.3f} ({result.confidence.reliability.value.upper()})")
                
                results.append({
                    'formula': formula,
                    'count': count,
                    'success': True,
                    'pubchem_cid': result.pubchem_cid,
                    'pubchem_name': result.pubchem_name,
                    'pubchem_formula': result.pubchem_formula,
                    'pubchem_smiles': result.pubchem_smiles,
                    'similarity': result.similarity_score.overall,
                    'confidence': result.confidence.confidence_score,
                    'reliability': result.confidence.reliability.value
                })
            else:
                print(f"    NO MATCH: {result.error_message}")
                results.append({
                    'formula': formula,
                    'count': count,
                    'success': False,
                    'error': result.error_message
                })
        
        except Exception as e:
            print(f"    ERROR: {str(e)}")
            results.append({
                'formula': formula,
                'count': count,
                'success': False,
                'error': str(e)
            })
    
    print()
    print("-" * 80)
    print()
    
    # Summary
    successful = sum(1 for r in results if r.get('success', False))
    print(f"[4] Summary:")
    print(f"    Total molecules: {len(results)}")
    print(f"    Successful matches: {successful}")
    print(f"    Failed: {len(results) - successful}")
    print()
    
    # Save results
    output_data = {
        'timestamp': datetime.now().isoformat(),
        'total_molecules': len(results),
        'successful_matches': successful,
        'matches': results
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Results saved to: {output_file}")
    print()
    
    # Print identified molecules
    if successful > 0:
        print("=" * 80)
        print("IDENTIFIED MOLECULES")
        print("=" * 80)
        for r in results:
            if r.get('success'):
                print(f"\n{r['pubchem_name']} (CID: {r['pubchem_cid']})")
                print(f"  Formula: {r['pubchem_formula']}")
                print(f"  SMILES: {r['pubchem_smiles']}")
                print(f"  Similarity: {r['similarity']:.3f}")
                print(f"  Confidence: {r['confidence']:.3f} ({r['reliability']})")
                print(f"  Occurrences: {r['count']}")
        print()
        print("=" * 80)
    
    return output_data

def main():
    filtered_file = Path("analysis/phase2b_miller_urey/batch_analysis_filtered.json")
    batch_file = Path("analysis/phase2b_miller_urey/batch_analysis.json")
    output_file = Path("analysis/phase2b_miller_urey/pubchem_matches.json")
    
    if not filtered_file.exists():
        print(f"ERROR: Filtered file not found: {filtered_file}")
        print("Run: python scripts/filter_real_molecules.py first")
        return
    
    match_top_molecules(filtered_file, batch_file, output_file, top_n=20)

if __name__ == "__main__":
    main()

