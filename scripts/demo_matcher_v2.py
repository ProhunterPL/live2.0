"""
Quick demo of PubChem Matcher v2
=================================

Demonstrates ML-based matching with confidence scoring.
"""

import sys
import json
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from matcher.matcher_v2 import MatcherV2
from matcher.confidence import generate_validation_report

# Test molecules
TEST_MOLECULES = {
    'water': {
        'formula': 'H2O',
        'atoms': ['O', 'H', 'H'],
        'bonds': [(0, 1), (0, 2)],
        'energy': -76.4
    },
    'formaldehyde': {
        'formula': 'CH2O',
        'atoms': ['C', 'H', 'H', 'O'],
        'bonds': [(0, 1), (0, 2), (0, 3, 2)],  # C=O double bond
        'energy': -114.2
    },
    'glycolic_acid': {
        'formula': 'C2H4O2',
        'atoms': ['C', 'C', 'O', 'O', 'H', 'H', 'H', 'H'],
        'bonds': [(0, 1), (0, 2), (1, 3, 2), (0, 4), (0, 5), (1, 6), (2, 7)],
        'energy': -150.5
    }
}


def demo_single_match():
    """Demo single molecule matching"""
    print("=" * 70)
    print("MATCHER V2 DEMO - Single Match")
    print("=" * 70)
    
    # Initialize matcher (without ML for speed)
    print("\n[1] Initializing MatcherV2...")
    matcher = MatcherV2(use_ml_classifier=False)
    print("    ✓ Matcher ready (ML disabled for demo)")
    
    # Test water molecule
    print("\n[2] Testing: Water (H2O)")
    print("-" * 70)
    cluster = TEST_MOLECULES['water']
    print(f"    Formula: {cluster['formula']}")
    print(f"    Atoms: {len(cluster['atoms'])}")
    print(f"    Bonds: {len(cluster['bonds'])}")
    
    try:
        print("\n[3] Matching to PubChem (requires internet)...")
        result = matcher.match_cluster(cluster, top_n=3, min_similarity=0.3)
        
        print("\n[4] Results:")
        print("-" * 70)
        
        if result.success:
            print(f"✅ Match found!")
            print(f"  PubChem CID: {result.pubchem_cid}")
            print(f"  Name: {result.pubchem_name}")
            print(f"  Formula: {result.pubchem_formula}")
            print(f"  SMILES: {result.pubchem_smiles}")
            print(f"")
            print(f"  Similarity: {result.similarity_score.overall:.3f}")
            print(f"    ├─ Topology:    {result.similarity_score.topology:.3f}")
            print(f"    ├─ Fingerprint: {result.similarity_score.fingerprint:.3f}")
            print(f"    ├─ Energy:      {result.similarity_score.energy:.3f}")
            print(f"    ├─ Spectral:    {result.similarity_score.spectral:.3f}")
            print(f"    └─ Geometric:   {result.similarity_score.geometric:.3f}")
            print(f"")
            print(f"  Confidence: {result.confidence.confidence_score:.3f}")
            print(f"  Reliability: {result.confidence.reliability.value.upper()}")
            print(f"  Status: {result.confidence.validation_status}")
            
            if result.confidence.warnings:
                print(f"\n  Warnings:")
                for warning in result.confidence.warnings:
                    print(f"    ⚠ {warning}")
        else:
            print(f"❌ Match failed: {result.error_message}")
        
        print("-" * 70)
        
    except Exception as e:
        print(f"❌ Error: {e}")
        print("(This demo requires internet connection for PubChem API)")


def demo_batch_match():
    """Demo batch matching"""
    print("\n\n")
    print("=" * 70)
    print("MATCHER V2 DEMO - Batch Match")
    print("=" * 70)
    
    print("\n[1] Initializing MatcherV2...")
    matcher = MatcherV2(use_ml_classifier=False)
    
    print("\n[2] Testing 3 molecules:")
    for name, cluster in TEST_MOLECULES.items():
        print(f"    - {name.replace('_', ' ').title()}: {cluster['formula']}")
    
    try:
        print("\n[3] Batch matching (requires internet)...")
        clusters = list(TEST_MOLECULES.values())
        results = matcher.match_batch(clusters, top_n=3, min_similarity=0.3)
        
        print("\n[4] Batch Results:")
        print("-" * 70)
        
        for i, (name, result) in enumerate(zip(TEST_MOLECULES.keys(), results)):
            print(f"\n{i+1}. {name.replace('_', ' ').title()}")
            if result.success:
                print(f"   ✅ {result.pubchem_name} (CID: {result.pubchem_cid})")
                print(f"   Similarity: {result.similarity_score.overall:.3f}, "
                      f"Confidence: {result.confidence.confidence_score:.3f}, "
                      f"Reliability: {result.confidence.reliability.value}")
            else:
                print(f"   ❌ Failed: {result.error_message}")
        
        # Summary
        successful = sum(1 for r in results if r.success)
        print(f"\n{'='*70}")
        print(f"Summary: {successful}/{len(results)} successful matches")
        print("=" * 70)
        
    except Exception as e:
        print(f"❌ Error: {e}")
        print("(This demo requires internet connection for PubChem API)")


def demo_confidence_checks():
    """Demo confidence evaluation (offline)"""
    print("\n\n")
    print("=" * 70)
    print("MATCHER V2 DEMO - Confidence Checks (Offline)")
    print("=" * 70)
    
    from matcher.confidence import MatchConfidenceEvaluator
    
    evaluator = MatchConfidenceEvaluator()
    
    # Test 1: Valid molecule
    print("\n[Test 1] Valid molecule (Water)")
    print("-" * 70)
    cluster = TEST_MOLECULES['water']
    is_plausible, issues = evaluator.is_chemically_plausible(cluster)
    
    print(f"Formula: {cluster['formula']}")
    print(f"Plausible: {is_plausible}")
    print(f"Issues: {issues if issues else 'None'}")
    print(f"✅ Chemical checks: valence={evaluator.check_valence(cluster)}, "
          f"charge={evaluator.check_charge_balance(cluster)}, "
          f"bonds={evaluator.check_bond_orders(cluster)}")
    
    # Test 2: Invalid molecule (too many bonds)
    print("\n\n[Test 2] Invalid molecule (Carbon with 5 bonds)")
    print("-" * 70)
    invalid_cluster = {
        'formula': 'CH5',
        'atoms': ['C', 'H', 'H', 'H', 'H', 'H'],
        'bonds': [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]  # 5 bonds!
    }
    
    is_plausible, issues = evaluator.is_chemically_plausible(invalid_cluster)
    
    print(f"Formula: {invalid_cluster['formula']}")
    print(f"Plausible: {is_plausible}")
    print(f"Issues: {issues}")
    print(f"❌ Chemical checks: valence={evaluator.check_valence(invalid_cluster)}, "
          f"charge={evaluator.check_charge_balance(invalid_cluster)}, "
          f"bonds={evaluator.check_bond_orders(invalid_cluster)}")
    
    print("\n" + "=" * 70)


def main():
    """Run all demos"""
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 15 + "PUBCHEM MATCHER V2 - DEMO" + " " * 28 + "║")
    print("╚" + "=" * 68 + "╝")
    print("")
    print("This demo showcases:")
    print("  1. Single molecule matching")
    print("  2. Batch matching (3 molecules)")
    print("  3. Confidence evaluation (offline)")
    print("")
    print("Note: Demos 1 & 2 require internet connection for PubChem API")
    print("")
    
    # Run confidence checks first (offline)
    demo_confidence_checks()
    
    # Ask before running online demos
    try:
        response = input("\nRun online demos (requires internet)? [y/N]: ")
        if response.lower() in ['y', 'yes']:
            demo_single_match()
            demo_batch_match()
        else:
            print("\nSkipping online demos.")
    except (KeyboardInterrupt, EOFError):
        print("\n\nDemo cancelled.")
    
    print("\n")
    print("=" * 70)
    print("DEMO COMPLETE!")
    print("=" * 70)
    print("\nFor more information, see:")
    print("  - docs/MATCHER_V2.md (documentation)")
    print("  - tests/test_matcher_v2.py (15 tests)")
    print("  - matcher/matcher_v2.py (source code)")
    print("")


if __name__ == "__main__":
    main()

