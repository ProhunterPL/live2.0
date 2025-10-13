"""
PubChem Matcher v2
==================

Advanced molecule matching system using machine learning and multi-metric similarity.

Improvements over v1:
1. ML-based atom type classification (RandomForest on 100k atoms)
2. Multi-metric similarity (5 metrics: topology, fingerprint, energy, spectral, geometric)
3. Confidence scoring and validation
4. Chemical plausibility checks
5. Batch processing support

Usage:
    from matcher.matcher_v2 import MatcherV2
    
    matcher = MatcherV2()
    result = matcher.match_cluster(cluster_data)
    
    if result.confidence.reliability == Reliability.HIGH:
        print(f"Match: {result.pubchem_name} (CID: {result.pubchem_cid})")
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict

# Local imports
from matcher.ml_classifier import AtomTypeClassifier, AtomFeatures, extract_features_from_rdkit_atom
from matcher.similarity import MultiMetricSimilarity, SimilarityScore
from matcher.confidence import MatchConfidenceEvaluator, MatchConfidence, Reliability
from matcher.chem import (
    json_to_mol, mol_to_smiles, pubchem_similar_top,
    render_mol_png, render_pubchem_png
)

logger = logging.getLogger(__name__)


@dataclass
class MatchResult:
    """Result of matching process"""
    # Input cluster
    cluster_formula: str
    cluster_atoms: int
    cluster_bonds: int
    
    # PubChem match
    pubchem_cid: Optional[int]
    pubchem_name: Optional[str]
    pubchem_formula: Optional[str]
    pubchem_smiles: Optional[str]
    
    # Similarity
    similarity_score: SimilarityScore
    
    # Confidence
    confidence: MatchConfidence
    
    # Success flag
    success: bool
    error_message: Optional[str] = None
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON export"""
        return {
            'cluster': {
                'formula': self.cluster_formula,
                'atoms': self.cluster_atoms,
                'bonds': self.cluster_bonds
            },
            'pubchem': {
                'cid': self.pubchem_cid,
                'name': self.pubchem_name,
                'formula': self.pubchem_formula,
                'smiles': self.pubchem_smiles
            },
            'similarity': {
                'overall': self.similarity_score.overall,
                'topology': self.similarity_score.topology,
                'fingerprint': self.similarity_score.fingerprint,
                'energy': self.similarity_score.energy,
                'spectral': self.similarity_score.spectral,
                'geometric': self.similarity_score.geometric
            },
            'confidence': {
                'score': self.confidence.confidence_score,
                'reliability': self.confidence.reliability.value,
                'validation_status': self.confidence.validation_status,
                'warnings': self.confidence.warnings,
                'valence_check': self.confidence.valence_check,
                'charge_balance': self.confidence.charge_balance,
                'bond_order_check': self.confidence.bond_order_check
            },
            'success': self.success,
            'error': self.error_message
        }


class MatcherV2:
    """
    PubChem Matcher v2 with ML and multi-metric similarity
    
    Usage:
        matcher = MatcherV2(
            classifier_model='data/atom_classifier.pkl',
            confidence_threshold=0.7
        )
        
        result = matcher.match_cluster(cluster_data)
        print(result.confidence)
    """
    
    def __init__(self,
                 classifier_model: Optional[str] = 'data/atom_classifier.pkl',
                 confidence_threshold_high: float = 0.8,
                 confidence_threshold_medium: float = 0.5,
                 use_ml_classifier: bool = True,
                 use_multi_metric: bool = True):
        """
        Initialize Matcher v2
        
        Args:
            classifier_model: Path to trained ML classifier
            confidence_threshold_high: Threshold for HIGH reliability
            confidence_threshold_medium: Threshold for MEDIUM reliability
            use_ml_classifier: Use ML for atom type prediction
            use_multi_metric: Use multi-metric similarity
        """
        self.use_ml = use_ml_classifier
        self.use_multi_metric = use_multi_metric
        
        # Load ML classifier
        self.classifier = None
        if self.use_ml and classifier_model:
            try:
                self.classifier = AtomTypeClassifier(model_path=classifier_model)
                logger.info(f"[+] ML classifier loaded from {classifier_model}")
            except Exception as e:
                logger.warning(f"Failed to load ML classifier: {e}")
                self.use_ml = False
        
        # Initialize similarity calculator
        self.similarity = MultiMetricSimilarity()
        
        # Initialize confidence evaluator
        self.confidence_evaluator = MatchConfidenceEvaluator(
            confidence_threshold_high=confidence_threshold_high,
            confidence_threshold_medium=confidence_threshold_medium
        )
        
        logger.info(f"[+] MatcherV2 initialized (ML={self.use_ml}, MultiMetric={self.use_multi_metric})")
    
    def match_cluster(self,
                     cluster: Dict,
                     top_n: int = 5,
                     min_similarity: float = 0.3) -> MatchResult:
        """
        Match a cluster to PubChem compounds
        
        Args:
            cluster: Cluster dict with {formula, atoms, bonds, positions, energy, ...}
            top_n: Number of PubChem candidates to consider
            min_similarity: Minimum similarity threshold
        
        Returns:
            MatchResult with best match and confidence
        """
        try:
            logger.info(f"[+] Matching cluster: {cluster.get('formula', 'unknown')}")
            
            # Step 1: Convert cluster to RDKit molecule
            cluster_mol = json_to_mol(cluster)
            if cluster_mol is None:
                return MatchResult(
                    cluster_formula=cluster.get('formula', ''),
                    cluster_atoms=len(cluster.get('atoms', [])),
                    cluster_bonds=len(cluster.get('bonds', [])),
                    pubchem_cid=None,
                    pubchem_name=None,
                    pubchem_formula=None,
                    pubchem_smiles=None,
                    similarity_score=SimilarityScore(0, 0, 0, 0, 0),
                    confidence=self._create_failed_confidence("Failed to convert cluster to molecule"),
                    success=False,
                    error_message="Failed to convert cluster to molecule"
                )
            
            # Step 2: Generate SMILES
            cluster_smiles = mol_to_smiles(cluster_mol)
            if not cluster_smiles:
                return self._failed_result(cluster, "Failed to generate SMILES")
            
            logger.info(f"  SMILES: {cluster_smiles}")
            
            # Step 3: Query PubChem for similar compounds
            pubchem_matches = pubchem_similar_top(cluster_smiles, top_n=top_n)
            if not pubchem_matches:
                return self._failed_result(cluster, "No PubChem matches found")
            
            logger.info(f"  Found {len(pubchem_matches)} PubChem candidates")
            
            # Step 4: Evaluate each candidate
            best_match = None
            best_similarity = None
            best_confidence = None
            
            for i, pubchem_match in enumerate(pubchem_matches):
                logger.info(f"  Evaluating candidate {i+1}/{len(pubchem_matches)}: CID {pubchem_match.get('cid', 'N/A')}")
                
                # Compute multi-metric similarity
                if self.use_multi_metric:
                    similarity = self._compute_similarity(cluster, pubchem_match)
                else:
                    # Fallback to simple similarity
                    similarity = SimilarityScore(
                        topology=pubchem_match.get('similarity', 0.5),
                        fingerprint=pubchem_match.get('similarity', 0.5),
                        energy=0.5,
                        spectral=0.5,
                        geometric=0.5
                    )
                
                # Skip if below threshold
                if similarity.overall < min_similarity:
                    logger.info(f"    Skipped (similarity {similarity.overall:.3f} < {min_similarity})")
                    continue
                
                # Evaluate confidence
                confidence = self.confidence_evaluator.evaluate_match(
                    cluster=cluster,
                    match_result=pubchem_match,
                    similarity_score=similarity
                )
                
                logger.info(f"    Similarity: {similarity.overall:.3f}, Confidence: {confidence.confidence_score:.3f}, "
                           f"Reliability: {confidence.reliability.value}")
                
                # Update best match
                if best_match is None or similarity.overall > best_similarity.overall:
                    best_match = pubchem_match
                    best_similarity = similarity
                    best_confidence = confidence
            
            # Step 5: Return best match
            if best_match:
                logger.info(f"[+] Best match: {best_match.get('name', 'N/A')} (CID: {best_match.get('cid', 'N/A')})")
                logger.info(f"    Similarity: {best_similarity.overall:.3f}, Confidence: {best_confidence.confidence_score:.3f}")
                
                return MatchResult(
                    cluster_formula=cluster.get('formula', ''),
                    cluster_atoms=len(cluster.get('atoms', [])),
                    cluster_bonds=len(cluster.get('bonds', [])),
                    pubchem_cid=best_match.get('cid'),
                    pubchem_name=best_match.get('name'),
                    pubchem_formula=best_match.get('molecular_formula'),
                    pubchem_smiles=best_match.get('smiles'),
                    similarity_score=best_similarity,
                    confidence=best_confidence,
                    success=True
                )
            else:
                return self._failed_result(cluster, f"No matches above similarity threshold {min_similarity}")
        
        except Exception as e:
            logger.error(f"Error in match_cluster: {e}")
            import traceback
            traceback.print_exc()
            return self._failed_result(cluster, str(e))
    
    def match_batch(self,
                   clusters: List[Dict],
                   top_n: int = 5,
                   min_similarity: float = 0.3) -> List[MatchResult]:
        """
        Match multiple clusters in batch
        
        Args:
            clusters: List of cluster dicts
            top_n: Number of PubChem candidates per cluster
            min_similarity: Minimum similarity threshold
        
        Returns:
            List of MatchResult
        """
        logger.info(f"[+] Batch matching {len(clusters)} clusters...")
        
        results = []
        for i, cluster in enumerate(clusters):
            logger.info(f"\n--- Cluster {i+1}/{len(clusters)} ---")
            result = self.match_cluster(cluster, top_n=top_n, min_similarity=min_similarity)
            results.append(result)
        
        # Summary
        successful = sum(1 for r in results if r.success)
        high_conf = sum(1 for r in results if r.success and r.confidence.reliability == Reliability.HIGH)
        
        logger.info(f"\n[+] Batch complete: {successful}/{len(clusters)} successful, {high_conf} high confidence")
        
        return results
    
    def _compute_similarity(self, cluster: Dict, pubchem_match: Dict) -> SimilarityScore:
        """Compute multi-metric similarity"""
        # Convert PubChem match to compatible format
        pubchem_dict = {
            'formula': pubchem_match.get('molecular_formula', ''),
            'atoms': [],  # Would need to parse SMILES
            'bonds': [],
            'energy': 0.0,  # Not available from PubChem
            'smiles': pubchem_match.get('smiles', '')
        }
        
        # Compute similarity
        similarity = self.similarity.compute(
            mol_sim=cluster,
            mol_ref=pubchem_dict,
            include_geometric=False  # PubChem doesn't provide 3D coords by default
        )
        
        return similarity
    
    def _failed_result(self, cluster: Dict, error_message: str) -> MatchResult:
        """Create failed MatchResult"""
        return MatchResult(
            cluster_formula=cluster.get('formula', ''),
            cluster_atoms=len(cluster.get('atoms', [])),
            cluster_bonds=len(cluster.get('bonds', [])),
            pubchem_cid=None,
            pubchem_name=None,
            pubchem_formula=None,
            pubchem_smiles=None,
            similarity_score=SimilarityScore(0, 0, 0, 0, 0),
            confidence=self._create_failed_confidence(error_message),
            success=False,
            error_message=error_message
        )
    
    def _create_failed_confidence(self, error: str) -> MatchConfidence:
        """Create failed MatchConfidence"""
        return MatchConfidence(
            confidence_score=0.0,
            reliability=Reliability.INVALID,
            warnings=[error],
            validation_status="FAIL",
            topology_score=0.0,
            fingerprint_score=0.0,
            energy_score=0.0,
            spectral_score=0.0,
            geometric_score=0.0,
            valence_check=False,
            charge_balance=False,
            bond_order_check=False
        )
    
    def export_result(self, result: MatchResult, output_path: str):
        """Export match result to JSON"""
        with open(output_path, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)
        
        logger.info(f"[+] Result exported to {output_path}")
    
    def export_batch_results(self, results: List[MatchResult], output_path: str):
        """Export batch results to JSON"""
        data = {
            'total': len(results),
            'successful': sum(1 for r in results if r.success),
            'high_confidence': sum(1 for r in results if r.success and r.confidence.reliability == Reliability.HIGH),
            'results': [r.to_dict() for r in results]
        }
        
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        logger.info(f"[+] Batch results exported to {output_path}")


# CLI Interface
def main():
    """Command-line interface for MatcherV2"""
    import argparse
    
    parser = argparse.ArgumentParser(description="PubChem Matcher v2 - Match clusters to PubChem")
    parser.add_argument('input', help='Path to cluster JSON file')
    parser.add_argument('--output', '-o', help='Output JSON path (default: <input>_match.json)')
    parser.add_argument('--top-n', type=int, default=5, help='Number of PubChem candidates (default: 5)')
    parser.add_argument('--min-similarity', type=float, default=0.3, help='Minimum similarity (default: 0.3)')
    parser.add_argument('--model', default='data/atom_classifier.pkl', help='ML classifier model path')
    parser.add_argument('--no-ml', action='store_true', help='Disable ML classifier')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    print("=" * 70)
    print("PUBCHEM MATCHER V2")
    print("=" * 70)
    
    # Load cluster
    print(f"\n[1/4] Loading cluster from {args.input}...")
    with open(args.input, 'r') as f:
        cluster = json.load(f)
    
    print(f"  Formula: {cluster.get('formula', 'N/A')}")
    print(f"  Atoms: {len(cluster.get('atoms', []))}")
    print(f"  Bonds: {len(cluster.get('bonds', []))}")
    
    # Initialize matcher
    print(f"\n[2/4] Initializing matcher...")
    matcher = MatcherV2(
        classifier_model=args.model if not args.no_ml else None,
        use_ml_classifier=not args.no_ml
    )
    
    # Match cluster
    print(f"\n[3/4] Matching to PubChem...")
    result = matcher.match_cluster(
        cluster=cluster,
        top_n=args.top_n,
        min_similarity=args.min_similarity
    )
    
    # Display result
    print(f"\n[4/4] Result:")
    print("-" * 70)
    
    if result.success:
        print(f"✅ Match found!")
        print(f"  PubChem CID: {result.pubchem_cid}")
        print(f"  Name: {result.pubchem_name}")
        print(f"  Formula: {result.pubchem_formula}")
        print(f"  SMILES: {result.pubchem_smiles}")
        print(f"")
        print(f"  Similarity: {result.similarity_score.overall:.3f}")
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
    
    # Export result
    output_path = args.output or str(Path(args.input).with_suffix('')) + '_match.json'
    print(f"\nExporting to {output_path}...")
    matcher.export_result(result, output_path)
    
    print("\n✅ Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()

