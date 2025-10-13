"""
Match Confidence Evaluator
===========================

Evaluates the confidence and reliability of molecule matches from simulation
to PubChem database.

Provides:
- Confidence scoring (0-1)
- Reliability classification (high/medium/low)
- Chemical plausibility checks
- Validation warnings

Used by PubChem Matcher v2 to assess match quality.
"""

import logging
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class Reliability(Enum):
    """Match reliability levels"""
    HIGH = "high"       # >0.8 confidence, chemically plausible
    MEDIUM = "medium"   # 0.5-0.8 confidence
    LOW = "low"         # <0.5 confidence or chemical issues
    INVALID = "invalid" # Chemically impossible


@dataclass
class MatchConfidence:
    """Result of confidence evaluation"""
    confidence_score: float  # 0-1
    reliability: Reliability
    warnings: List[str]
    validation_status: str  # "PASS", "WARNING", "FAIL"
    
    # Detailed metrics
    topology_score: float
    fingerprint_score: float
    energy_score: float
    spectral_score: float
    geometric_score: float
    
    # Chemical plausibility
    valence_check: bool
    charge_balance: bool
    bond_order_check: bool
    
    def __repr__(self):
        return (f"MatchConfidence(score={self.confidence_score:.3f}, "
                f"reliability={self.reliability.value}, "
                f"status={self.validation_status}, "
                f"warnings={len(self.warnings)})")


class MatchConfidenceEvaluator:
    """
    Evaluates confidence of molecule matches
    
    Usage:
        evaluator = MatchConfidenceEvaluator()
        confidence = evaluator.evaluate_match(cluster, match_result)
        
        if confidence.reliability == Reliability.HIGH:
            print(f"High confidence match: {confidence.confidence_score:.3f}")
    """
    
    def __init__(self,
                 confidence_threshold_high: float = 0.8,
                 confidence_threshold_medium: float = 0.5):
        """
        Initialize evaluator
        
        Args:
            confidence_threshold_high: Threshold for HIGH reliability
            confidence_threshold_medium: Threshold for MEDIUM reliability
        """
        self.threshold_high = confidence_threshold_high
        self.threshold_medium = confidence_threshold_medium
    
    def evaluate_match(self,
                      cluster: Dict,
                      match_result: Dict,
                      similarity_score: Optional['SimilarityScore'] = None) -> MatchConfidence:
        """
        Evaluate confidence of a match
        
        Args:
            cluster: Cluster dict from simulation {formula, atoms, bonds, ...}
            match_result: PubChem match result {cid, name, smiles, ...}
            similarity_score: SimilarityScore from MultiMetricSimilarity
        
        Returns:
            MatchConfidence with detailed evaluation
        """
        warnings = []
        
        # 1. Check chemical plausibility of cluster
        valence_ok = self.check_valence(cluster)
        charge_ok = self.check_charge_balance(cluster)
        bond_order_ok = self.check_bond_orders(cluster)
        
        if not valence_ok:
            warnings.append("Valence violations detected in cluster")
        if not charge_ok:
            warnings.append("Charge imbalance detected")
        if not bond_order_ok:
            warnings.append("Unusual bond orders detected")
        
        # 2. Extract similarity metrics
        if similarity_score:
            topo = similarity_score.topology
            fp = similarity_score.fingerprint
            energy = similarity_score.energy
            spectral = similarity_score.spectral
            geo = similarity_score.geometric
            overall = similarity_score.overall
        else:
            # Fallback to match_result scores if available
            topo = match_result.get('topology_score', 0.5)
            fp = match_result.get('fingerprint_score', 0.5)
            energy = match_result.get('energy_score', 0.5)
            spectral = match_result.get('spectral_score', 0.5)
            geo = match_result.get('geometric_score', 0.5)
            overall = match_result.get('similarity', 0.5)
        
        # 3. Check consistency between metrics
        metric_variance = self._compute_metric_variance([topo, fp, energy, spectral, geo])
        if metric_variance > 0.3:
            warnings.append(f"High variance between metrics ({metric_variance:.2f})")
        
        # 4. Check formula match
        cluster_formula = cluster.get('formula', '')
        match_formula = match_result.get('molecular_formula', '')
        if cluster_formula and match_formula and cluster_formula != match_formula:
            warnings.append(f"Formula mismatch: {cluster_formula} vs {match_formula}")
            overall *= 0.5  # Heavy penalty for formula mismatch
        
        # 5. Check atom count
        cluster_atoms = len(cluster.get('atoms', []))
        match_atoms = match_result.get('heavy_atom_count', 0)
        if cluster_atoms > 0 and match_atoms > 0:
            atom_diff = abs(cluster_atoms - match_atoms) / max(cluster_atoms, match_atoms)
            if atom_diff > 0.2:
                warnings.append(f"Atom count differs by {atom_diff*100:.1f}%")
        
        # 6. Determine reliability
        if not valence_ok or not charge_ok:
            reliability = Reliability.INVALID
            validation_status = "FAIL"
        elif overall >= self.threshold_high and len(warnings) <= 1:
            reliability = Reliability.HIGH
            validation_status = "PASS"
        elif overall >= self.threshold_medium:
            reliability = Reliability.MEDIUM
            validation_status = "WARNING" if warnings else "PASS"
        else:
            reliability = Reliability.LOW
            validation_status = "WARNING"
        
        return MatchConfidence(
            confidence_score=overall,
            reliability=reliability,
            warnings=warnings,
            validation_status=validation_status,
            topology_score=topo,
            fingerprint_score=fp,
            energy_score=energy,
            spectral_score=spectral,
            geometric_score=geo,
            valence_check=valence_ok,
            charge_balance=charge_ok,
            bond_order_check=bond_order_ok
        )
    
    def check_valence(self, cluster: Dict) -> bool:
        """
        Check for valence violations
        
        Returns True if all atoms have reasonable valence
        """
        atoms = cluster.get('atoms', [])
        bonds = cluster.get('bonds', [])
        
        if not atoms or not bonds:
            return True  # Can't validate
        
        # Count bonds per atom
        bond_count = {i: 0 for i in range(len(atoms))}
        for bond in bonds:
            i, j = bond[0], bond[1]
            bond_order = bond[2] if len(bond) > 2 else 1
            bond_count[i] += bond_order
            bond_count[j] += bond_order
        
        # Check valence limits
        max_valence = {
            'H': 1, 'C': 4, 'N': 3, 'O': 2,
            'P': 5, 'S': 6, 'F': 1, 'Cl': 1,
            'Br': 1, 'I': 1
        }
        
        for i, atom in enumerate(atoms):
            element = atom if isinstance(atom, str) else atom.get('element', 'C')
            max_val = max_valence.get(element, 4)
            
            if bond_count[i] > max_val:
                logger.warning(f"Valence violation: {element} has {bond_count[i]} bonds (max {max_val})")
                return False
        
        return True
    
    def check_charge_balance(self, cluster: Dict) -> bool:
        """
        Check for charge balance
        
        Returns True if total charge is reasonable
        """
        atoms = cluster.get('atoms', [])
        
        total_charge = 0
        for atom in atoms:
            if isinstance(atom, dict):
                charge = atom.get('charge', 0)
                total_charge += charge
        
        # Allow ±2 total charge
        if abs(total_charge) > 2:
            logger.warning(f"High total charge: {total_charge}")
            return False
        
        return True
    
    def check_bond_orders(self, cluster: Dict) -> bool:
        """
        Check for unusual bond orders
        
        Returns True if all bond orders are 1, 2, or 3
        """
        bonds = cluster.get('bonds', [])
        
        for bond in bonds:
            if len(bond) > 2:
                bond_order = bond[2]
                if bond_order < 1 or bond_order > 3:
                    logger.warning(f"Unusual bond order: {bond_order}")
                    return False
        
        return True
    
    def is_chemically_plausible(self, cluster: Dict) -> Tuple[bool, List[str]]:
        """
        Check if cluster is chemically plausible
        
        Returns:
            (is_plausible, list_of_issues)
        """
        issues = []
        
        # Check valence
        if not self.check_valence(cluster):
            issues.append("Valence violations")
        
        # Check charge
        if not self.check_charge_balance(cluster):
            issues.append("Charge imbalance")
        
        # Check bond orders
        if not self.check_bond_orders(cluster):
            issues.append("Invalid bond orders")
        
        # Check for disconnected atoms
        if self._has_disconnected_atoms(cluster):
            issues.append("Disconnected atoms")
        
        return (len(issues) == 0, issues)
    
    def _compute_metric_variance(self, metrics: List[float]) -> float:
        """Compute variance of metric scores"""
        import numpy as np
        return float(np.var(metrics))
    
    def _has_disconnected_atoms(self, cluster: Dict) -> bool:
        """Check if any atoms are disconnected"""
        atoms = cluster.get('atoms', [])
        bonds = cluster.get('bonds', [])
        
        if len(atoms) <= 1:
            return False
        
        if not bonds:
            return True  # Multiple atoms, no bonds
        
        # Build adjacency list
        adj = {i: [] for i in range(len(atoms))}
        for bond in bonds:
            i, j = bond[0], bond[1]
            adj[i].append(j)
            adj[j].append(i)
        
        # BFS to check connectivity
        visited = set()
        queue = [0]
        visited.add(0)
        
        while queue:
            node = queue.pop(0)
            for neighbor in adj[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        
        # All atoms should be reachable
        return len(visited) != len(atoms)


def generate_validation_report(confidence: MatchConfidence,
                               cluster: Dict,
                               match_result: Dict) -> str:
    """
    Generate human-readable validation report
    
    Args:
        confidence: MatchConfidence result
        cluster: Original cluster
        match_result: PubChem match
    
    Returns:
        Formatted report string
    """
    lines = []
    lines.append("=" * 70)
    lines.append("MATCH VALIDATION REPORT")
    lines.append("=" * 70)
    lines.append("")
    
    # Status
    status_icon = {
        "PASS": "✅",
        "WARNING": "⚠️",
        "FAIL": "❌"
    }
    lines.append(f"Status: {status_icon.get(confidence.validation_status, '?')} {confidence.validation_status}")
    lines.append(f"Reliability: {confidence.reliability.value.upper()}")
    lines.append(f"Confidence Score: {confidence.confidence_score:.3f}")
    lines.append("")
    
    # Cluster info
    lines.append("Cluster:")
    lines.append(f"  Formula: {cluster.get('formula', 'N/A')}")
    lines.append(f"  Atoms: {len(cluster.get('atoms', []))}")
    lines.append(f"  Bonds: {len(cluster.get('bonds', []))}")
    lines.append("")
    
    # Match info
    lines.append("PubChem Match:")
    lines.append(f"  CID: {match_result.get('cid', 'N/A')}")
    lines.append(f"  Name: {match_result.get('name', 'N/A')}")
    lines.append(f"  Formula: {match_result.get('molecular_formula', 'N/A')}")
    lines.append(f"  SMILES: {match_result.get('smiles', 'N/A')}")
    lines.append("")
    
    # Similarity metrics
    lines.append("Similarity Metrics:")
    lines.append(f"  Topology:    {confidence.topology_score:.3f}")
    lines.append(f"  Fingerprint: {confidence.fingerprint_score:.3f}")
    lines.append(f"  Energy:      {confidence.energy_score:.3f}")
    lines.append(f"  Spectral:    {confidence.spectral_score:.3f}")
    lines.append(f"  Geometric:   {confidence.geometric_score:.3f}")
    lines.append("")
    
    # Chemical checks
    lines.append("Chemical Plausibility:")
    lines.append(f"  Valence: {'✓' if confidence.valence_check else '✗'}")
    lines.append(f"  Charge:  {'✓' if confidence.charge_balance else '✗'}")
    lines.append(f"  Bonds:   {'✓' if confidence.bond_order_check else '✗'}")
    lines.append("")
    
    # Warnings
    if confidence.warnings:
        lines.append("Warnings:")
        for warning in confidence.warnings:
            lines.append(f"  ⚠ {warning}")
        lines.append("")
    
    lines.append("=" * 70)
    
    return "\n".join(lines)


# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("=" * 70)
    print("MATCH CONFIDENCE EVALUATOR - EXAMPLE")
    print("=" * 70)
    
    # Example 1: High confidence match
    print("\n\nExample 1: HIGH CONFIDENCE MATCH")
    print("-" * 70)
    
    cluster_good = {
        'formula': 'C2H4O2',
        'atoms': ['C', 'C', 'O', 'O', 'H', 'H', 'H', 'H'],
        'bonds': [(0, 1), (0, 2), (1, 3), (0, 4), (0, 5), (1, 6), (1, 7)]
    }
    
    match_good = {
        'cid': 757,
        'name': 'Glycolic acid',
        'molecular_formula': 'C2H4O2',
        'smiles': 'C(C(=O)O)O',
        'heavy_atom_count': 4,
        'similarity': 0.92,
        'topology_score': 0.95,
        'fingerprint_score': 0.90,
        'energy_score': 0.88,
        'spectral_score': 0.92,
        'geometric_score': 0.85
    }
    
    evaluator = MatchConfidenceEvaluator()
    confidence = evaluator.evaluate_match(cluster_good, match_good)
    
    report = generate_validation_report(confidence, cluster_good, match_good)
    print(report)
    
    # Example 2: Low confidence match
    print("\n\nExample 2: LOW CONFIDENCE MATCH")
    print("-" * 70)
    
    cluster_bad = {
        'formula': 'C2H4O2',
        'atoms': ['C', 'C', 'O', 'O', 'H', 'H', 'H', 'H'],
        'bonds': [(0, 1), (0, 2), (1, 3)]  # Missing many bonds!
    }
    
    match_bad = {
        'cid': 123,
        'name': 'Unknown compound',
        'molecular_formula': 'C3H6O',  # Different formula!
        'smiles': 'CCC=O',
        'heavy_atom_count': 4,
        'similarity': 0.35,
        'topology_score': 0.45,
        'fingerprint_score': 0.30,
        'energy_score': 0.40,
        'spectral_score': 0.25,
        'geometric_score': 0.35
    }
    
    confidence_bad = evaluator.evaluate_match(cluster_bad, match_bad)
    report_bad = generate_validation_report(confidence_bad, cluster_bad, match_bad)
    print(report_bad)
    
    print("\n" + "=" * 70)

