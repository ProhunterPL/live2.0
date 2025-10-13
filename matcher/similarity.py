"""
Multi-Metric Molecular Similarity
==================================

Computes molecular similarity using 5 complementary metrics:

1. **Topology**: Graph similarity (isomorphism, subgraph)
2. **Fingerprint**: Chemical fingerprint similarity (Tanimoto)
3. **Energy**: Relative energy difference
4. **Spectral**: Graph spectral similarity (eigenvalues)
5. **Geometric**: 3D structure similarity (RMSD)

Combined score provides robust molecule matching for PubChem Matcher v2.
"""

import logging
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class SimilarityScore:
    """Multi-metric similarity score"""
    topology: float  # 0-1
    fingerprint: float  # 0-1
    energy: float  # 0-1
    spectral: float  # 0-1
    geometric: float  # 0-1
    
    @property
    def overall(self) -> float:
        """Weighted average similarity"""
        weights = {
            'topology': 0.30,
            'fingerprint': 0.30,
            'energy': 0.15,
            'spectral': 0.15,
            'geometric': 0.10
        }
        return (
            weights['topology'] * self.topology +
            weights['fingerprint'] * self.fingerprint +
            weights['energy'] * self.energy +
            weights['spectral'] * self.spectral +
            weights['geometric'] * self.geometric
        )
    
    def __repr__(self):
        return (f"SimilarityScore(overall={self.overall:.3f}, "
                f"topo={self.topology:.3f}, fp={self.fingerprint:.3f}, "
                f"E={self.energy:.3f}, spec={self.spectral:.3f}, geo={self.geometric:.3f})")


class MultiMetricSimilarity:
    """
    Computes molecular similarity using 5 metrics
    
    Usage:
        similarity = MultiMetricSimilarity()
        score = similarity.compute(mol_sim, mol_ref)
        
        if score.overall > 0.8:
            print(f"Good match! {score}")
    """
    
    def __init__(self):
        """Initialize similarity calculator"""
        self.method_available = {
            'rdkit': self._check_rdkit(),
            'networkx': self._check_networkx(),
            'scipy': self._check_scipy()
        }
    
    def _check_rdkit(self) -> bool:
        """Check if RDKit is available"""
        try:
            import rdkit
            return True
        except ImportError:
            logger.warning("RDKit not available - fingerprint metric disabled")
            return False
    
    def _check_networkx(self) -> bool:
        """Check if NetworkX is available"""
        try:
            import networkx
            return True
        except ImportError:
            logger.warning("NetworkX not available - topology metric may be limited")
            return False
    
    def _check_scipy(self) -> bool:
        """Check if SciPy is available"""
        try:
            import scipy
            return True
        except ImportError:
            logger.warning("SciPy not available - spectral metric disabled")
            return False
    
    def compute(self,
               mol_sim: Dict,
               mol_ref: Dict,
               include_geometric: bool = True) -> SimilarityScore:
        """
        Compute multi-metric similarity
        
        Args:
            mol_sim: Simulated molecule dict {formula, atoms, bonds, positions, energy}
            mol_ref: Reference molecule dict (same format or RDKit mol)
            include_geometric: Include geometric (RMSD) metric (requires 3D coords)
        
        Returns:
            SimilarityScore with all metrics
        """
        # Compute each metric
        topo_score = self.topology_similarity(mol_sim, mol_ref)
        fp_score = self.fingerprint_similarity(mol_sim, mol_ref)
        energy_score = self.energy_similarity(mol_sim, mol_ref)
        spectral_score = self.spectral_similarity(mol_sim, mol_ref)
        
        if include_geometric and 'positions' in mol_sim and 'positions' in mol_ref:
            geo_score = self.geometric_similarity(mol_sim, mol_ref)
        else:
            geo_score = 0.5  # Neutral if not available
        
        return SimilarityScore(
            topology=topo_score,
            fingerprint=fp_score,
            energy=energy_score,
            spectral=spectral_score,
            geometric=geo_score
        )
    
    def topology_similarity(self, mol_sim: Dict, mol_ref: Dict) -> float:
        """
        Metric 1: Topology Similarity
        
        Compares molecular graphs (atoms as nodes, bonds as edges).
        Uses graph isomorphism / subgraph matching.
        
        Returns similarity [0, 1]
        """
        # Simple version: compare formula and bond count
        # Full version would use graph isomorphism (NetworkX)
        
        formula_sim = mol_sim.get('formula', '')
        formula_ref = mol_ref.get('formula', '')
        
        if formula_sim == formula_ref:
            formula_match = 1.0
        else:
            # Partial match based on atom counts
            formula_match = self._formula_similarity(formula_sim, formula_ref)
        
        # Compare number of bonds
        n_bonds_sim = len(mol_sim.get('bonds', []))
        n_bonds_ref = len(mol_ref.get('bonds', []))
        
        if n_bonds_ref > 0:
            bond_ratio = min(n_bonds_sim, n_bonds_ref) / max(n_bonds_sim, n_bonds_ref)
        else:
            bond_ratio = 1.0 if n_bonds_sim == 0 else 0.0
        
        # Combine
        topology_score = 0.7 * formula_match + 0.3 * bond_ratio
        
        return topology_score
    
    def fingerprint_similarity(self, mol_sim: Dict, mol_ref: Dict) -> float:
        """
        Metric 2: Chemical Fingerprint Similarity
        
        Uses molecular fingerprints (Morgan/ECFP4) and Tanimoto similarity.
        
        Returns similarity [0, 1]
        """
        if not self.method_available['rdkit']:
            # Fallback: Use simple atom/bond counts
            return self._simple_fingerprint(mol_sim, mol_ref)
        
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from rdkit import DataStructs
            
            # Convert to RDKit mols (if not already)
            rdkit_sim = self._to_rdkit_mol(mol_sim)
            rdkit_ref = self._to_rdkit_mol(mol_ref)
            
            if rdkit_sim is None or rdkit_ref is None:
                return self._simple_fingerprint(mol_sim, mol_ref)
            
            # Generate Morgan fingerprints (ECFP4)
            fp_sim = AllChem.GetMorganFingerprint(rdkit_sim, 2)
            fp_ref = AllChem.GetMorganFingerprint(rdkit_ref, 2)
            
            # Tanimoto similarity
            tanimoto = DataStructs.TanimotoSimilarity(fp_sim, fp_ref)
            
            return tanimoto
        
        except Exception as e:
            logger.warning(f"Fingerprint calculation failed: {e}")
            return self._simple_fingerprint(mol_sim, mol_ref)
    
    def energy_similarity(self, mol_sim: Dict, mol_ref: Dict) -> float:
        """
        Metric 3: Energy Similarity
        
        Compares relative energies. Similar molecules should have similar energies.
        
        Returns similarity [0, 1]
        """
        energy_sim = mol_sim.get('energy', 0.0)
        energy_ref = mol_ref.get('energy', 0.0)
        
        # Normalize by number of atoms
        n_atoms_sim = len(mol_sim.get('atoms', []))
        n_atoms_ref = len(mol_ref.get('atoms', []))
        
        if n_atoms_sim > 0 and n_atoms_ref > 0:
            energy_per_atom_sim = energy_sim / n_atoms_sim
            energy_per_atom_ref = energy_ref / n_atoms_ref
            
            # Relative difference
            if abs(energy_per_atom_ref) > 1e-6:
                rel_diff = abs(energy_per_atom_sim - energy_per_atom_ref) / abs(energy_per_atom_ref)
            else:
                rel_diff = abs(energy_per_atom_sim - energy_per_atom_ref)
            
            # Convert to similarity (closer energies → higher score)
            # Allow up to 20% difference for "similar"
            energy_score = max(0.0, 1.0 - rel_diff / 0.2)
        else:
            energy_score = 0.5  # Neutral if no energy data
        
        return energy_score
    
    def spectral_similarity(self, mol_sim: Dict, mol_ref: Dict) -> float:
        """
        Metric 4: Graph Spectral Similarity
        
        Compares eigenvalues of molecular graph Laplacian.
        Graph spectrum is a molecular "fingerprint".
        
        Returns similarity [0, 1]
        """
        if not self.method_available['scipy']:
            return 0.5  # Neutral fallback
        
        try:
            # Build adjacency matrices
            adj_sim = self._build_adjacency_matrix(mol_sim)
            adj_ref = self._build_adjacency_matrix(mol_ref)
            
            # Compute Laplacian eigenvalues
            eigs_sim = self._compute_graph_spectrum(adj_sim)
            eigs_ref = self._compute_graph_spectrum(adj_ref)
            
            # Compare spectra
            spectral_score = self._compare_spectra(eigs_sim, eigs_ref)
            
            return spectral_score
        
        except Exception as e:
            logger.warning(f"Spectral similarity failed: {e}")
            return 0.5
    
    def geometric_similarity(self, mol_sim: Dict, mol_ref: Dict) -> float:
        """
        Metric 5: Geometric Similarity
        
        Compares 3D structures using RMSD (Root Mean Square Deviation).
        
        Returns similarity [0, 1]
        """
        positions_sim = mol_sim.get('positions', None)
        positions_ref = mol_ref.get('positions', None)
        
        if positions_sim is None or positions_ref is None:
            return 0.5  # Neutral if no 3D coords
        
        # Convert to numpy arrays
        pos_sim = np.array(positions_sim)
        pos_ref = np.array(positions_ref)
        
        # Must have same number of atoms
        if pos_sim.shape[0] != pos_ref.shape[0]:
            return 0.0
        
        # Compute RMSD (after optimal alignment)
        rmsd = self._compute_rmsd(pos_sim, pos_ref)
        
        # Convert RMSD to similarity
        # RMSD < 0.5 Å → very similar (score ~1.0)
        # RMSD > 2.0 Å → different (score ~0.0)
        geometric_score = max(0.0, 1.0 - rmsd / 2.0)
        
        return geometric_score
    
    # Helper methods
    
    def _formula_similarity(self, formula1: str, formula2: str) -> float:
        """Compare chemical formulas (simple version)"""
        # Count atoms in each formula
        def parse_formula(f):
            # Simplified parser: counts C, H, N, O
            counts = {}
            for element in ['C', 'H', 'N', 'O', 'S', 'P']:
                counts[element] = f.count(element)
            return counts
        
        counts1 = parse_formula(formula1)
        counts2 = parse_formula(formula2)
        
        # Jaccard similarity
        all_elements = set(counts1.keys()) | set(counts2.keys())
        if not all_elements:
            return 1.0
        
        intersection = sum(min(counts1.get(e, 0), counts2.get(e, 0)) for e in all_elements)
        union = sum(max(counts1.get(e, 0), counts2.get(e, 0)) for e in all_elements)
        
        if union == 0:
            return 1.0
        
        return intersection / union
    
    def _simple_fingerprint(self, mol_sim: Dict, mol_ref: Dict) -> float:
        """Fallback fingerprint without RDKit"""
        # Use atom counts and bond types
        n_atoms_sim = len(mol_sim.get('atoms', []))
        n_atoms_ref = len(mol_ref.get('atoms', []))
        n_bonds_sim = len(mol_sim.get('bonds', []))
        n_bonds_ref = len(mol_ref.get('bonds', []))
        
        atom_similarity = 1.0 - abs(n_atoms_sim - n_atoms_ref) / max(n_atoms_sim, n_atoms_ref, 1)
        bond_similarity = 1.0 - abs(n_bonds_sim - n_bonds_ref) / max(n_bonds_sim, n_bonds_ref, 1)
        
        return (atom_similarity + bond_similarity) / 2.0
    
    def _to_rdkit_mol(self, mol_dict: Dict):
        """Convert molecule dict to RDKit Mol (if needed)"""
        # Placeholder - full implementation would build RDKit mol from atoms/bonds
        # For now, assume SMILES is provided
        smiles = mol_dict.get('smiles', None)
        if smiles:
            from rdkit import Chem
            return Chem.MolFromSmiles(smiles)
        return None
    
    def _build_adjacency_matrix(self, mol_dict: Dict) -> np.ndarray:
        """Build adjacency matrix from bonds"""
        n_atoms = len(mol_dict.get('atoms', []))
        adj = np.zeros((n_atoms, n_atoms))
        
        for bond in mol_dict.get('bonds', []):
            i, j = bond[0], bond[1]
            adj[i, j] = 1
            adj[j, i] = 1
        
        return adj
    
    def _compute_graph_spectrum(self, adj: np.ndarray) -> np.ndarray:
        """Compute Laplacian eigenvalues"""
        from scipy import linalg
        
        # Degree matrix
        deg = np.diag(np.sum(adj, axis=1))
        
        # Laplacian
        L = deg - adj
        
        # Eigenvalues
        eigenvalues = linalg.eigvalsh(L)
        eigenvalues = np.sort(eigenvalues)
        
        return eigenvalues
    
    def _compare_spectra(self, eigs1: np.ndarray, eigs2: np.ndarray) -> float:
        """Compare two graph spectra"""
        # Pad to same length
        max_len = max(len(eigs1), len(eigs2))
        eigs1_pad = np.pad(eigs1, (0, max_len - len(eigs1)), constant_values=0)
        eigs2_pad = np.pad(eigs2, (0, max_len - len(eigs2)), constant_values=0)
        
        # Normalized difference
        diff = np.linalg.norm(eigs1_pad - eigs2_pad)
        norm = np.linalg.norm(eigs1_pad) + np.linalg.norm(eigs2_pad)
        
        if norm > 1e-10:
            similarity = 1.0 - diff / norm
        else:
            similarity = 1.0
        
        return max(0.0, similarity)
    
    def _compute_rmsd(self, pos1: np.ndarray, pos2: np.ndarray) -> float:
        """Compute RMSD between two sets of positions"""
        # Center both structures
        pos1_centered = pos1 - np.mean(pos1, axis=0)
        pos2_centered = pos2 - np.mean(pos2, axis=0)
        
        # Compute RMSD (without rotation alignment for simplicity)
        # Full version would use Kabsch algorithm
        rmsd = np.sqrt(np.mean(np.sum((pos1_centered - pos2_centered)**2, axis=1)))
        
        return rmsd


# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("=" * 70)
    print("MULTI-METRIC SIMILARITY - EXAMPLE")
    print("=" * 70)
    
    # Example molecules (simplified dicts)
    mol_sim = {
        'formula': 'C2H4O2',
        'atoms': ['C', 'C', 'O', 'O', 'H', 'H', 'H', 'H'],
        'bonds': [(0, 1), (0, 2), (1, 3), (0, 4), (0, 5), (1, 6), (1, 7)],
        'energy': -150.5,
        'smiles': 'C(C(=O)O)O'  # glycolic acid
    }
    
    mol_ref_same = {
        'formula': 'C2H4O2',
        'atoms': ['C', 'C', 'O', 'O', 'H', 'H', 'H', 'H'],
        'bonds': [(0, 1), (0, 2), (1, 3), (0, 4), (0, 5), (1, 6), (1, 7)],
        'energy': -148.2,
        'smiles': 'C(C(=O)O)O'
    }
    
    mol_ref_different = {
        'formula': 'C2H6O',
        'atoms': ['C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H'],
        'bonds': [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5), (1, 6), (1, 7), (2, 8)],
        'energy': -120.0,
        'smiles': 'CCO'  # ethanol
    }
    
    # Compute similarities
    similarity = MultiMetricSimilarity()
    
    print("\nTest 1: Same molecule (glycolic acid vs glycolic acid)")
    score_same = similarity.compute(mol_sim, mol_ref_same, include_geometric=False)
    print(f"  {score_same}")
    print(f"  Overall: {score_same.overall:.3f} (expect ~1.0)")
    
    print("\nTest 2: Different molecules (glycolic acid vs ethanol)")
    score_diff = similarity.compute(mol_sim, mol_ref_different, include_geometric=False)
    print(f"  {score_diff}")
    print(f"  Overall: {score_diff.overall:.3f} (expect <0.5)")
    
    print("\n" + "=" * 70)

