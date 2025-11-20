"""
Enhanced Molecule Extractor with Atom Types
============================================

Extracts molecules WITH atom type information for PubChem matching.

Author: Live 2.0 Team
Date: 2025-11-20
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional
from collections import Counter
import numpy as np

logger = logging.getLogger(__name__)

# Periodic table mapping
ELEMENT_SYMBOLS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 
    16: 'S', 17: 'Cl', 18: 'Ar'
}

def infer_atom_type_from_mass(mass: float) -> str:
    """
    Infer atom type from mass
    
    Approximate masses:
    H: 1.0, C: 12.0, N: 14.0, O: 16.0, S: 32.0
    """
    if 0.5 < mass < 1.5:
        return 'H'
    elif 11 < mass < 13:
        return 'C'
    elif 13 < mass < 15:
        return 'N'
    elif 15 < mass < 17:
        return 'O'
    elif 31 < mass < 33:
        return 'S'
    elif 14 < mass < 16:
        return 'N'  # Could be N or O, default to N
    else:
        return 'X'  # Unknown

def calculate_molecular_formula(atoms: List[str]) -> str:
    """Calculate molecular formula from atom list"""
    atom_counts = Counter(atoms)
    
    # Order: C, H, N, O, S, others
    formula_order = ['C', 'H', 'N', 'O', 'S']
    formula = ""
    
    for element in formula_order:
        if element in atom_counts:
            count = atom_counts[element]
            if count == 1:
                formula += element
            else:
                formula += f"{element}{count}"
            del atom_counts[element]
    
    # Add remaining elements
    for element in sorted(atom_counts.keys()):
        count = atom_counts[element]
        if count == 1:
            formula += element
        else:
            formula += f"{element}{count}"
    
    return formula or "X"

class EnhancedMoleculeExtractor:
    """
    Enhanced extractor that includes atom types
    
    Usage:
        extractor = EnhancedMoleculeExtractor(results_dir)
        molecules = extractor.extract_molecules_with_types()
    """
    
    def __init__(self, results_dir: Path):
        self.results_dir = Path(results_dir)
        logger.info(f"EnhancedMoleculeExtractor initialized for: {self.results_dir}")
    
    def extract_molecules_with_types(self) -> List[Dict]:
        """
        Extract molecules WITH atom type information
        
        Returns:
            List of molecules with structure:
            {
                'formula': 'H2O',           # Real formula!
                'formula_hash': '...',      # Original hash
                'atoms': ['O', 'H', 'H'],   # Atom types!
                'bonds': [[0, 1], [0, 2]],  # Bond pairs
                'positions': [...],          # 3D positions (if available)
                'size': 3,
                'count': 10
            }
        """
        logger.info("Extracting molecules with atom types...")
        
        # Load results
        results_file = self.results_dir / "results.json"
        if not results_file.exists():
            logger.error(f"Results file not found: {results_file}")
            return []
        
        with open(results_file, 'r') as f:
            results = json.load(f)
        
        molecules_raw = results.get('molecules_detected', [])
        if not molecules_raw:
            logger.warning("No molecules found in results")
            return []
        
        # Try to get snapshot with full data
        snapshot_dir = self.results_dir / "snapshots"
        snapshot_data = None
        
        if snapshot_dir.exists():
            # Get last snapshot
            snapshots = sorted(snapshot_dir.glob("step_*.json"))
            if snapshots:
                with open(snapshots[-1], 'r') as f:
                    snapshot_data = json.load(f)
                logger.info(f"Loaded snapshot: {snapshots[-1].name}")
        
        # Extract molecules with types
        molecules_typed = []
        
        for mol in molecules_raw:
            try:
                mol_typed = self._extract_molecule_types(mol, snapshot_data)
                if mol_typed:
                    molecules_typed.append(mol_typed)
            except Exception as e:
                logger.warning(f"Failed to extract types for molecule: {e}")
        
        logger.info(f"Extracted {len(molecules_typed)} molecules with atom types")
        return molecules_typed
    
    def _extract_molecule_types(self, mol: Dict, snapshot: Optional[Dict]) -> Optional[Dict]:
        """Extract atom types for a single molecule"""
        
        cluster = mol.get('cluster', [])
        bonds = mol.get('bonds', [])
        size = mol.get('size', len(cluster))
        count = mol.get('count', 1)
        formula_hash = mol.get('formula', 'unknown')
        
        if not cluster:
            return None
        
        # Try to get atom types from snapshot
        atoms = []
        positions = []
        
        if snapshot and 'attributes' in snapshot:
            attributes = snapshot['attributes']
            snapshot_positions = snapshot.get('positions', [])
            
            for particle_idx in cluster:
                if particle_idx < len(attributes):
                    # attributes = [mass, charge_x, charge_y, charge_z]
                    mass = attributes[particle_idx][0]
                    atom_type = infer_atom_type_from_mass(mass)
                    atoms.append(atom_type)
                    
                    if particle_idx < len(snapshot_positions):
                        positions.append(snapshot_positions[particle_idx])
        
        # If couldn't get from snapshot, fallback to generic
        if not atoms:
            # Create placeholder atoms
            atoms = ['X'] * size
            logger.debug(f"Could not determine atom types, using placeholders")
        
        # Calculate real formula
        formula = calculate_molecular_formula(atoms)
        
        # Convert bonds to atom indices (0-indexed)
        bond_pairs = []
        cluster_map = {pid: i for i, pid in enumerate(cluster)}
        
        for bond in bonds:
            if len(bond) >= 2:
                i, j = bond[0], bond[1]
                if i in cluster_map and j in cluster_map:
                    bond_pairs.append([cluster_map[i], cluster_map[j]])
        
        return {
            'formula': formula,
            'formula_hash': formula_hash,
            'atoms': atoms,
            'bonds': bond_pairs,
            'positions': positions if positions else None,
            'size': size,
            'count': count
        }
    
    def save_for_matching(self, molecules: List[Dict], output_file: Path):
        """Save molecules in format ready for PubChem matching"""
        
        # Group by formula
        formula_groups = {}
        for mol in molecules:
            formula = mol['formula']
            if formula not in formula_groups:
                formula_groups[formula] = []
            formula_groups[formula].append(mol)
        
        # Create summary
        summary = []
        for formula, mols in sorted(formula_groups.items(), key=lambda x: -sum(m['count'] for m in x[1])):
            total_count = sum(m['count'] for m in mols)
            # Take first instance as representative
            representative = mols[0].copy()
            representative['count'] = total_count
            representative['instances'] = len(mols)
            summary.append(representative)
        
        output_data = {
            'total_molecules': len(molecules),
            'unique_formulas': len(formula_groups),
            'molecules': summary
        }
        
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        logger.info(f"Saved {len(summary)} unique molecules to: {output_file}")
        return output_data


def extract_run_with_types(run_dir: Path) -> List[Dict]:
    """
    Convenience function to extract molecules with types from a run
    
    Args:
        run_dir: Path to run directory (e.g., results/.../run_1/)
    
    Returns:
        List of molecules with atom types
    """
    extractor = EnhancedMoleculeExtractor(run_dir)
    molecules = extractor.extract_molecules_with_types()
    
    # Save for matching
    output_file = run_dir / "molecules_with_types.json"
    extractor.save_for_matching(molecules, output_file)
    
    return molecules


if __name__ == "__main__":
    # Test extraction
    import sys
    
    if len(sys.argv) > 1:
        run_dir = Path(sys.argv[1])
        print(f"Extracting molecules from: {run_dir}")
        molecules = extract_run_with_types(run_dir)
        print(f"Extracted {len(molecules)} molecules with atom types")
    else:
        print("Usage: python molecule_extractor_enhanced.py <run_directory>")

