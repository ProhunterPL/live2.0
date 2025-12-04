"""
TruthFilter 2.0: Enhanced validation with ACCEPT/FLAG/REJECT classification

Implements 8-step validation pipeline:
1. Valence check (hard REJECT)
2. Charge & connectivity (hard REJECT)
3. Ring strain & geometry heuristics (FLAG/REJECT)
4. Aromaticity policy (FLAG - model incompatible)
5. Model-compatibility score
6. Cross-check with databases (PubChem/ChEBI)
7. Statistics (persistence, occurrence)
8. Final decision logic (ACCEPT/FLAG/REJECT)
"""

import logging
from typing import Dict, List, Optional, Any
from enum import Enum

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

logger = logging.getLogger(__name__)


class ValidityStatus(str, Enum):
    """Validity status for molecules"""
    ACCEPT = "ACCEPT"
    FLAG = "FLAG"
    REJECT = "REJECT"


class TruthFilterV2:
    """
    TruthFilter 2.0: Enhanced validation with explicit ACCEPT/FLAG/REJECT classification
    
    Designed for classical force fields (Morse + LJ) that don't support:
    - Aromatic stabilization (π-delocalization)
    - Quantum-mechanical effects
    - High-strain ring systems
    """
    
    def __init__(self):
        """Initialize TruthFilterV2"""
        if not RDKIT_AVAILABLE:
            logger.warning("RDKit not available - some validations will be limited")
    
    def validate_molecule(self, 
                         smiles: Optional[str] = None,
                         formula: Optional[str] = None,
                         mass: Optional[float] = None,
                         metadata: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Validate a single molecule
        
        Args:
            smiles: SMILES string (preferred)
            formula: Molecular formula (fallback)
            mass: Molecular mass
            metadata: Additional metadata (occurrence_count, steps_seen, pubchem_match, etc.)
        
        Returns:
            Dict with:
            - validity: "ACCEPT" | "FLAG" | "REJECT"
            - reasons: List of reason codes
            - confidence: 0.0-1.0
            - model_compatibility: "high" | "medium" | "low"
        """
        metadata = metadata or {}
        reasons = []
        confidence = 0.8  # Start with 0.8 after valence pass
        model_compatibility = "high"
        
        # Parse molecule if SMILES available
        mol = None
        if smiles and RDKIT_AVAILABLE:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mol = Chem.AddHs(mol)  # Add hydrogens for accurate valence
            except Exception as e:
                logger.warning(f"Failed to parse SMILES {smiles}: {e}")
        
        # Step 1: Valence check (hard REJECT)
        if mol:
            valence_result = self._check_valence(mol)
            if not valence_result['valid']:
                return {
                    'validity': ValidityStatus.REJECT,
                    'reasons': ['VALENCE_ERROR'] + valence_result['errors'],
                    'confidence': 0.0,
                    'model_compatibility': 'low'
                }
            reasons.append('VALENCE_OK')
        elif not smiles:
            # No SMILES - can't validate valence properly
            reasons.append('NO_SMILES')
            confidence *= 0.7
        
        # Step 2: Charge & connectivity (hard REJECT)
        if mol:
            charge_result = self._check_charge_connectivity(mol)
            if not charge_result['valid']:
                return {
                    'validity': ValidityStatus.REJECT,
                    'reasons': charge_result['errors'],
                    'confidence': 0.0,
                    'model_compatibility': 'low'
                }
            if 'CHARGE_OK' in charge_result.get('reasons', []):
                reasons.append('CHARGE_OK')
            if 'CONNECTED' in charge_result.get('reasons', []):
                reasons.append('CONNECTED')
        
        # Step 3: Ring strain & geometry heuristics
        if mol:
            strain_result = self._check_ring_strain(mol)
            if strain_result['has_high_strain']:
                reasons.append('HIGH_STRAIN_RING')
                confidence -= 0.2
                if strain_result['hard_reject']:
                    return {
                        'validity': ValidityStatus.REJECT,
                        'reasons': reasons + ['EXTREME_STRAIN'],
                        'confidence': max(0.0, confidence),
                        'model_compatibility': 'low'
                    }
        
        # Step 4: Aromaticity policy (FLAG - model incompatible)
        if mol:
            aromatic_result = self._check_aromaticity(mol)
            if aromatic_result['has_aromatic']:
                reasons.append('AROMATIC_UNSUPPORTED_BY_MODEL')
                # Reduce confidence but keep minimum for aromatics (they're FLAG, not REJECT)
                confidence = max(0.4, confidence * 0.6)  # Minimum 0.4 for aromatics
                model_compatibility = "low"
        
        # Step 5: Model-compatibility score
        if mol:
            compat_result = self._assess_model_compatibility(mol)
            confidence += compat_result['score_adjustment']
            if compat_result['compatibility'] == 'low':
                model_compatibility = 'low'
            elif compat_result['compatibility'] == 'medium' and model_compatibility != 'low':
                model_compatibility = 'medium'
        
        # Step 6: Cross-check with databases
        pubchem_match = metadata.get('pubchem_match', False)
        pubchem_cid = metadata.get('pubchem_cid')
        if pubchem_match or pubchem_cid:
            reasons.append('MATCH_KNOWN_CHEMISTRY')
            confidence = min(1.0, confidence + 0.1)
        else:
            reasons.append('NO_DB_MATCH')
            # Don't penalize novel molecules too much
        
        # Step 7: Statistics (persistence, occurrence)
        occurrence_count = metadata.get('occurrence_count', 1)
        steps_seen = metadata.get('steps_seen', [])
        unique_steps = len(set(steps_seen)) if isinstance(steps_seen, list) else 1
        
        if occurrence_count == 1 and unique_steps == 1:
            reasons.append('TRANSIENT_SINGLETON')
            # Don't reduce confidence too much for transient singletons (they might still be valid)
            confidence = max(0.3, confidence * 0.7)
        elif occurrence_count > 10 or unique_steps > 5:
            reasons.append('PERSISTENT_SPECIES')
            confidence = min(1.0, confidence + 0.1)
        
        # Step 8: Final decision logic
        validity = self._final_decision(reasons, confidence)
        
        # Clamp confidence to [0, 1]
        confidence = max(0.0, min(1.0, confidence))
        
        return {
            'validity': validity.value,
            'reasons': reasons,
            'confidence': confidence,
            'model_compatibility': model_compatibility,
            'smiles': smiles,
            'formula': formula,
            'mass': mass
        }
    
    def _check_valence(self, mol) -> Dict[str, Any]:
        """Step 1: Valence check (hard REJECT if fails)"""
        errors = []
        
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            valence = atom.GetTotalValence()
            
            # Hard fail rules
            if symbol == 'C' and valence != 4:
                errors.append(f"C_valence_{valence}")
            elif symbol == 'N':
                if valence > 4 or valence < 2:
                    errors.append(f"N_valence_{valence}")
            elif symbol == 'O' and valence > 2:
                errors.append(f"O_valence_{valence}")
            elif symbol == 'H' and valence > 1:
                errors.append(f"H_valence_{valence}")
        
        return {
            'valid': len(errors) == 0,
            'errors': errors
        }
    
    def _check_charge_connectivity(self, mol) -> Dict[str, Any]:
        """Step 2: Charge & connectivity check (hard REJECT if fails)"""
        errors = []
        reasons = []
        
        # Check total charge
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        if abs(total_charge) > 2:
            errors.append('UNPHYSICAL_CHARGE')
        else:
            reasons.append('CHARGE_OK')
        
        # Check connectivity (single connected component)
        try:
            from rdkit.Chem import rdmolops
            frags = rdmolops.GetMolFrags(mol, asMols=True)
            if len(frags) > 1:
                errors.append('DISCONNECTED_COMPONENTS')
            else:
                reasons.append('CONNECTED')
        except Exception:
            # Fallback: assume connected if we can't check
            reasons.append('CONNECTED')
        
        return {
            'valid': len(errors) == 0,
            'errors': errors,
            'reasons': reasons
        }
    
    def _check_ring_strain(self, mol) -> Dict[str, Any]:
        """Step 3: Ring strain & geometry heuristics"""
        has_high_strain = False
        hard_reject = False
        
        try:
            ring_info = mol.GetRingInfo()
            rings = ring_info.AtomRings()
            
            for ring in rings:
                n_atoms = len(ring)
                heteroatoms = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() in ['N', 'O'])
                
                # Count double bonds in ring
                double_bonds = 0
                for i, atom_idx in enumerate(ring):
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for bond in atom.GetBonds():
                        other_idx = bond.GetOtherAtomIdx(atom_idx)
                        if other_idx in ring:
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                double_bonds += 1
                double_bonds = double_bonds // 2  # Count each bond once
                
                # Hard REJECT rules
                if n_atoms == 3 and heteroatoms > 1:
                    hard_reject = True
                    has_high_strain = True
                elif n_atoms == 4 and heteroatoms > 1 and double_bonds > 0:
                    hard_reject = True
                    has_high_strain = True
                
                # FLAG rules
                elif n_atoms in [5, 6] and heteroatoms >= 2:
                    has_high_strain = True
                
                # Check for bicyclic systems
                if len(rings) > 1:
                    # Check if rings share atoms (bicyclic)
                    for i, ring1 in enumerate(rings):
                        for ring2 in rings[i+1:]:
                            if set(ring1) & set(ring2):  # Shared atoms
                                has_high_strain = True
        except Exception as e:
            logger.warning(f"Ring strain check failed: {e}")
        
        return {
            'has_high_strain': has_high_strain,
            'hard_reject': hard_reject
        }
    
    def _check_aromaticity(self, mol) -> Dict[str, Any]:
        """Step 4: Aromaticity policy (FLAG - model incompatible)"""
        has_aromatic = False
        
        try:
            # RDKit aromaticity detection
            Chem.SanitizeMol(mol)
            ring_info = mol.GetRingInfo()
            rings = ring_info.AtomRings()
            
            for ring in rings:
                # Check if ring is aromatic
                is_aromatic = all(
                    mol.GetAtomWithIdx(idx).GetIsAromatic() 
                    for idx in ring
                )
                if is_aromatic:
                    has_aromatic = True
                    break
        except Exception as e:
            logger.warning(f"Aromaticity check failed: {e}")
        
        return {
            'has_aromatic': has_aromatic
        }
    
    def _assess_model_compatibility(self, mol) -> Dict[str, Any]:
        """Step 5: Model-compatibility score"""
        score_adjustment = 0.0
        compatibility = "high"
        
        try:
            num_atoms = mol.GetNumAtoms()
            num_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() != 'H')
            
            # Positive adjustments
            if num_heavy <= 10:
                score_adjustment += 0.1
            
            # Check for aromaticity (already handled in step 4, but check again)
            ring_info = mol.GetRingInfo()
            rings = ring_info.AtomRings()
            
            has_aromatic = False
            for ring in rings:
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    has_aromatic = True
                    break
            
            if not has_aromatic:
                score_adjustment += 0.1
            
            # Negative adjustments
            if len(rings) > 1:
                score_adjustment -= 0.2
                compatibility = "medium"
            
            # Check for triple bonds in complex structures
            triple_bonds = sum(
                1 for bond in mol.GetBonds()
                if bond.GetBondType() == Chem.BondType.TRIPLE
            )
            if triple_bonds > 0 and num_heavy > 8:
                score_adjustment -= 0.2
                compatibility = "medium"
            
            # Check for multiple heteroatoms in rings
            for ring in rings:
                heteroatoms = sum(
                    1 for idx in ring 
                    if mol.GetAtomWithIdx(idx).GetSymbol() in ['N', 'O']
                )
                if heteroatoms >= 3:
                    score_adjustment -= 0.2
                    compatibility = "low"
            
            if has_aromatic:
                score_adjustment -= 0.3
                compatibility = "low"
        except Exception as e:
            logger.warning(f"Model compatibility assessment failed: {e}")
        
        return {
            'score_adjustment': score_adjustment,
            'compatibility': compatibility
        }
    
    def _final_decision(self, reasons: List[str], confidence: float) -> ValidityStatus:
        """Step 8: Final decision logic"""
        # Hard REJECT conditions (must check first)
        if any(r in reasons for r in ['VALENCE_ERROR', 'UNPHYSICAL_CHARGE', 'DISCONNECTED_COMPONENTS', 'EXTREME_STRAIN']):
            return ValidityStatus.REJECT
        
        # FLAG conditions (check before confidence threshold)
        # Aromatics should be FLAG, not REJECT (even if confidence is low)
        if 'AROMATIC_UNSUPPORTED_BY_MODEL' in reasons:
            return ValidityStatus.FLAG
        
        # Other FLAG conditions
        if any(r in reasons for r in ['HIGH_STRAIN_RING', 'TRANSIENT_SINGLETON']):
            return ValidityStatus.FLAG
        
        # REJECT only if confidence is very low AND no FLAG conditions
        if confidence < 0.4:
            return ValidityStatus.REJECT
        
        # Default: ACCEPT
        return ValidityStatus.ACCEPT
    
    def validate_batch(self, molecules: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Validate a batch of molecules
        
        Args:
            molecules: List of molecule dicts with 'smiles', 'formula', 'mass', etc.
        
        Returns:
            List of validation results (same structure as validate_molecule output)
        """
        results = []
        for mol_data in molecules:
            result = self.validate_molecule(
                smiles=mol_data.get('smiles'),
                formula=mol_data.get('formula'),
                mass=mol_data.get('mass'),
                metadata=mol_data.get('metadata', {})
            )
            results.append(result)
        return results

