"""
Molecule Filter: Filtrowanie molekuł pod kątem wiarygodności chemicznej
"""
import logging
from typing import List, Dict, Any, Tuple
from backend.validation.types import FilterResult, FilterStatus

logger = logging.getLogger(__name__)


def is_real_molecule(mol: Dict[str, Any]) -> bool:
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


def check_valence(mol: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Check for valence violations
    
    Returns:
        (is_valid, list_of_violations)
    """
    atoms = mol.get('atoms', [])
    bonds = mol.get('bonds', [])
    
    if not atoms or not bonds:
        return True, []  # Can't validate
    
    violations = []
    
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
            violations.append(f"{element} atom {i} has {bond_count[i]} bonds (max {max_val})")
    
    return len(violations) == 0, violations


def check_charge_balance(mol: Dict[str, Any]) -> Tuple[bool, float]:
    """
    Check for charge balance
    
    Returns:
        (is_valid, total_charge)
    """
    atoms = mol.get('atoms', [])
    
    total_charge = 0.0
    for atom in atoms:
        if isinstance(atom, dict):
            charge = atom.get('charge', 0.0)
            if isinstance(charge, (list, tuple)):
                # Charge as vector - compute magnitude
                charge = sum(c**2 for c in charge)**0.5
            total_charge += charge
    
    # Allow ±2 total charge
    is_valid = abs(total_charge) <= 2.0
    
    return is_valid, total_charge


def check_bond_orders(mol: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Check for unusual bond orders
    
    Returns:
        (is_valid, list_of_issues)
    """
    bonds = mol.get('bonds', [])
    issues = []
    
    for bond in bonds:
        if len(bond) > 2:
            bond_order = bond[2]
            if bond_order < 1 or bond_order > 3:
                issues.append(f"Bond {bond[:2]} has invalid order: {bond_order}")
    
    return len(issues) == 0, issues


class MoleculeFilter:
    """
    Filter molecules based on chemical plausibility
    """
    
    def __init__(self, strict_valence: bool = True, strict_charge: bool = True):
        """
        Initialize molecule filter
        
        Args:
            strict_valence: Reject molecules with valence violations
            strict_charge: Reject molecules with charge imbalance > 2
        """
        self.strict_valence = strict_valence
        self.strict_charge = strict_charge
    
    def filter(self, molecules: List[Dict[str, Any]]) -> tuple[FilterResult, List[Dict[str, Any]]]:
        """
        Filter molecules
        
        Args:
            molecules: List of molecule dicts
            
        Returns:
            Tuple of (FilterResult, filtered_molecules)
        """
        warnings = []
        errors = []
        
        real_molecules = []
        clusters_removed = []
        
        valence_violations = 0
        charge_issues = 0
        bond_order_issues = 0
        
        for mol in molecules:
            # Check if real molecule (not cluster)
            if not is_real_molecule(mol):
                clusters_removed.append(mol)
                continue
            
            # Check valence
            valence_ok, valence_viols = check_valence(mol)
            if not valence_ok:
                valence_violations += 1
                if self.strict_valence:
                    errors.append(f"Molecule {mol.get('formula', 'unknown')}: valence violations")
                    clusters_removed.append(mol)
                    continue
                else:
                    warnings.append(f"Molecule {mol.get('formula', 'unknown')}: {len(valence_viols)} valence violations")
            
            # Check charge
            charge_ok, total_charge = check_charge_balance(mol)
            if not charge_ok:
                charge_issues += 1
                if self.strict_charge:
                    errors.append(f"Molecule {mol.get('formula', 'unknown')}: charge imbalance {total_charge:.2f}")
                    clusters_removed.append(mol)
                    continue
                else:
                    warnings.append(f"Molecule {mol.get('formula', 'unknown')}: charge {total_charge:.2f}")
            
            # Check bond orders
            bond_ok, bond_issues = check_bond_orders(mol)
            if not bond_ok:
                bond_order_issues += 1
                warnings.append(f"Molecule {mol.get('formula', 'unknown')}: {len(bond_issues)} bond order issues")
            
            # Passed all checks
            real_molecules.append(mol)
        
        # Determine status
        if len(errors) > 0:
            status = FilterStatus.FAIL
        elif len(warnings) > 0 or valence_violations > 0 or charge_issues > 0:
            status = FilterStatus.WARNING
        else:
            status = FilterStatus.PASS
        
        # Calculate statistics
        total_original = len(molecules)
        total_filtered = len(real_molecules)
        retention_rate = total_filtered / total_original if total_original > 0 else 0.0
        
        details = {
            'total_original': total_original,
            'total_filtered': total_filtered,
            'clusters_removed': len(clusters_removed),
            'retention_rate': retention_rate,
            'valence_violations': valence_violations,
            'charge_issues': charge_issues,
            'bond_order_issues': bond_order_issues
        }
        
        return FilterResult(
            name="molecule_filter",
            status=status,
            details=details,
            warnings=warnings,
            errors=errors
        ), real_molecules
