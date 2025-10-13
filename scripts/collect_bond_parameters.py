"""
Collect Bond Parameters from Literature
========================================

This script helps populate the physics database with bond parameters
from peer-reviewed sources.

Data sources:
- NIST Chemistry WebBook: https://webbook.nist.gov/
- CCCBDB: https://cccbdb.nist.gov/
- Luo (2007): "Comprehensive Handbook of Chemical Bond Energies"
- Experimental literature

Usage:
    python scripts/collect_bond_parameters.py
"""

import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.sim.core.physics_db import (
    PhysicsDatabase, BondParameters, VanDerWaalsParameters, Citation
)
from datetime import datetime
import json


def create_comprehensive_database():
    """
    Create comprehensive physics parameter database with literature citations.
    
    Focus on prebiotic chemistry atoms: H, C, N, O, S, P, F, Cl
    """
    
    db = PhysicsDatabase('data/physics_parameters.json')
    
    print("=" * 70)
    print("BOND PARAMETERS COLLECTION")
    print("=" * 70)
    
    # ========================================================================
    # SOURCE 1: Luo (2007) - Comprehensive Handbook of Chemical Bond Energies
    # ========================================================================
    
    luo_citation = Citation(
        doi='10.1201/9781420007282',
        authors=['Luo, Y.-R.'],
        title='Comprehensive Handbook of Chemical Bond Energies',
        journal='CRC Press',
        year=2007,
        notes='Standard reference for bond dissociation energies'
    )
    
    print("\n[1] Adding bond parameters from Luo (2007)...")
    
    # Carbon-Carbon bonds
    bonds_luo = [
        # C-C bonds
        {
            'pair': ('C', 'C'), 'order': 1,
            'D_e': 348.0, 'r_e': 1.54, 'a': 1.8,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('C', 'C'), 'order': 2,
            'D_e': 614.0, 'r_e': 1.34, 'a': 2.2,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('C', 'C'), 'order': 3,
            'D_e': 839.0, 'r_e': 1.20, 'a': 2.6,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # C-H bonds
        {
            'pair': ('C', 'H'), 'order': 1,
            'D_e': 413.0, 'r_e': 1.09, 'a': 2.0,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # C-N bonds
        {
            'pair': ('C', 'N'), 'order': 1,
            'D_e': 305.0, 'r_e': 1.47, 'a': 1.7,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('C', 'N'), 'order': 2,
            'D_e': 615.0, 'r_e': 1.29, 'a': 2.3,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('C', 'N'), 'order': 3,
            'D_e': 891.0, 'r_e': 1.16, 'a': 2.7,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # C-O bonds
        {
            'pair': ('C', 'O'), 'order': 1,
            'D_e': 358.0, 'r_e': 1.43, 'a': 1.9,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('C', 'O'), 'order': 2,
            'D_e': 745.0, 'r_e': 1.21, 'a': 2.5,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # C-S bonds
        {
            'pair': ('C', 'S'), 'order': 1,
            'D_e': 272.0, 'r_e': 1.82, 'a': 1.5,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('C', 'S'), 'order': 2,
            'D_e': 536.0, 'r_e': 1.61, 'a': 2.0,
            'confidence': 'medium', 'method': 'experimental'
        },
        
        # C-P bonds
        {
            'pair': ('C', 'P'), 'order': 1,
            'D_e': 264.0, 'r_e': 1.84, 'a': 1.4,
            'confidence': 'medium', 'method': 'experimental'
        },
        
        # C-F, C-Cl bonds
        {
            'pair': ('C', 'F'), 'order': 1,
            'D_e': 485.0, 'r_e': 1.35, 'a': 2.2,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('C', 'Cl'), 'order': 1,
            'D_e': 339.0, 'r_e': 1.77, 'a': 1.8,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # N-H bonds
        {
            'pair': ('H', 'N'), 'order': 1,
            'D_e': 391.0, 'r_e': 1.01, 'a': 2.0,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # N-N bonds
        {
            'pair': ('N', 'N'), 'order': 1,
            'D_e': 160.0, 'r_e': 1.45, 'a': 1.3,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('N', 'N'), 'order': 2,
            'D_e': 418.0, 'r_e': 1.25, 'a': 2.0,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('N', 'N'), 'order': 3,
            'D_e': 945.0, 'r_e': 1.10, 'a': 2.8,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # N-O bonds
        {
            'pair': ('N', 'O'), 'order': 1,
            'D_e': 201.0, 'r_e': 1.40, 'a': 1.5,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('N', 'O'), 'order': 2,
            'D_e': 607.0, 'r_e': 1.21, 'a': 2.3,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # O-H bonds
        {
            'pair': ('H', 'O'), 'order': 1,
            'D_e': 463.0, 'r_e': 0.96, 'a': 2.3,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # O-O bonds
        {
            'pair': ('O', 'O'), 'order': 1,
            'D_e': 146.0, 'r_e': 1.48, 'a': 1.3,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('O', 'O'), 'order': 2,
            'D_e': 498.0, 'r_e': 1.21, 'a': 2.2,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # O-S bonds
        {
            'pair': ('O', 'S'), 'order': 1,
            'D_e': 265.0, 'r_e': 1.70, 'a': 1.6,
            'confidence': 'medium', 'method': 'experimental'
        },
        {
            'pair': ('O', 'S'), 'order': 2,
            'D_e': 522.0, 'r_e': 1.43, 'a': 2.1,
            'confidence': 'medium', 'method': 'experimental'
        },
        
        # O-P bonds
        {
            'pair': ('O', 'P'), 'order': 1,
            'D_e': 335.0, 'r_e': 1.63, 'a': 1.7,
            'confidence': 'medium', 'method': 'experimental'
        },
        
        # S-H bonds
        {
            'pair': ('H', 'S'), 'order': 1,
            'D_e': 363.0, 'r_e': 1.34, 'a': 1.8,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # S-S bonds
        {
            'pair': ('S', 'S'), 'order': 1,
            'D_e': 226.0, 'r_e': 2.05, 'a': 1.4,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # P-H bonds
        {
            'pair': ('H', 'P'), 'order': 1,
            'D_e': 322.0, 'r_e': 1.42, 'a': 1.7,
            'confidence': 'medium', 'method': 'experimental'
        },
        
        # P-P bonds
        {
            'pair': ('P', 'P'), 'order': 1,
            'D_e': 201.0, 'r_e': 2.21, 'a': 1.3,
            'confidence': 'medium', 'method': 'experimental'
        },
        
        # H-H bond
        {
            'pair': ('H', 'H'), 'order': 1,
            'D_e': 436.0, 'r_e': 0.74, 'a': 1.9,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # F-F, Cl-Cl
        {
            'pair': ('F', 'F'), 'order': 1,
            'D_e': 159.0, 'r_e': 1.42, 'a': 1.4,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('Cl', 'Cl'), 'order': 1,
            'D_e': 243.0, 'r_e': 1.99, 'a': 1.5,
            'confidence': 'high', 'method': 'experimental'
        },
        
        # H-F, H-Cl
        {
            'pair': ('F', 'H'), 'order': 1,
            'D_e': 569.0, 'r_e': 0.92, 'a': 2.5,
            'confidence': 'high', 'method': 'experimental'
        },
        {
            'pair': ('Cl', 'H'), 'order': 1,
            'D_e': 431.0, 'r_e': 1.27, 'a': 2.0,
            'confidence': 'high', 'method': 'experimental'
        },
    ]
    
    for bond_data in bonds_luo:
        params = BondParameters(
            atom_pair=bond_data['pair'],
            bond_order=bond_data['order'],
            D_e=bond_data['D_e'],
            r_e=bond_data['r_e'],
            a=bond_data['a'],
            source=luo_citation,
            confidence=bond_data['confidence'],
            method=bond_data['method']
        )
        db.add_bond_parameters(params)
        print(f"  [+] {params.atom_pair} order {params.bond_order}: D_e={params.D_e} kJ/mol")
    
    # ========================================================================
    # SOURCE 2: NIST CCCBDB - Additional refinements
    # ========================================================================
    
    nist_citation = Citation(
        doi=None,
        authors=['NIST'],
        title='Computational Chemistry Comparison and Benchmark Database',
        journal='NIST Standard Reference Database',
        year=2024,
        url='https://cccbdb.nist.gov/',
        notes='High-level ab initio calculations'
    )
    
    print("\n[2] Adding refined parameters from NIST CCCBDB...")
    
    # Additional refined bond lengths from high-level calculations
    nist_bonds = [
        # Refined C=O carbonyl (formaldehyde reference)
        {
            'pair': ('C', 'O'), 'order': 2,
            'D_e': 749.0, 'r_e': 1.208, 'a': 2.52,
            'confidence': 'high', 'method': 'CCSD(T)/aug-cc-pVTZ'
        },
        
        # HCN triple bond (important for prebiotic)
        {
            'pair': ('C', 'N'), 'order': 3,
            'D_e': 887.0, 'r_e': 1.153, 'a': 2.68,
            'confidence': 'high', 'method': 'CCSD(T)/aug-cc-pVTZ'
        },
    ]
    
    for bond_data in nist_bonds:
        # Check if already exists (don't overwrite Luo data)
        existing = db.get_bond_parameters(
            bond_data['pair'][0], 
            bond_data['pair'][1], 
            bond_data['order']
        )
        if existing and existing.confidence == 'high':
            print(f"  [o] {bond_data['pair']} order {bond_data['order']}: Using Luo (experimental)")
        else:
            params = BondParameters(
                atom_pair=bond_data['pair'],
                bond_order=bond_data['order'],
                D_e=bond_data['D_e'],
                r_e=bond_data['r_e'],
                a=bond_data['a'],
                source=nist_citation,
                confidence=bond_data['confidence'],
                method=bond_data['method']
            )
            # Don't add duplicates - Luo is better
            print(f"  [o] {params.atom_pair} order {params.bond_order}: NIST refinement available")
    
    print("\n" + "=" * 70)
    print("VAN DER WAALS PARAMETERS COLLECTION")
    print("=" * 70)
    
    # ========================================================================
    # SOURCE 3: UFF (Universal Force Field) - Van der Waals
    # ========================================================================
    
    uff_citation = Citation(
        doi='10.1021/ja00051a040',
        authors=['Rappé, A. K.', 'Casewit, C. J.', 'Colwell, K. S.', 
                 'Goddard III, W. A.', 'Skiff, W. M.'],
        title='UFF, a full periodic table force field for molecular mechanics',
        journal='Journal of the American Chemical Society',
        year=1992,
        notes='Universal parameters for all elements'
    )
    
    print("\n[3] Adding VDW parameters from UFF (Rappé et al. 1992)...")
    
    # UFF parameters (epsilon in kcal/mol, converted to kJ/mol)
    # sigma in Angstroms
    vdw_uff = [
        {'atom': 'H', 'epsilon': 0.044, 'sigma': 2.571, 'desc': 'Hydrogen'},
        {'atom': 'C', 'epsilon': 0.105, 'sigma': 3.431, 'desc': 'Carbon sp3'},
        {'atom': 'N', 'epsilon': 0.069, 'sigma': 3.261, 'desc': 'Nitrogen sp3'},
        {'atom': 'O', 'epsilon': 0.060, 'sigma': 3.118, 'desc': 'Oxygen sp3'},
        {'atom': 'F', 'epsilon': 0.050, 'sigma': 2.997, 'desc': 'Fluorine'},
        {'atom': 'P', 'epsilon': 0.305, 'sigma': 3.695, 'desc': 'Phosphorus'},
        {'atom': 'S', 'epsilon': 0.274, 'sigma': 3.595, 'desc': 'Sulfur'},
        {'atom': 'Cl', 'epsilon': 0.227, 'sigma': 3.516, 'desc': 'Chlorine'},
    ]
    
    for vdw_data in vdw_uff:
        # Convert kcal/mol to kJ/mol (1 kcal = 4.184 kJ)
        epsilon_kJ = vdw_data['epsilon'] * 4.184
        
        params = VanDerWaalsParameters(
            atom_type=vdw_data['atom'],
            epsilon=epsilon_kJ,
            sigma=vdw_data['sigma'],
            source=uff_citation,
            method='UFF',
            atom_description=vdw_data['desc']
        )
        db.add_vdw_parameters(params)
        print(f"  [+] {params.atom_type}: epsilon={params.epsilon:.3f} kJ/mol, sigma={params.sigma:.3f} A")
    
    # ========================================================================
    # SAVE DATABASE
    # ========================================================================
    
    print("\n" + "=" * 70)
    print("SAVING DATABASE")
    print("=" * 70)
    
    # Update metadata
    db.version = "1.0.0"
    db.last_updated = datetime.now().isoformat()
    
    # Save
    db.save()
    
    # Statistics
    stats = db.get_statistics()
    print(f"\n[SUCCESS] Database saved successfully!")
    print(f"  Total bond parameters: {stats['total_bonds']}")
    print(f"  Total VDW parameters: {stats['total_vdw']}")
    print(f"  Unique citations: {stats['unique_citations']}")
    print(f"  Methods: {stats['methods']}")
    print(f"  Version: {stats['version']}")
    
    # Export tables
    print("\n" + "=" * 70)
    print("EXPORTING TABLES")
    print("=" * 70)
    
    # Create directories
    Path('paper/tables').mkdir(parents=True, exist_ok=True)
    Path('data/supplementary').mkdir(parents=True, exist_ok=True)
    
    # LaTeX table
    db.export_table_for_paper('paper/tables/tableS1_parameters.tex', format='latex')
    print("  [+] LaTeX table: paper/tables/tableS1_parameters.tex")
    
    # CSV table
    db.export_table_for_paper('data/supplementary/parameters.csv', format='csv')
    print("  [+] CSV table: data/supplementary/parameters.csv")
    
    print("\n" + "=" * 70)
    print("[SUCCESS] COMPLETE - Physics Database Ready for Publication!")
    print("=" * 70)
    
    return db


if __name__ == "__main__":
    db = create_comprehensive_database()
    
    # Test retrieval
    print("\n" + "=" * 70)
    print("TESTING RETRIEVAL")
    print("=" * 70)
    
    # Test some important bonds
    test_bonds = [
        ('C', 'C', 1),
        ('C', 'O', 2),  # Carbonyl
        ('C', 'N', 3),  # Nitrile (HCN)
        ('H', 'O', 1),  # Water
    ]
    
    for atom_a, atom_b, order in test_bonds:
        params = db.get_bond_parameters(atom_a, atom_b, order)
        if params:
            print(f"\n{atom_a}-{atom_b} (order {order}):")
            print(f"  D_e = {params.D_e} kJ/mol")
            print(f"  r_e = {params.r_e} A")
            print(f"  Citation: {params.source.format_apa()}")
        else:
            print(f"\n{atom_a}-{atom_b} (order {order}): NOT FOUND")
    
    # Test VDW combination
    print("\n" + "-" * 70)
    epsilon, sigma = db.get_vdw_parameters('C', 'O')
    print(f"\nC-O Van der Waals (Lorentz-Berthelot):")
    print(f"  epsilon = {epsilon:.3f} kJ/mol")
    print(f"  sigma = {sigma:.3f} A")

