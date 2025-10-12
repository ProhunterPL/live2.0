"""
Parameter Validation Script
===========================

Validates physics_parameters.json against schema and checks for:
- Missing DOIs
- Duplicate entries
- Physical plausibility
- Completeness
"""

import json
import sys
from pathlib import Path
from typing import List, Dict, Tuple

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'backend'))

from sim.core.physics_db import PhysicsDatabase, BondParameters, VanDerWaalsParameters  # type: ignore


class ParameterValidator:
    """Validates physics parameters database"""
    
    def __init__(self, db_path: str):
        self.db = PhysicsDatabase(db_path)
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.info: List[str] = []
    
    def validate_all(self) -> bool:
        """Run all validation checks"""
        print("=" * 60)
        print("PHYSICS PARAMETERS VALIDATION")
        print("=" * 60)
        
        self.check_citations()
        self.check_physical_plausibility()
        self.check_completeness()
        self.check_duplicates()
        self.check_consistency()
        
        return self.print_report()
    
    def check_citations(self):
        """Check that all parameters have proper citations"""
        print("\n[1/5] Checking citations...")
        
        missing_doi = []
        incomplete_citations = []
        
        # Check bonds
        for (pair, order), params in self.db.bonds.items():
            key = f"{pair[0]}-{pair[1]} (order {order})"
            
            if not params.source.doi:
                missing_doi.append(key)
            
            if not params.source.authors or not params.source.year:
                incomplete_citations.append(key)
        
        # Check VDW
        for atom_type, params in self.db.vdw.items():
            if not params.source.doi:
                missing_doi.append(f"VDW: {atom_type}")
            
            if not params.source.authors or not params.source.year:
                incomplete_citations.append(f"VDW: {atom_type}")
        
        if missing_doi:
            self.warnings.append(f"Missing DOIs for {len(missing_doi)} parameters")
            for item in missing_doi[:5]:  # Show first 5
                self.warnings.append(f"  - {item}")
        
        if incomplete_citations:
            self.errors.append(f"Incomplete citations for {len(incomplete_citations)} parameters")
    
    def check_physical_plausibility(self):
        """Check that parameters are physically reasonable"""
        print("[2/5] Checking physical plausibility...")
        
        # Bond parameter ranges (reasonable values)
        D_e_range = (50, 1000)    # kJ/mol
        r_e_range = (0.5, 3.0)    # Angstroms
        a_range = (0.5, 5.0)      # 1/Angstrom
        
        # VDW parameter ranges
        epsilon_range = (0.01, 2.0)  # kJ/mol
        sigma_range = (2.0, 5.0)     # Angstroms
        
        for (pair, order), params in self.db.bonds.items():
            key = f"{pair[0]}-{pair[1]} (order {order})"
            
            if not (D_e_range[0] <= params.D_e <= D_e_range[1]):
                self.warnings.append(f"{key}: D_e = {params.D_e} kJ/mol outside typical range")
            
            if not (r_e_range[0] <= params.r_e <= r_e_range[1]):
                self.errors.append(f"{key}: r_e = {params.r_e} A outside reasonable range")
            
            if not (a_range[0] <= params.a <= a_range[1]):
                self.warnings.append(f"{key}: a = {params.a} 1/A outside typical range")
            
            # Check that higher bond orders have higher D_e
            single_bond_key = (pair, 1)
            if order > 1 and single_bond_key in self.db.bonds:
                single_params = self.db.bonds[single_bond_key]
                if params.D_e <= single_params.D_e:
                    self.warnings.append(
                        f"{key}: D_e ({params.D_e}) should be > single bond ({single_params.D_e})"
                    )
        
        for atom_type, params in self.db.vdw.items():
            if not (epsilon_range[0] <= params.epsilon <= epsilon_range[1]):
                self.warnings.append(
                    f"VDW {atom_type}: epsilon = {params.epsilon} kJ/mol outside typical range"
                )
            
            if not (sigma_range[0] <= params.sigma <= sigma_range[1]):
                self.warnings.append(
                    f"VDW {atom_type}: sigma = {params.sigma} A outside typical range"
                )
    
    def check_completeness(self):
        """Check for missing important parameters"""
        print("[3/5] Checking completeness...")
        
        # Important bonds for prebiotic chemistry
        important_bonds = [
            ('C', 'C'), ('C', 'H'), ('C', 'N'), ('C', 'O'),
            ('N', 'H'), ('O', 'H'), ('C', 'S'), ('N', 'N'),
            ('O', 'O'), ('S', 'S')
        ]
        
        # Important atoms
        important_atoms = ['H', 'C', 'N', 'O', 'S', 'P']
        
        missing_bonds = []
        for pair in important_bonds:
            normalized_pair = tuple(sorted(pair))
            if (normalized_pair, 1) not in self.db.bonds:
                missing_bonds.append(f"{pair[0]}-{pair[1]}")
        
        missing_vdw = []
        for atom in important_atoms:
            if atom not in self.db.vdw:
                missing_vdw.append(atom)
        
        if missing_bonds:
            self.info.append(f"Missing {len(missing_bonds)} important bonds: {', '.join(missing_bonds[:5])}")
        
        if missing_vdw:
            self.warnings.append(f"Missing VDW parameters for: {', '.join(missing_vdw)}")
    
    def check_duplicates(self):
        """Check for duplicate entries"""
        print("[4/5] Checking for duplicates...")
        
        # This should never happen with current implementation, but check anyway
        bond_keys = list(self.db.bonds.keys())
        if len(bond_keys) != len(set(bond_keys)):
            self.errors.append("Duplicate bond entries detected!")
        
        vdw_keys = list(self.db.vdw.keys())
        if len(vdw_keys) != len(set(vdw_keys)):
            self.errors.append("Duplicate VDW entries detected!")
    
    def check_consistency(self):
        """Check internal consistency"""
        print("[5/5] Checking consistency...")
        
        for (pair, order), params in self.db.bonds.items():
            # Check k_spring calculation
            expected_k = 2 * params.D_e * (params.a ** 2)
            if params.k_spring:
                diff = abs(params.k_spring - expected_k) / expected_k
                if diff > 0.01:  # 1% tolerance
                    self.warnings.append(
                        f"{pair[0]}-{pair[1]} (order {order}): "
                        f"k_spring inconsistent with Morse parameters"
                    )
            
            # Check r_0 vs r_e
            if params.r_0 and abs(params.r_0 - params.r_e) > 0.01:
                self.info.append(
                    f"{pair[0]}-{pair[1]} (order {order}): "
                    f"r_0 ({params.r_0}) != r_e ({params.r_e})"
                )
    
    def print_report(self) -> bool:
        """Print validation report"""
        print("\n" + "=" * 60)
        print("VALIDATION REPORT")
        print("=" * 60)
        
        stats = self.db.get_statistics()
        print(f"\nDatabase Statistics:")
        print(f"  Total bonds: {stats['total_bonds']}")
        print(f"  Total VDW: {stats['total_vdw']}")
        print(f"  Unique citations: {stats['unique_citations']}")
        print(f"  Methods: {stats['methods']}")
        
        print(f"\nValidation Results:")
        print(f"  Errors: {len(self.errors)}")
        print(f"  Warnings: {len(self.warnings)}")
        print(f"  Info: {len(self.info)}")
        
        if self.errors:
            print("\n[X] ERRORS:")
            for error in self.errors:
                print(f"  - {error}")
        
        if self.warnings:
            print("\n[!] WARNINGS:")
            for warning in self.warnings:
                print(f"  - {warning}")
        
        if self.info:
            print("\n[i] INFO:")
            for info_msg in self.info:
                print(f"  - {info_msg}")
        
        print("\n" + "=" * 60)
        
        if self.errors:
            print("[X] VALIDATION FAILED")
            return False
        elif self.warnings:
            print("[!] VALIDATION PASSED WITH WARNINGS")
            return True
        else:
            print("[OK] VALIDATION PASSED")
            return True


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Validate physics parameters database")
    parser.add_argument(
        '--db',
        default='data/physics_parameters.json',
        help='Path to physics_parameters.json'
    )
    parser.add_argument(
        '--strict',
        action='store_true',
        help='Treat warnings as errors'
    )
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.db).exists():
        print(f"[X] Database file not found: {args.db}")
        print("   Use the example: data/physics_parameters_example.json")
        print("   Or create a new database using physics_db.py")
        sys.exit(1)
    
    # Validate
    validator = ParameterValidator(args.db)
    success = validator.validate_all()
    
    if args.strict and validator.warnings:
        print("\n[!] Strict mode: treating warnings as errors")
        success = False
    
    sys.exit(0 if success else 1)

