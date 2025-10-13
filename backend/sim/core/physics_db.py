"""
Physics Parameters Database
===========================

Database of physical parameters (bond energies, VDW parameters, etc.) 
with citations from scientific literature.

This ensures all simulation parameters are traceable to peer-reviewed sources
rather than arbitrary values.

References:
-----------
- NIST Chemistry WebBook: https://webbook.nist.gov/
- CCCBDB: https://cccbdb.nist.gov/
- Luo (2007): "Comprehensive Handbook of Chemical Bond Energies"
  DOI: 10.1201/9781420007282
- Rappé et al. (1992): "UFF force field"
  DOI: 10.1021/ja00051a040
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Any
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class Citation:
    """Scientific citation for parameter source"""
    doi: Optional[str] = None
    authors: List[str] = field(default_factory=list)
    title: str = ""
    journal: str = ""
    year: int = 0
    url: Optional[str] = None
    notes: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            'doi': self.doi,
            'authors': self.authors,
            'title': self.title,
            'journal': self.journal,
            'year': self.year,
            'url': self.url,
            'notes': self.notes
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Citation':
        """Create from dictionary"""
        return cls(**data)
    
    def format_apa(self) -> str:
        """Format citation in APA style"""
        authors_str = ', '.join(self.authors[:3])
        if len(self.authors) > 3:
            authors_str += ', et al.'
        
        citation = f"{authors_str} ({self.year}). {self.title}. {self.journal}."
        if self.doi:
            citation += f" https://doi.org/{self.doi}"
        
        return citation


@dataclass
class BondParameters:
    """
    Bond parameters for specific atom pair and bond order.
    
    Morse potential: V(r) = D_e * (1 - exp(-a*(r - r_e)))^2
    Harmonic (spring): V(r) = 0.5 * k * (r - r_0)^2
    """
    atom_pair: Tuple[str, str]  # e.g., ('C', 'C')
    bond_order: int  # 1, 2, or 3 (single, double, triple)
    
    # Morse potential parameters
    D_e: float  # Dissociation energy (kJ/mol)
    r_e: float  # Equilibrium bond length (Angstroms)
    a: float    # Width parameter (1/Angstrom)
    
    # Harmonic approximation (optional, for small displacements)
    k_spring: Optional[float] = None  # Force constant (kJ/mol/Å²)
    r_0: Optional[float] = None       # Equilibrium length (Å)
    
    # Metadata
    source: Citation = field(default_factory=Citation)
    confidence: str = 'medium'  # 'high', 'medium', 'low'
    method: str = 'unknown'     # 'experimental', 'DFT', 'MP2', 'CCSD(T)', 'fitted'
    temperature: float = 298.15  # Temperature (K) for experimental values
    
    def __post_init__(self):
        """Validate and normalize data"""
        # Normalize atom pair (alphabetical order for consistency)
        self.atom_pair = tuple(sorted(self.atom_pair))
        
        # Validate bond order
        if self.bond_order not in [1, 2, 3]:
            raise ValueError(f"Bond order must be 1, 2, or 3, got {self.bond_order}")
        
        # Calculate k_spring from Morse if not provided
        if self.k_spring is None and self.D_e > 0 and self.a > 0:
            # k = 2 * D_e * a^2 (at equilibrium)
            self.k_spring = 2 * self.D_e * (self.a ** 2)
        
        # Set r_0 = r_e if not provided
        if self.r_0 is None:
            self.r_0 = self.r_e
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            'atom_pair': list(self.atom_pair),
            'bond_order': self.bond_order,
            'D_e': self.D_e,
            'r_e': self.r_e,
            'a': self.a,
            'k_spring': self.k_spring,
            'r_0': self.r_0,
            'source': self.source.to_dict(),
            'confidence': self.confidence,
            'method': self.method,
            'temperature': self.temperature
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'BondParameters':
        """Create from dictionary"""
        data = data.copy()
        data['atom_pair'] = tuple(data['atom_pair'])
        data['source'] = Citation.from_dict(data['source'])
        return cls(**data)


@dataclass
class VanDerWaalsParameters:
    """
    Van der Waals parameters for atom type.
    
    Lennard-Jones potential: V(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
    """
    atom_type: str  # e.g., 'C', 'N', 'O'
    
    # Lennard-Jones parameters
    epsilon: float  # Well depth (kJ/mol)
    sigma: float    # Zero-crossing distance (Angstroms)
    
    # Alternative representations
    r_min: Optional[float] = None  # Distance at minimum (Å), r_min = 2^(1/6) * sigma
    
    # Metadata
    source: Citation = field(default_factory=Citation)
    method: str = 'UFF'  # 'UFF', 'OPLS', 'AMBER', 'CHARMM', etc.
    atom_description: str = ''  # e.g., 'sp3 carbon', 'carbonyl oxygen'
    
    def __post_init__(self):
        """Validate and compute derived parameters"""
        if self.r_min is None and self.sigma > 0:
            # r_min = 2^(1/6) * sigma ≈ 1.122 * sigma
            self.r_min = 2.0 ** (1.0/6.0) * self.sigma
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'atom_type': self.atom_type,
            'epsilon': self.epsilon,
            'sigma': self.sigma,
            'r_min': self.r_min,
            'source': self.source.to_dict(),
            'method': self.method,
            'atom_description': self.atom_description
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'VanDerWaalsParameters':
        """Create from dictionary"""
        data = data.copy()
        data['source'] = Citation.from_dict(data['source'])
        return cls(**data)


class PhysicsDatabase:
    """
    Central database for all physics parameters with literature citations.
    
    Usage:
        db = PhysicsDatabase('data/physics_parameters.json')
        params = db.get_bond_parameters('C', 'C', order=1)
        print(params.source.format_apa())
    """
    
    def __init__(self, db_path: str = 'data/physics_parameters.json'):
        self.db_path = Path(db_path)
        
        # Storage
        self.bonds: Dict[Tuple[Tuple[str, str], int], BondParameters] = {}
        self.vdw: Dict[str, VanDerWaalsParameters] = {}
        
        # Metadata
        self.version = "1.0.0"
        self.last_updated = None
        self.total_citations = 0
        
        # Load if exists
        if self.db_path.exists():
            self.load()
        else:
            logger.warning(f"Database file not found: {db_path}")
            logger.info("Starting with empty database. Use add_* methods to populate.")
    
    def add_bond_parameters(self, params: BondParameters):
        """Add bond parameters to database"""
        key = (params.atom_pair, params.bond_order)
        self.bonds[key] = params
        logger.info(f"Added bond parameters: {params.atom_pair} (order {params.bond_order})")
    
    def add_vdw_parameters(self, params: VanDerWaalsParameters):
        """Add VDW parameters to database"""
        self.vdw[params.atom_type] = params
        logger.info(f"Added VDW parameters: {params.atom_type}")
    
    def get_bond_parameters(self, atom_a: str, atom_b: str, order: int = 1) -> Optional[BondParameters]:
        """
        Get bond parameters for atom pair with citation.
        
        Args:
            atom_a: First atom type
            atom_b: Second atom type
            order: Bond order (1, 2, or 3)
        
        Returns:
            BondParameters with citation, or None if not found
        """
        # Normalize to alphabetical order
        pair = tuple(sorted([atom_a, atom_b]))
        key = (pair, order)
        
        params = self.bonds.get(key)
        
        if params is None:
            logger.warning(f"No bond parameters found for {pair} (order {order})")
        
        return params
    
    def get_vdw_parameters(self, atom_a: str, atom_b: str) -> Tuple[float, float]:
        """
        Get VDW parameters for atom pair using Lorentz-Berthelot combination rules.
        
        Combination rules:
        - epsilon_ij = sqrt(epsilon_i * epsilon_j)
        - sigma_ij = (sigma_i + sigma_j) / 2
        
        Returns:
            (epsilon, sigma) tuple
        """
        params_a = self.vdw.get(atom_a)
        params_b = self.vdw.get(atom_b)
        
        if params_a is None or params_b is None:
            logger.warning(f"Missing VDW parameters for {atom_a}-{atom_b} pair")
            # Return generic values as fallback
            return (0.5, 3.4)  # Typical C-C values
        
        # Lorentz-Berthelot rules
        epsilon = (params_a.epsilon * params_b.epsilon) ** 0.5
        sigma = (params_a.sigma + params_b.sigma) / 2.0
        
        return (epsilon, sigma)
    
    def get_single_vdw(self, atom_type: str) -> Optional[VanDerWaalsParameters]:
        """Get VDW parameters for single atom type"""
        return self.vdw.get(atom_type)
    
    def save(self, path: Optional[Path] = None):
        """Save database to JSON file"""
        if path is None:
            path = self.db_path
        
        # Ensure directory exists
        path.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert to dict
        data = {
            'version': self.version,
            'last_updated': self.last_updated,
            'bonds': {
                f"{pair[0]}_{pair[1]}_order{order}": params.to_dict()
                for (pair, order), params in self.bonds.items()
            },
            'vdw': {
                atom_type: params.to_dict()
                for atom_type, params in self.vdw.items()
            }
        }
        
        # Write JSON
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)
        
        logger.info(f"Saved database to {path}")
        logger.info(f"  Bonds: {len(self.bonds)}")
        logger.info(f"  VDW: {len(self.vdw)}")
    
    def load(self, path: Optional[Path] = None):
        """Load database from JSON file"""
        if path is None:
            path = self.db_path
        
        if not path.exists():
            logger.error(f"Database file not found: {path}")
            return
        
        with open(path, 'r') as f:
            data = json.load(f)
        
        self.version = data.get('version', '1.0.0')
        self.last_updated = data.get('last_updated')
        
        # Load bonds
        self.bonds = {}
        for key_str, params_dict in data.get('bonds', {}).items():
            params = BondParameters.from_dict(params_dict)
            key = (params.atom_pair, params.bond_order)
            self.bonds[key] = params
        
        # Load VDW
        self.vdw = {}
        for atom_type, params_dict in data.get('vdw', {}).items():
            self.vdw[atom_type] = VanDerWaalsParameters.from_dict(params_dict)
        
        logger.info(f"Loaded database from {path}")
        logger.info(f"  Bonds: {len(self.bonds)}")
        logger.info(f"  VDW: {len(self.vdw)}")
    
    def export_table_for_paper(self, output_path: str, format: str = 'latex'):
        """
        Generate publication-ready table of parameters.
        
        Args:
            output_path: Where to save table
            format: 'latex' or 'csv'
        """
        if format == 'latex':
            self._export_latex_table(output_path)
        elif format == 'csv':
            self._export_csv_table(output_path)
        else:
            raise ValueError(f"Unknown format: {format}")
    
    def _export_latex_table(self, output_path: str):
        """Export as LaTeX table for paper"""
        lines = [
            r"\begin{table}[h]",
            r"\centering",
            r"\caption{Physical Parameters from Literature}",
            r"\label{tab:parameters}",
            r"\begin{tabular}{llccccl}",
            r"\hline",
            r"Atoms & Order & $D_e$ (kJ/mol) & $r_e$ (\AA) & $a$ (\AA$^{-1}$) & Method & Reference \\",
            r"\hline"
        ]
        
        # Bond parameters
        for (pair, order), params in sorted(self.bonds.items()):
            atom_str = f"{pair[0]}--{pair[1]}"
            ref_str = f"\\cite{{{params.source.doi or 'unknown'}}}"
            
            line = (f"{atom_str} & {order} & {params.D_e:.1f} & "
                   f"{params.r_e:.3f} & {params.a:.3f} & "
                   f"{params.method} & {ref_str} \\\\")
            lines.append(line)
        
        lines.extend([
            r"\hline",
            r"\end{tabular}",
            r"\end{table}"
        ])
        
        with open(output_path, 'w') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Exported LaTeX table to {output_path}")
    
    def _export_csv_table(self, output_path: str):
        """Export as CSV for supplementary materials"""
        import csv
        
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Header
            writer.writerow([
                'Atom 1', 'Atom 2', 'Order', 'D_e (kJ/mol)', 'r_e (Å)', 
                'a (1/Å)', 'Method', 'DOI', 'Citation'
            ])
            
            # Data
            for (pair, order), params in sorted(self.bonds.items()):
                writer.writerow([
                    pair[0], pair[1], order,
                    params.D_e, params.r_e, params.a,
                    params.method,
                    params.source.doi or '',
                    params.source.format_apa()
                ])
        
        logger.info(f"Exported CSV table to {output_path}")
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics"""
        unique_dois = set()
        for params in self.bonds.values():
            if params.source.doi:
                unique_dois.add(params.source.doi)
        for params in self.vdw.values():
            if params.source.doi:
                unique_dois.add(params.source.doi)
        
        methods = {}
        for params in self.bonds.values():
            methods[params.method] = methods.get(params.method, 0) + 1
        
        return {
            'total_bonds': len(self.bonds),
            'total_vdw': len(self.vdw),
            'unique_citations': len(unique_dois),
            'methods': methods,
            'version': self.version
        }


# Example usage and testing
if __name__ == "__main__":
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    
    # Create example database
    db = PhysicsDatabase('data/physics_parameters_example.json')
    
    # Example: C-C single bond (from Luo 2007)
    cc_single = BondParameters(
        atom_pair=('C', 'C'),
        bond_order=1,
        D_e=348.0,  # kJ/mol
        r_e=1.54,   # Angstroms
        a=1.8,      # 1/Angstrom (typical)
        source=Citation(
            doi='10.1201/9781420007282',
            authors=['Luo, Y.-R.'],
            title='Comprehensive Handbook of Chemical Bond Energies',
            journal='CRC Press',
            year=2007
        ),
        confidence='high',
        method='experimental'
    )
    db.add_bond_parameters(cc_single)
    
    # Example: Carbon VDW (from UFF)
    c_vdw = VanDerWaalsParameters(
        atom_type='C',
        epsilon=0.105,  # kcal/mol = 0.439 kJ/mol
        sigma=3.431,    # Angstroms
        source=Citation(
            doi='10.1021/ja00051a040',
            authors=['Rappé, A. K.', 'Casewit, C. J.', 'Colwell, K. S.', 'Goddard III, W. A.', 'Skiff, W. M.'],
            title='UFF, a full periodic table force field',
            journal='Journal of the American Chemical Society',
            year=1992
        ),
        method='UFF',
        atom_description='Generic carbon (sp3)'
    )
    db.add_vdw_parameters(c_vdw)
    
    # Save
    db.save()
    
    # Test retrieval
    params = db.get_bond_parameters('C', 'C', order=1)
    if params:
        print(f"\nC-C bond parameters:")
        print(f"  D_e = {params.D_e} kJ/mol")
        print(f"  r_e = {params.r_e} Angstrom")
        print(f"  Citation: {params.source.format_apa()}")
    
    # Statistics
    stats = db.get_statistics()
    print(f"\nDatabase statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value}")

