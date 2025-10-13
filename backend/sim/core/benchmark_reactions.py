"""
Benchmark Reactions Module
===========================

Loads and validates literature data for benchmark prebiotic chemistry reactions.

Provides easy access to reaction conditions, expected yields, and reference data.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ReactionConditions:
    """Reaction conditions from literature"""
    temperature: float  # K
    pH: Optional[float] = None
    pH_min: Optional[float] = None
    pH_max: Optional[float] = None
    solvent: str = "H2O"
    catalyst: Optional[str] = None
    
    def __repr__(self):
        pH_str = f"{self.pH}" if self.pH else f"{self.pH_min}-{self.pH_max}"
        return f"T={self.temperature}K, pH={pH_str}, solvent={self.solvent}"


@dataclass
class ProductData:
    """Expected product data from literature"""
    name: str
    formula: str
    smiles: Optional[str] = None
    yield_min: float = 0.0
    yield_max: float = 1.0
    detection_time_min: Optional[float] = None
    detection_time_max: Optional[float] = None
    note: Optional[str] = None
    
    @property
    def yield_range(self) -> Tuple[float, float]:
        return (self.yield_min, self.yield_max)
    
    @property
    def yield_mean(self) -> float:
        return (self.yield_min + self.yield_max) / 2.0
    
    def __repr__(self):
        yield_pct = f"{self.yield_mean*100:.1f}%"
        return f"{self.name} ({self.formula}): {yield_pct} yield"


@dataclass
class ReferenceData:
    """Literature reference"""
    doi: str
    authors: List[str]
    title: str
    journal: str
    year: int
    confidence: str = "medium"
    note: Optional[str] = None
    
    def __repr__(self):
        authors_str = ", ".join(self.authors[:2])
        if len(self.authors) > 2:
            authors_str += " et al."
        return f"{authors_str} ({self.year}). {self.title}. {self.journal}. DOI: {self.doi}"


class BenchmarkReactionDatabase:
    """
    Database of benchmark reactions for validation
    
    Loads literature data for:
    - Formose reaction
    - Strecker synthesis
    - HCN polymerization
    - Phosphorylation reactions
    """
    
    def __init__(self, data_path: str = "data/benchmark_reactions.json"):
        """
        Load benchmark reaction database
        
        Args:
            data_path: Path to JSON file with reaction data
        """
        self.data_path = Path(data_path)
        self.data = self._load_data()
        
        logger.info(f"Loaded benchmark reactions from {data_path}")
        logger.info(f"  Version: {self.data.get('version', 'unknown')}")
        logger.info(f"  Reactions: {len(self.data.get('reactions', {}))}")
    
    def _load_data(self) -> Dict:
        """Load JSON data"""
        if not self.data_path.exists():
            raise FileNotFoundError(f"Benchmark reactions data not found: {self.data_path}")
        
        with open(self.data_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    def get_reaction(self, name: str) -> Dict:
        """
        Get reaction data by name
        
        Args:
            name: Reaction name ('formose', 'strecker', 'hcn_polymerization', etc.)
        
        Returns:
            Dict with reaction data
        """
        reactions = self.data.get('reactions', {})
        if name not in reactions:
            available = list(reactions.keys())
            raise KeyError(f"Reaction '{name}' not found. Available: {available}")
        
        return reactions[name]
    
    def get_conditions(self, reaction_name: str) -> ReactionConditions:
        """Get reaction conditions as dataclass"""
        reaction = self.get_reaction(reaction_name)
        conditions_dict = reaction.get('conditions', {})
        
        return ReactionConditions(
            temperature=conditions_dict.get('temperature', 298),
            pH=conditions_dict.get('pH'),
            pH_min=conditions_dict.get('pH_min'),
            pH_max=conditions_dict.get('pH_max'),
            solvent=conditions_dict.get('solvent', 'H2O'),
            catalyst=conditions_dict.get('catalyst')
        )
    
    def get_products(self, reaction_name: str) -> List[ProductData]:
        """
        Get expected products for a reaction
        
        Returns list of ProductData objects
        """
        reaction = self.get_reaction(reaction_name)
        products_dict = reaction.get('products', {})
        
        products = []
        for product_name, product_info in products_dict.items():
            if 'formula' not in product_info:
                # Skip entries without formula (e.g. oligomers summary)
                continue
            
            product = ProductData(
                name=product_name,
                formula=product_info['formula'],
                smiles=product_info.get('smiles'),
                yield_min=product_info.get('yield_min', 0.0),
                yield_max=product_info.get('yield_max', 1.0),
                detection_time_min=product_info.get('detection_time_min'),
                detection_time_max=product_info.get('detection_time_max'),
                note=product_info.get('note')
            )
            products.append(product)
        
        return products
    
    def get_references(self, reaction_name: str) -> List[ReferenceData]:
        """Get literature references for a reaction"""
        reaction = self.get_reaction(reaction_name)
        refs_list = reaction.get('references', [])
        
        references = []
        for ref in refs_list:
            reference = ReferenceData(
                doi=ref['doi'],
                authors=ref['authors'],
                title=ref['title'],
                journal=ref['journal'],
                year=ref['year'],
                confidence=ref.get('confidence', 'medium'),
                note=ref.get('note')
            )
            references.append(reference)
        
        return references
    
    def get_mechanism(self, reaction_name: str) -> Dict[str, str]:
        """Get reaction mechanism steps"""
        reaction = self.get_reaction(reaction_name)
        return reaction.get('mechanism', {})
    
    def is_autocatalytic(self, reaction_name: str) -> bool:
        """Check if reaction is autocatalytic"""
        reaction = self.get_reaction(reaction_name)
        characteristics = reaction.get('characteristics', {})
        return characteristics.get('autocatalytic', False)
    
    def get_all_reactions(self) -> List[str]:
        """Get list of all available reaction names"""
        return list(self.data.get('reactions', {}).keys())
    
    def summary(self) -> str:
        """Generate human-readable summary"""
        lines = []
        lines.append("=" * 70)
        lines.append("BENCHMARK REACTIONS DATABASE")
        lines.append("=" * 70)
        lines.append(f"Version: {self.data.get('version', 'unknown')}")
        lines.append(f"Last updated: {self.data.get('last_updated', 'unknown')}")
        lines.append("")
        
        for reaction_name in self.get_all_reactions():
            reaction = self.get_reaction(reaction_name)
            lines.append(f"[{reaction_name.upper()}]")
            lines.append(f"  Name: {reaction['name']}")
            lines.append(f"  Reaction: {reaction['reaction']}")
            
            conditions = self.get_conditions(reaction_name)
            lines.append(f"  Conditions: {conditions}")
            
            products = self.get_products(reaction_name)
            lines.append(f"  Products ({len(products)}):")
            for product in products[:3]:  # Show first 3
                lines.append(f"    - {product}")
            
            refs = self.get_references(reaction_name)
            lines.append(f"  References: {len(refs)}")
            lines.append("")
        
        lines.append("=" * 70)
        return "\n".join(lines)


def validate_against_literature(
    product_name: str,
    observed_yield: float,
    expected_yield_range: Tuple[float, float],
    tolerance: float = 0.30
) -> bool:
    """
    Validate observed yield against literature range
    
    Args:
        product_name: Name of product
        observed_yield: Observed yield from simulation
        expected_yield_range: (min, max) from literature
        tolerance: Allowed deviation (default ±30%)
    
    Returns:
        True if observed yield is within tolerance of expected range
    """
    yield_min, yield_max = expected_yield_range
    yield_mean = (yield_min + yield_max) / 2.0
    
    # Allow ±tolerance deviation
    lower_bound = yield_min * (1 - tolerance)
    upper_bound = yield_max * (1 + tolerance)
    
    is_valid = lower_bound <= observed_yield <= upper_bound
    
    if not is_valid:
        logger.warning(
            f"Yield validation FAILED for {product_name}: "
            f"observed={observed_yield:.3f}, expected={yield_mean:.3f} "
            f"(range: {yield_min:.3f}-{yield_max:.3f}, tolerance=±{tolerance*100:.0f}%)"
        )
    else:
        logger.info(
            f"Yield validation PASSED for {product_name}: "
            f"observed={observed_yield:.3f} within expected range "
            f"{yield_min:.3f}-{yield_max:.3f}"
        )
    
    return is_valid


# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Load database
    db = BenchmarkReactionDatabase()
    
    # Print summary
    print(db.summary())
    
    # Example: Get formose reaction data
    print("\n" + "="*70)
    print("EXAMPLE: Formose Reaction")
    print("="*70)
    
    conditions = db.get_conditions('formose')
    print(f"Conditions: {conditions}")
    
    products = db.get_products('formose')
    print(f"\nExpected products:")
    for product in products:
        print(f"  {product}")
    
    refs = db.get_references('formose')
    print(f"\nReferences:")
    for ref in refs:
        print(f"  {ref}")
    
    # Test validation
    print("\n" + "="*70)
    print("EXAMPLE: Yield Validation")
    print("="*70)
    
    # Simulate a result
    observed_yield = 0.22  # 22% glycolaldehyde
    expected = products[0].yield_range  # Should be (0.15, 0.30)
    
    is_valid = validate_against_literature(
        product_name="glycolaldehyde",
        observed_yield=observed_yield,
        expected_yield_range=expected,
        tolerance=0.30
    )
    
    print(f"\nValidation result: {'PASS' if is_valid else 'FAIL'}")

