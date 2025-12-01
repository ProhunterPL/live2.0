"""
Literature Validator: Weryfikacja zgodności z literaturą
"""
import logging
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
from backend.validation.types import FilterResult, FilterStatus, ValidationLevel
from backend.sim.core.benchmark_reactions import BenchmarkReactionDatabase, validate_against_literature

logger = logging.getLogger(__name__)


# Mapowanie scenariuszy na reakcje benchmarkowe
SCENARIO_TO_REACTIONS = {
    'miller_urey_extended': ['miller_urey', 'formose', 'strecker'],
    'miller_urey': ['miller_urey', 'formose'],
    'hydrothermal_extended': ['hydrothermal', 'strecker'],
    'hydrothermal': ['hydrothermal'],
    'formamide_extended': ['formamide', 'formose'],
    'formamide': ['formamide'],
}


class LiteratureValidator:
    """
    Validate simulation results against literature benchmarks
    """
    
    def __init__(self, 
                 validation_level: ValidationLevel = ValidationLevel.MEDIUM,
                 data_path: Optional[str] = None):
        """
        Initialize literature validator
        
        Args:
            validation_level: Level of validation (STRICT/MEDIUM/LENIENT)
            data_path: Path to benchmark reactions JSON (optional)
        """
        self.validation_level = validation_level
        
        # Set tolerance based on validation level
        if validation_level == ValidationLevel.STRICT:
            self.yield_tolerance = 0.20  # ±20%
            self.min_expected_products_rate = 0.5  # 50% of expected products
        elif validation_level == ValidationLevel.MEDIUM:
            self.yield_tolerance = 0.30  # ±30%
            self.min_expected_products_rate = 0.4  # 40% of expected products
        else:  # LENIENT
            self.yield_tolerance = 0.50  # ±50%
            self.min_expected_products_rate = 0.3  # 30% of expected products
        
        # Load benchmark database
        if data_path is None:
            # Try to find data file
            possible_paths = [
                "data/benchmark_reactions.json",
                "backend/sim/core/data/benchmark_reactions.json",
                Path(__file__).parent.parent / "sim" / "core" / "data" / "benchmark_reactions.json"
            ]
            data_path = None
            for path in possible_paths:
                if Path(path).exists():
                    data_path = path
                    break
        
        if data_path and Path(data_path).exists():
            try:
                self.db = BenchmarkReactionDatabase(data_path)
            except Exception as e:
                logger.warning(f"Could not load benchmark database from {data_path}: {e}")
                self.db = None
        else:
            logger.warning("Benchmark reactions database not found, literature validation will be limited")
            self.db = None
    
    def _get_expected_products(self, scenario: str) -> List[Dict[str, Any]]:
        """
        Get expected products for scenario
        
        Args:
            scenario: Scenario name (e.g., 'miller_urey_extended')
            
        Returns:
            List of expected product dicts
        """
        if self.db is None:
            return []
        
        # Get reactions for scenario
        reaction_names = SCENARIO_TO_REACTIONS.get(scenario, [])
        if not reaction_names:
            # Try to extract base scenario name
            base_scenario = scenario.replace('_extended', '').replace('_extended', '')
            reaction_names = SCENARIO_TO_REACTIONS.get(base_scenario, [])
        
        expected_products = []
        for reaction_name in reaction_names:
            try:
                products = self.db.get_products(reaction_name)
                for product in products:
                    expected_products.append({
                        'name': product.name,
                        'formula': product.formula,
                        'yield_range': product.yield_range,
                        'yield_mean': product.yield_mean,
                        'reaction': reaction_name
                    })
            except Exception as e:
                logger.warning(f"Could not get products for reaction {reaction_name}: {e}")
        
        return expected_products
    
    def _detect_products(self, molecules: List[Dict[str, Any]], 
                         expected_products: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Detect which expected products are present in molecules
        
        Args:
            molecules: List of detected molecules
            expected_products: List of expected products
            
        Returns:
            Dict with detection results
        """
        detected = []
        not_detected = []
        
        # Create formula lookup
        molecule_formulas = {}
        for mol in molecules:
            formula = mol.get('formula', '')
            if formula:
                molecule_formulas[formula.upper()] = mol
        
        # Check each expected product
        for expected in expected_products:
            formula = expected['formula'].upper()
            if formula in molecule_formulas:
                detected.append({
                    'name': expected['name'],
                    'formula': expected['formula'],
                    'reaction': expected['reaction'],
                    'molecule': molecule_formulas[formula]
                })
            else:
                not_detected.append(expected)
        
        return {
            'detected': detected,
            'not_detected': not_detected,
            'detection_rate': len(detected) / len(expected_products) if expected_products else 0.0
        }
    
    def _validate_yields(self, detected_products: List[Dict[str, Any]], 
                        total_molecules: int) -> Dict[str, Any]:
        """
        Validate yields against literature
        
        Args:
            detected_products: List of detected products
            total_molecules: Total number of molecules in simulation
            
        Returns:
            Dict with yield validation results
        """
        yield_validations = []
        valid_yields = 0
        invalid_yields = 0
        
        for product in detected_products:
            expected = product.get('expected', {})
            yield_range = expected.get('yield_range', (0.0, 1.0))
            
            # Calculate observed yield (simplified - would need actual counts)
            observed_yield = 1.0 / total_molecules if total_molecules > 0 else 0.0
            
            # Validate against literature
            is_valid = validate_against_literature(
                product['name'],
                observed_yield,
                yield_range,
                tolerance=self.yield_tolerance
            )
            
            yield_validations.append({
                'product': product['name'],
                'observed_yield': observed_yield,
                'expected_range': yield_range,
                'is_valid': is_valid
            })
            
            if is_valid:
                valid_yields += 1
            else:
                invalid_yields += 1
        
        yield_match_rate = valid_yields / len(yield_validations) if yield_validations else 0.0
        
        return {
            'yield_validations': yield_validations,
            'valid_yields': valid_yields,
            'invalid_yields': invalid_yields,
            'yield_match_rate': yield_match_rate
        }
    
    def validate(self, molecules: List[Dict[str, Any]], 
                 scenario: str) -> FilterResult:
        """
        Validate molecules against literature
        
        Args:
            molecules: List of detected molecules
            scenario: Scenario name (e.g., 'miller_urey_extended')
            
        Returns:
            FilterResult with validation results
        """
        warnings = []
        errors = []
        
        if self.db is None:
            warnings.append("Benchmark database not available, skipping literature validation")
            return FilterResult(
                name="literature_validator",
                status=FilterStatus.WARNING,
                details={'available': False},
                warnings=warnings
            )
        
        # Get expected products
        expected_products = self._get_expected_products(scenario)
        
        if not expected_products:
            warnings.append(f"No expected products found for scenario: {scenario}")
            return FilterResult(
                name="literature_validator",
                status=FilterStatus.WARNING,
                details={'expected_products_count': 0},
                warnings=warnings
            )
        
        # Detect products
        detection_results = self._detect_products(molecules, expected_products)
        detected_count = len(detection_results['detected'])
        detection_rate = detection_results['detection_rate']
        
        # Validate yields (if we have molecule counts)
        total_molecules = len(molecules)
        yield_results = self._validate_yields(detection_results['detected'], total_molecules)
        
        # Check if meets minimum requirements
        if detection_rate < self.min_expected_products_rate:
            errors.append(
                f"Detection rate {detection_rate:.1%} below minimum "
                f"{self.min_expected_products_rate:.1%} for {self.validation_level.value} level"
            )
        
        # Determine status
        if len(errors) > 0:
            status = FilterStatus.FAIL
        elif detection_rate < self.min_expected_products_rate * 1.2:  # 20% buffer
            status = FilterStatus.WARNING
            warnings.append(f"Detection rate {detection_rate:.1%} is close to minimum")
        else:
            status = FilterStatus.PASS
        
        details = {
            'expected_products_count': len(expected_products),
            'detected_products_count': detected_count,
            'detection_rate': detection_rate,
            'yield_match_rate': yield_results['yield_match_rate'],
            'detected_products': [p['name'] for p in detection_results['detected']],
            'not_detected_products': [p['name'] for p in detection_results['not_detected']]
        }
        
        return FilterResult(
            name="literature_validator",
            status=status,
            details=details,
            warnings=warnings,
            errors=errors
        )

