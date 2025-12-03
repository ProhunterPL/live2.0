"""
Truth-Filter: Główna klasa systemu walidacji wyników
"""
import json
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime

from backend.validation.types import (
    ValidationLevel, TruthResult, FilterResult, FilterStatus
)
from backend.validation.molecule_filter import MoleculeFilter
from backend.validation.literature_validator import LiteratureValidator
from backend.validation.simulation_quality import SimulationQualityFilter
from backend.validation.truth_report import TruthReportGenerator

logger = logging.getLogger(__name__)

# Import for match confidence (optional)
try:
    from matcher.confidence import MatchConfidenceEvaluator, Reliability
    MATCHER_AVAILABLE = True
except ImportError:
    MATCHER_AVAILABLE = False
    logger.warning("matcher.confidence not available, match confidence filtering disabled")


class TruthFilter:
    """
    Main Truth-Filter class for validating simulation results
    """
    
    def __init__(self,
                 validation_level: ValidationLevel = ValidationLevel.MEDIUM,
                 output_dir: Optional[str] = None):
        """
        Initialize Truth-Filter
        
        Args:
            validation_level: Level of validation (STRICT/MEDIUM/LENIENT)
            output_dir: Output directory for reports (optional)
        """
        self.validation_level = validation_level
        self.output_dir = Path(output_dir) if output_dir else None
        
        # Initialize filters
        self.molecule_filter = MoleculeFilter(
            strict_valence=(validation_level == ValidationLevel.STRICT),
            strict_charge=(validation_level == ValidationLevel.STRICT)
        )
        
        self.literature_validator = LiteratureValidator(
            validation_level=validation_level
        )
        
        # Set quality thresholds based on level
        if validation_level == ValidationLevel.STRICT:
            self.quality_filter = SimulationQualityFilter(
                min_completion_rate=1.0,
                max_energy_drift=0.005,  # 0.5%
                min_performance=1.0
            )
            self.max_energy_drift = 0.005
            self.max_momentum_drift = 0.001
        elif validation_level == ValidationLevel.MEDIUM:
            self.quality_filter = SimulationQualityFilter(
                min_completion_rate=0.95,
                max_energy_drift=0.01,  # 1%
                min_performance=0.5
            )
            self.max_energy_drift = 0.01
            self.max_momentum_drift = 0.005
        else:  # LENIENT
            self.quality_filter = SimulationQualityFilter(
                min_completion_rate=0.90,
                max_energy_drift=0.02,  # 2%
                min_performance=0.1
            )
            self.max_energy_drift = 0.02
            self.max_momentum_drift = 0.01
        
        # Initialize match confidence evaluator if available
        if MATCHER_AVAILABLE:
            if validation_level == ValidationLevel.STRICT:
                self.confidence_threshold = 0.8
            elif validation_level == ValidationLevel.MEDIUM:
                self.confidence_threshold = 0.6
            else:
                self.confidence_threshold = 0.4
            self.match_evaluator = MatchConfidenceEvaluator(
                confidence_threshold_high=self.confidence_threshold,
                confidence_threshold_medium=self.confidence_threshold * 0.75
            )
        else:
            self.match_evaluator = None
        
        # Initialize report generator
        self.report_generator = TruthReportGenerator(output_dir=output_dir)
    
    def _load_molecules(self, run_path: Path) -> List[Dict[str, Any]]:
        """
        Load molecules from run directory
        
        Args:
            run_path: Path to run directory
            
        Returns:
            List of molecule dicts
        """
        molecules = []
        
        # Try molecules.json
        molecules_path = run_path / "molecules.json"
        if molecules_path.exists():
            try:
                with open(molecules_path, 'r') as f:
                    data = json.load(f)
                    if isinstance(data, list):
                        molecules = data
                    elif isinstance(data, dict) and 'molecules' in data:
                        molecules = data['molecules']
            except Exception as e:
                logger.warning(f"Could not load molecules.json: {e}")
        
        # Try results.json
        if not molecules:
            results_path = run_path / "results.json"
            if results_path.exists():
                try:
                    with open(results_path, 'r') as f:
                        data = json.load(f)
                        molecules = data.get('molecules', data.get('molecules_detected', []))
                except Exception as e:
                    logger.warning(f"Could not load molecules from results.json: {e}")
        
        return molecules
    
    def _extract_scenario(self, run_path: Path) -> str:
        """
        Extract scenario name from run path
        
        Args:
            run_path: Path to run directory
            
        Returns:
            Scenario name
        """
        # Try to extract from path: results/phase2b_additional/miller_urey_extended/run_1
        parts = run_path.parts
        for part in reversed(parts):
            if 'miller_urey' in part.lower():
                return 'miller_urey_extended' if 'extended' in part else 'miller_urey'
            elif 'hydrothermal' in part.lower():
                return 'hydrothermal_extended' if 'extended' in part else 'hydrothermal'
            elif 'formamide' in part.lower():
                return 'formamide_extended' if 'extended' in part else 'formamide'
        
        # Default
        return 'unknown'
    
    def _validate_thermodynamics(self, run_path: Path) -> FilterResult:
        """
        Validate thermodynamics from log file
        
        Args:
            run_path: Path to run directory
            
        Returns:
            FilterResult
        """
        log_path = run_path / "simulation.log"
        warnings = []
        errors = []
        details = {}
        
        if not log_path.exists():
            warnings.append("simulation.log not found, skipping thermodynamic validation")
            return FilterResult(
                name="thermodynamics",
                status=FilterStatus.WARNING,
                details=details,
                warnings=warnings
            )
        
        # Try to parse energy/momentum drift from log
        # This is simplified - full implementation would parse more carefully
        try:
            with open(log_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Look for energy drift in last lines
            energy_drift = None
            momentum_drift = None
            
            for line in reversed(lines[-100:]):  # Check last 100 lines
                if 'energy' in line.lower() and 'drift' in line.lower():
                    # Try to extract number (look for percentage or small decimal)
                    import re
                    # Look for percentage: "X.XX%" or "X.XX %"
                    match = re.search(r'([\d.]+)\s*%', line)
                    if match:
                        energy_drift = float(match.group(1)) / 100.0
                        if energy_drift < 1.0:  # Sanity check: drift should be < 1.0
                            break
                    # Look for small decimal: "0.XXXX" or "X.XXe-X"
                    match = re.search(r'(0\.\d+|[\d.]+e-?\d+)', line)
                    if match:
                        energy_drift = float(match.group(1))
                        if energy_drift < 1.0:  # Sanity check
                            break
                
                if 'momentum' in line.lower() and 'drift' in line.lower():
                    import re
                    # Same pattern matching for momentum
                    match = re.search(r'([\d.]+)\s*%', line)
                    if match:
                        momentum_drift = float(match.group(1)) / 100.0
                        if momentum_drift < 1.0:
                            break
                    match = re.search(r'(0\.\d+|[\d.]+e-?\d+)', line)
                    if match:
                        momentum_drift = float(match.group(1))
                        if momentum_drift < 1.0:
                            break
            
            details['energy_drift'] = energy_drift
            details['momentum_drift'] = momentum_drift
            
            # Check thresholds
            if energy_drift is not None:
                if energy_drift > self.max_energy_drift:
                    errors.append(
                        f"Energy drift {energy_drift:.4f} exceeds threshold "
                        f"{self.max_energy_drift:.4f}"
                    )
                elif energy_drift > self.max_energy_drift * 0.8:
                    warnings.append(
                        f"Energy drift {energy_drift:.4f} is close to threshold"
                    )
            
            if momentum_drift is not None:
                if momentum_drift > self.max_momentum_drift:
                    errors.append(
                        f"Momentum drift {momentum_drift:.4f} exceeds threshold "
                        f"{self.max_momentum_drift:.4f}"
                    )
        
        except Exception as e:
            warnings.append(f"Error parsing thermodynamic data: {e}")
        
        # Determine status
        if len(errors) > 0:
            status = FilterStatus.FAIL
        elif len(warnings) > 0:
            status = FilterStatus.WARNING
        else:
            status = FilterStatus.PASS
        
        return FilterResult(
            name="thermodynamics",
            status=status,
            details=details,
            warnings=warnings,
            errors=errors
        )
    
    def _filter_matches(self, molecules: List[Dict[str, Any]]) -> FilterResult:
        """
        Filter molecules based on match confidence
        
        Args:
            molecules: List of molecules (with matches if available)
            
        Returns:
            FilterResult
        """
        if not MATCHER_AVAILABLE or not self.match_evaluator:
            return FilterResult(
                name="match_confidence",
                status=FilterStatus.PASS,
                details={'available': False},
                warnings=["Match confidence evaluator not available"]
            )
        
        high_confidence = []
        medium_confidence = []
        low_confidence = []
        warnings = []
        
        for mol in molecules:
            # Check if molecule has match data
            match_data = mol.get('match', mol.get('pubchem_match'))
            if not match_data:
                continue
            
            # Evaluate confidence
            try:
                confidence = self.match_evaluator.evaluate_match(mol, match_data)
                
                if confidence.reliability == Reliability.HIGH:
                    high_confidence.append(mol)
                elif confidence.reliability == Reliability.MEDIUM:
                    medium_confidence.append(mol)
                else:
                    low_confidence.append(mol)
                
                if confidence.warnings:
                    warnings.extend(confidence.warnings)
            
            except Exception as e:
                logger.warning(f"Error evaluating match confidence: {e}")
        
        details = {
            'high_confidence_matches': len(high_confidence),
            'medium_confidence_matches': len(medium_confidence),
            'low_confidence_matches': len(low_confidence),
            'total_matches': len(high_confidence) + len(medium_confidence) + len(low_confidence)
        }
        
        # Determine status
        if len(high_confidence) == 0 and details['total_matches'] > 0:
            status = FilterStatus.WARNING
            warnings.append("No high-confidence matches found")
        else:
            status = FilterStatus.PASS
        
        return FilterResult(
            name="match_confidence",
            status=status,
            details=details,
            warnings=warnings
        )
    
    def filter_run(self, run_path: str) -> TruthResult:
        """
        Filter single run
        
        Args:
            run_path: Path to run directory
            
        Returns:
            TruthResult with validation results
        """
        run_dir = Path(run_path)
        run_id = run_dir.name
        
        logger.info(f"Filtering run: {run_id}")
        
        filters = {}
        filtered_results = {
            'molecules': [],
            'reactions': [],
            'matches': []
        }
        
        # 1. Simulation Quality Check
        logger.info("  Checking simulation quality...")
        quality_result = self.quality_filter.check(str(run_dir))
        filters['simulation_quality'] = quality_result
        
        if quality_result.status == FilterStatus.FAIL:
            logger.warning(f"  Simulation quality check FAILED: {quality_result.errors}")
        
        # 2. Thermodynamic Validation
        logger.info("  Validating thermodynamics...")
        thermo_result = self._validate_thermodynamics(run_dir)
        filters['thermodynamics'] = thermo_result
        
        # 3. Load Molecules
        logger.info("  Loading molecules...")
        molecules = self._load_molecules(run_dir)
        logger.info(f"  Loaded {len(molecules)} molecules")
        
        # 4. Molecule Filtering
        logger.info("  Filtering molecules...")
        mol_result, filtered_molecules = self.molecule_filter.filter(molecules)
        filters['molecule_filter'] = mol_result
        filtered_results['molecules'] = filtered_molecules
        logger.info(f"  Filtered to {len(filtered_molecules)} real molecules")
        
        # 5. Literature Validation
        logger.info("  Validating against literature...")
        scenario = self._extract_scenario(run_dir)
        lit_result = self.literature_validator.validate(filtered_molecules, scenario)
        filters['literature_validator'] = lit_result
        
        # 6. Match Confidence Filter
        logger.info("  Filtering match confidence...")
        match_result = self._filter_matches(filtered_molecules)
        filters['match_confidence'] = match_result
        
        # Determine overall status
        if any(f.status == FilterStatus.FAIL for f in filters.values()):
            overall_status = FilterStatus.FAIL
        elif any(f.status == FilterStatus.WARNING for f in filters.values()):
            overall_status = FilterStatus.WARNING
        else:
            overall_status = FilterStatus.PASS
        
        # Create summary
        summary = {
            'total_molecules': len(molecules),
            'filtered_molecules': len(filtered_molecules),
            'retention_rate': len(filtered_molecules) / len(molecules) if molecules else 0.0,
            'scenario': scenario
        }
        
        result = TruthResult(
            run_id=run_id,
            run_path=str(run_dir),
            validation_level=self.validation_level,
            overall_status=overall_status,
            filters=filters,
            filtered_results=filtered_results,
            summary=summary,
            timestamp=datetime.now().timestamp()
        )
        
        logger.info(f"  Validation complete: {overall_status.value}")
        
        return result
    
    def filter_batch(self, results_dir: str, 
                     scenario: Optional[str] = None) -> List[TruthResult]:
        """
        Filter batch of runs
        
        Args:
            results_dir: Path to results directory
            scenario: Optional scenario filter (e.g., 'miller_urey_extended')
            
        Returns:
            List of TruthResult objects
        """
        results_path = Path(results_dir)
        results = []
        
        # Find all run directories
        if scenario:
            scenario_path = results_path / scenario
            if scenario_path.exists():
                run_dirs = [d for d in scenario_path.iterdir() if d.is_dir() and d.name.startswith('run_')]
            else:
                logger.warning(f"Scenario directory not found: {scenario_path}")
                return []
        else:
            # Search all scenarios
            run_dirs = []
            for scenario_dir in results_path.iterdir():
                if scenario_dir.is_dir():
                    for run_dir in scenario_dir.iterdir():
                        if run_dir.is_dir() and run_dir.name.startswith('run_'):
                            run_dirs.append(run_dir)
        
        logger.info(f"Found {len(run_dirs)} runs to filter")
        
        for run_dir in sorted(run_dirs):
            try:
                result = self.filter_run(str(run_dir))
                results.append(result)
            except Exception as e:
                logger.error(f"Error filtering {run_dir}: {e}")
        
        return results
    
    def generate_report(self, result: TruthResult, 
                       format: str = "both") -> Dict[str, Path]:
        """
        Generate report for result
        
        Args:
            result: TruthResult
            format: "json", "markdown", or "both"
            
        Returns:
            Dict with format -> path mapping
        """
        report_data = self.report_generator.generate_report(result)
        output_files = {}
        
        run_id = result.run_id.replace('/', '_')
        
        if format in ("json", "both"):
            json_path = self.report_generator.save_json(
                report_data,
                f"truth_report_{run_id}.json"
            )
            output_files['json'] = json_path
        
        if format in ("markdown", "both"):
            md_path = self.report_generator.save_markdown(
                result,
                f"truth_report_{run_id}.md"
            )
            output_files['markdown'] = md_path
        
        return output_files
    
    def generate_summary(self, results: List[TruthResult]) -> Dict[str, Any]:
        """
        Generate summary for multiple results
        
        Args:
            results: List of TruthResult objects
            
        Returns:
            Summary dict
        """
        return self.report_generator.generate_summary(results)

