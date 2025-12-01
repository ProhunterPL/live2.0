"""
Type definitions for Truth-Filter
"""
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from enum import Enum


class ValidationLevel(Enum):
    """Poziomy walidacji"""
    STRICT = "STRICT"    # Dla publikacji - najwyższe wymagania
    MEDIUM = "MEDIUM"    # Dla analizy - umiarkowane wymagania
    LENIENT = "LENIENT"  # Dla eksploracji - podstawowe wymagania


class FilterStatus(Enum):
    """Status filtru"""
    PASS = "PASS"
    WARNING = "WARNING"
    FAIL = "FAIL"


@dataclass
class FilterResult:
    """Wynik pojedynczego filtru"""
    name: str
    status: FilterStatus
    details: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    
    def is_pass(self) -> bool:
        """Czy filtr przeszedł"""
        return self.status == FilterStatus.PASS
    
    def has_warnings(self) -> bool:
        """Czy są ostrzeżenia"""
        return len(self.warnings) > 0


@dataclass
class TruthResult:
    """Kompletny wynik walidacji runa"""
    run_id: str
    run_path: str
    validation_level: ValidationLevel
    overall_status: FilterStatus
    filters: Dict[str, FilterResult] = field(default_factory=dict)
    filtered_results: Dict[str, Any] = field(default_factory=dict)
    summary: Dict[str, Any] = field(default_factory=dict)
    timestamp: Optional[float] = None
    
    def is_pass(self) -> bool:
        """Czy run przeszedł walidację"""
        return self.overall_status == FilterStatus.PASS
    
    def get_warnings(self) -> List[str]:
        """Zbierz wszystkie ostrzeżenia"""
        warnings = []
        for filter_result in self.filters.values():
            warnings.extend(filter_result.warnings)
        return warnings
    
    def get_errors(self) -> List[str]:
        """Zbierz wszystkie błędy"""
        errors = []
        for filter_result in self.filters.values():
            errors.extend(filter_result.errors)
        return errors

