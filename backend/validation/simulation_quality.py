"""
Simulation Quality Filter: Sprawdzanie jakości symulacji
"""
import json
import logging
import re
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
from backend.validation.types import FilterResult, FilterStatus

logger = logging.getLogger(__name__)


class SimulationQualityFilter:
    """
    Check simulation quality (completion, stability, performance)
    """
    
    def __init__(self, 
                 min_completion_rate: float = 1.0,
                 max_energy_drift: float = 0.01,
                 min_performance: float = 1.0):
        """
        Initialize quality filter
        
        Args:
            min_completion_rate: Minimum completion rate (1.0 = 100%)
            max_energy_drift: Maximum allowed energy drift
            min_performance: Minimum steps per second
        """
        self.min_completion_rate = min_completion_rate
        self.max_energy_drift = max_energy_drift
        self.min_performance = min_performance
    
    def _parse_log(self, log_path: Path) -> Dict[str, Any]:
        """
        Parse simulation log for quality metrics
        
        Returns:
            Dict with parsed metrics
        """
        metrics = {
            'completion_rate': 0.0,
            'final_step': 0,
            'max_steps': 0,
            'crashes': 0,
            'errors': [],
            'energy_drift': None,
            'performance': None
        }
        
        if not log_path.exists():
            return metrics
        
        try:
            with open(log_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Find final step
            for line in reversed(lines):
                # Look for step information
                step_match = re.search(r'step[:\s]+(\d+)', line, re.IGNORECASE)
                if step_match:
                    metrics['final_step'] = int(step_match.group(1))
                    break
            
            # Find max steps from config or log
            for line in lines:
                max_steps_match = re.search(r'max[_\s]?steps[:\s]+(\d+)', line, re.IGNORECASE)
                if max_steps_match:
                    metrics['max_steps'] = int(max_steps_match.group(1))
                    break
            
            # Count crashes/errors
            for line in lines:
                if 'error' in line.lower() or 'exception' in line.lower():
                    metrics['errors'].append(line.strip())
                    if 'crash' in line.lower() or 'fatal' in line.lower():
                        metrics['crashes'] += 1
            
            # Calculate completion rate
            if metrics['max_steps'] > 0:
                metrics['completion_rate'] = metrics['final_step'] / metrics['max_steps']
            
            # Try to find energy drift
            for line in reversed(lines):
                energy_match = re.search(r'energy[_\s]?drift[:\s]+([\d.]+)', line, re.IGNORECASE)
                if energy_match:
                    metrics['energy_drift'] = float(energy_match.group(1))
                    break
            
            # Try to find performance
            for line in reversed(lines):
                perf_match = re.search(r'steps?[_\s]?per[_\s]?sec[:\s]+([\d.]+)', line, re.IGNORECASE)
                if perf_match:
                    metrics['performance'] = float(perf_match.group(1))
                    break
        
        except Exception as e:
            logger.warning(f"Error parsing log {log_path}: {e}")
        
        return metrics
    
    def _parse_results_json(self, results_path: Path) -> Dict[str, Any]:
        """
        Parse results.json for quality metrics
        
        Returns:
            Dict with metrics from results.json
        """
        metrics = {}
        
        if not results_path.exists():
            return metrics
        
        try:
            with open(results_path, 'r') as f:
                data = json.load(f)
            
            # Extract relevant metrics
            metrics['final_step'] = data.get('final_step', data.get('step', 0))
            metrics['max_steps'] = data.get('max_steps', data.get('simulation', {}).get('max_steps', 0))
            metrics['completion_rate'] = (
                metrics['final_step'] / metrics['max_steps'] 
                if metrics['max_steps'] > 0 else 0.0
            )
            
            # Check for stability metrics
            if 'stability' in data:
                metrics['stability'] = data['stability']
            
            # Check for performance metrics
            if 'performance' in data:
                metrics['performance'] = data['performance']
            
            # Check for energy metrics
            if 'energy' in data:
                energy_data = data['energy']
                if 'drift' in energy_data:
                    metrics['energy_drift'] = energy_data['drift']
        
        except Exception as e:
            logger.warning(f"Error parsing results.json {results_path}: {e}")
        
        return metrics
    
    def check(self, run_path: str) -> FilterResult:
        """
        Check simulation quality
        
        Args:
            run_path: Path to run directory
            
        Returns:
            FilterResult with quality check results
        """
        run_dir = Path(run_path)
        warnings = []
        errors = []
        
        # Check if run directory exists
        if not run_dir.exists():
            return FilterResult(
                name="simulation_quality",
                status=FilterStatus.FAIL,
                details={'error': 'Run directory does not exist'},
                errors=[f"Run directory not found: {run_path}"]
            )
        
        # Parse log
        log_path = run_dir / "simulation.log"
        log_metrics = self._parse_log(log_path)
        
        # Parse results.json
        results_path = run_dir / "results.json"
        results_metrics = self._parse_results_json(results_path)
        
        # Merge metrics (results.json takes precedence)
        metrics = {**log_metrics, **results_metrics}
        
        # Check completion
        completion_rate = metrics.get('completion_rate', 0.0)
        if completion_rate < self.min_completion_rate:
            errors.append(
                f"Completion rate {completion_rate:.1%} below minimum "
                f"{self.min_completion_rate:.1%}"
            )
        
        # Check crashes
        crashes = metrics.get('crashes', 0)
        if crashes > 0:
            errors.append(f"Simulation crashed {crashes} time(s)")
        
        # Check energy drift
        energy_drift = metrics.get('energy_drift')
        if energy_drift is not None and energy_drift > self.max_energy_drift:
            warnings.append(
                f"Energy drift {energy_drift:.3f} exceeds threshold "
                f"{self.max_energy_drift:.3f}"
            )
        
        # Check performance
        performance = metrics.get('performance')
        if performance is not None and performance < self.min_performance:
            warnings.append(
                f"Performance {performance:.2f} steps/sec below minimum "
                f"{self.min_performance:.2f}"
            )
        
        # Determine status
        if len(errors) > 0:
            status = FilterStatus.FAIL
        elif len(warnings) > 0:
            status = FilterStatus.WARNING
        else:
            status = FilterStatus.PASS
        
        details = {
            'completion_rate': completion_rate,
            'final_step': metrics.get('final_step', 0),
            'max_steps': metrics.get('max_steps', 0),
            'crashes': crashes,
            'energy_drift': energy_drift,
            'performance': performance,
            'stability': metrics.get('stability', 'unknown')
        }
        
        return FilterResult(
            name="simulation_quality",
            status=status,
            details=details,
            warnings=warnings,
            errors=errors
        )

