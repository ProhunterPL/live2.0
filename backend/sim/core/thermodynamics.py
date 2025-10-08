"""
Thermodynamic validation module for Live 2.0
Validates fundamental physical laws in simulation
"""

import numpy as np
import taichi as ti
import time
import logging
from typing import Dict, List, Optional, Any, Tuple
from scipy import stats
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class ValidationResult:
    """Result of a thermodynamic validation test"""
    passed: bool
    error: float
    details: Dict[str, Any]
    timestamp: float
    step: int

class ThermodynamicValidator:
    """Walidator praw termodynamiki i mechaniki statystycznej"""
    
    def __init__(self, config):
        self.tolerance_energy = getattr(config, 'energy_tolerance', 1e-3)  # 0.1%
        self.tolerance_momentum = getattr(config, 'momentum_tolerance', 1e-4)  # 0.01%
        self.boltzmann_bins = 50
        self.validation_log = []
        
        # Physical constants
        self.k_B = 1.380649e-23  # Boltzmann constant (J/K)
        self.N_A = 6.02214076e23  # Avogadro's number
        
        logger.info(f"ThermodynamicValidator initialized with energy_tolerance={self.tolerance_energy}")
    
    def validate_energy_conservation(self, state_before, state_after, 
                                   energy_injected: float, energy_dissipated: float,
                                   step: int) -> ValidationResult:
        """
        Sprawdź: E_after = E_before + E_injected - E_dissipated ± ε
        
        Args:
            state_before: Simulation state before step
            state_after: Simulation state after step  
            energy_injected: Energy added during step
            energy_dissipated: Energy lost during step
            step: Current step number
            
        Returns:
            ValidationResult with energy conservation test results
        """
        try:
            E_before = self.compute_total_energy(state_before)
            E_after = self.compute_total_energy(state_after)
            
            expected = E_before + energy_injected - energy_dissipated
            actual = E_after
            
            # Relative error
            if abs(expected) > 1e-10:
                error = abs(actual - expected) / abs(expected)
            else:
                error = abs(actual - expected)
            
            passed = error < self.tolerance_energy
            
            details = {
                'E_before': E_before,
                'E_after': E_after,
                'E_expected': expected,
                'E_injected': energy_injected,
                'E_dissipated': energy_dissipated,
                'absolute_error': abs(actual - expected),
                'relative_error': error,
                'tolerance': self.tolerance_energy
            }
            
            result = ValidationResult(
                passed=passed,
                error=error,
                details=details,
                timestamp=time.time(),
                step=step
            )
            
            if not passed:
                logger.warning(f"Energy conservation violation at step {step}: "
                             f"error={error:.2e} > tolerance={self.tolerance_energy:.2e}")
            
            return result
            
        except Exception as e:
            logger.error(f"Energy conservation validation failed: {e}")
            return ValidationResult(
                passed=False,
                error=float('inf'),
                details={'error': str(e)},
                timestamp=time.time(),
                step=step
            )
    
    def validate_momentum_conservation(self, state_before, state_after, step: int) -> ValidationResult:
        """
        W izolowanym systemie: Σ(m·v) = const
        
        Args:
            state_before: Simulation state before step
            state_after: Simulation state after step
            step: Current step number
            
        Returns:
            ValidationResult with momentum conservation test results
        """
        try:
            p_before = self.compute_total_momentum(state_before)
            p_after = self.compute_total_momentum(state_after)
            
            dp = np.linalg.norm(p_after - p_before)
            p_total = np.linalg.norm(p_before) + 1e-10
            
            relative_error = dp / p_total
            passed = relative_error < self.tolerance_momentum
            
            details = {
                'p_before': p_before.tolist(),
                'p_after': p_after.tolist(),
                'dp_magnitude': dp,
                'p_total': p_total,
                'relative_error': relative_error,
                'tolerance': self.tolerance_momentum
            }
            
            result = ValidationResult(
                passed=passed,
                error=relative_error,
                details=details,
                timestamp=time.time(),
                step=step
            )
            
            if not passed:
                logger.warning(f"Momentum conservation violation at step {step}: "
                             f"error={relative_error:.2e} > tolerance={self.tolerance_momentum:.2e}")
            
            return result
            
        except Exception as e:
            logger.error(f"Momentum conservation validation failed: {e}")
            return ValidationResult(
                passed=False,
                error=float('inf'),
                details={'error': str(e)},
                timestamp=time.time(),
                step=step
            )
    
    def validate_maxwell_boltzmann(self, velocities: np.ndarray, temperature: float, 
                                 step: int) -> ValidationResult:
        """
        Sprawdź, czy rozkład prędkości odpowiada rozkładowi M-B
        
        Args:
            velocities: Array of particle velocities (N x 2)
            temperature: System temperature (K)
            step: Current step number
            
        Returns:
            ValidationResult with M-B distribution test results
        """
        try:
            if len(velocities) < 10:
                # Not enough particles for statistical test
                return ValidationResult(
                    passed=True,
                    error=0.0,
                    details={'note': 'Insufficient particles for M-B test'},
                    timestamp=time.time(),
                    step=step
                )
            
            # Calculate speeds
            speeds = np.linalg.norm(velocities, axis=1)
            
            # Remove outliers (speeds > 5σ)
            mean_speed = np.mean(speeds)
            std_speed = np.std(speeds)
            valid_speeds = speeds[speeds < mean_speed + 5 * std_speed]
            
            if len(valid_speeds) < 10:
                return ValidationResult(
                    passed=True,
                    error=0.0,
                    details={'note': 'Too many outliers for M-B test'},
                    timestamp=time.time(),
                    step=step
                )
            
            # Histogram empirical
            hist, bin_edges = np.histogram(valid_speeds, bins=self.boltzmann_bins, 
                                          density=True)
            
            # Theoretical M-B distribution
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
            # Maxwell-Boltzmann PDF: f(v) = (m/(2πkT))^(3/2) * 4πv² * exp(-mv²/(2kT))
            # For 2D: f(v) = (m/(2πkT)) * v * exp(-mv²/(2kT))
            mass = 1.0  # Assume unit mass for now
            theoretical = self.maxwell_boltzmann_pdf_2d(bin_centers, temperature, mass)
            
            # Normalize theoretical to match histogram
            theoretical = theoretical * np.sum(hist) / np.sum(theoretical)
            
            # Chi-square goodness of fit
            # Remove bins with zero theoretical probability
            valid_bins = theoretical > 1e-10
            if np.sum(valid_bins) < 5:
                return ValidationResult(
                    passed=True,
                    error=0.0,
                    details={'note': 'Insufficient valid bins for chi-square test'},
                    timestamp=time.time(),
                    step=step
                )
            
            chi2, p_value = stats.chisquare(hist[valid_bins], theoretical[valid_bins])
            
            passed = p_value > 0.05  # 95% confidence
            
            details = {
                'temperature': temperature,
                'n_particles': len(valid_speeds),
                'mean_speed': np.mean(valid_speeds),
                'std_speed': np.std(valid_speeds),
                'chi2': chi2,
                'p_value': p_value,
                'histogram': hist.tolist(),
                'theoretical': theoretical.tolist(),
                'bin_centers': bin_centers.tolist()
            }
            
            result = ValidationResult(
                passed=passed,
                error=1.0 - p_value,  # Error = 1 - p_value
                details=details,
                timestamp=time.time(),
                step=step
            )
            
            if not passed:
                logger.warning(f"Maxwell-Boltzmann violation at step {step}: "
                             f"p_value={p_value:.3f} < 0.05")
            
            return result
            
        except Exception as e:
            logger.error(f"Maxwell-Boltzmann validation failed: {e}")
            return ValidationResult(
                passed=False,
                error=1.0,
                details={'error': str(e)},
                timestamp=time.time(),
                step=step
            )
    
    def validate_second_law(self, state_before, state_after, step: int) -> ValidationResult:
        """
        II zasada termodynamiki: ΔS ≥ 0 (dla izolowanego systemu)
        
        Args:
            state_before: Simulation state before step
            state_after: Simulation state after step
            step: Current step number
            
        Returns:
            ValidationResult with entropy change test results
        """
        try:
            S_before = self.compute_entropy(state_before)
            S_after = self.compute_entropy(state_after)
            
            delta_S = S_after - S_before
            
            # Allow small negative changes due to numerical errors
            passed = delta_S >= -1e-6
            
            details = {
                'S_before': S_before,
                'S_after': S_after,
                'delta_S': delta_S,
                'tolerance': -1e-6
            }
            
            result = ValidationResult(
                passed=passed,
                error=max(0, -delta_S),  # Error is magnitude of negative change
                details=details,
                timestamp=time.time(),
                step=step
            )
            
            if not passed:
                logger.warning(f"Second law violation at step {step}: "
                             f"ΔS={delta_S:.2e} < 0")
            
            return result
            
        except Exception as e:
            logger.error(f"Second law validation failed: {e}")
            return ValidationResult(
                passed=False,
                error=float('inf'),
                details={'error': str(e)},
                timestamp=time.time(),
                step=step
            )
    
    def compute_total_energy(self, state) -> float:
        """Calculate total energy of the system"""
        try:
            total_energy = 0.0
            
            # Energy field energy
            if hasattr(state, 'energy_field'):
                energy_field = state.energy_field.to_numpy()
                total_energy += float(energy_field.sum())
            
            # Particle kinetic energy
            if hasattr(state, 'velocities') and hasattr(state, 'attributes'):
                velocities = state.velocities.to_numpy()
                attributes = state.attributes.to_numpy()
                active = state.active.to_numpy()
                
                for i in range(len(velocities)):
                    if active[i] == 1:
                        mass = attributes[i][0]  # First component is mass
                        vx, vy = velocities[i]
                        kinetic = 0.5 * mass * (vx * vx + vy * vy)
                        total_energy += kinetic
            
            # Bond potential energy (if available)
            if hasattr(state, 'bond_energy'):
                total_energy += float(state.bond_energy)
            
            return total_energy
            
        except Exception as e:
            logger.error(f"Energy calculation failed: {e}")
            return 0.0
    
    def compute_total_momentum(self, state) -> np.ndarray:
        """Calculate total momentum of the system"""
        try:
            total_momentum = np.zeros(2)
            
            if hasattr(state, 'velocities') and hasattr(state, 'attributes'):
                velocities = state.velocities.to_numpy()
                attributes = state.attributes.to_numpy()
                active = state.active.to_numpy()
                
                for i in range(len(velocities)):
                    if active[i] == 1:
                        mass = attributes[i][0]
                        vx, vy = velocities[i]
                        total_momentum[0] += mass * vx
                        total_momentum[1] += mass * vy
            
            return total_momentum
            
        except Exception as e:
            logger.error(f"Momentum calculation failed: {e}")
            return np.zeros(2)
    
    def compute_entropy(self, state) -> float:
        """
        Entropia konfiguracyjna (Shannon) + wkład kinetyczny
        """
        try:
            entropy = 0.0
            
            # Configurational entropy (grid-based)
            if hasattr(state, 'positions'):
                positions = state.positions.to_numpy()
                active = state.active.to_numpy()
                
                # Create density grid
                grid_size = 32  # Coarse grid for entropy calculation
                density_grid = np.zeros((grid_size, grid_size))
                
                for i in range(len(positions)):
                    if active[i] == 1:
                        x, y = positions[i]
                        # Map to grid coordinates
                        grid_x = int(x * grid_size / 128) % grid_size  # Assuming 128x128 world
                        grid_y = int(y * grid_size / 128) % grid_size
                        density_grid[grid_y, grid_x] += 1
                
                # Normalize to probabilities
                total_particles = np.sum(density_grid)
                if total_particles > 0:
                    p = density_grid / total_particles
                    p = p[p > 0]  # Remove zeros
                    S_config = -np.sum(p * np.log(p))
                    entropy += S_config
            
            # Kinetic entropy (from velocity distribution)
            if hasattr(state, 'velocities'):
                velocities = state.velocities.to_numpy()
                active = state.active.to_numpy()
                
                valid_velocities = velocities[active == 1]
                if len(valid_velocities) > 0:
                    T = self.compute_temperature(valid_velocities)
                    if T > 0:
                        # S_kinetic = N * (3/2 * ln(T) + const) for 2D
                        N = len(valid_velocities)
                        S_kinetic = N * (1.0 * np.log(T) + 1.0)  # Simplified for 2D
                        entropy += S_kinetic
            
            return entropy
            
        except Exception as e:
            logger.error(f"Entropy calculation failed: {e}")
            return 0.0
    
    def compute_temperature(self, velocities: np.ndarray) -> float:
        """Compute temperature from velocity distribution"""
        try:
            if len(velocities) < 2:
                return 0.0
            
            # For 2D: T = m * <v²> / (2 * k_B)
            # Using units where k_B = 1
            speeds_squared = np.sum(velocities**2, axis=1)
            mean_speed_squared = np.mean(speeds_squared)
            
            # Assume unit mass
            temperature = mean_speed_squared / 2.0
            
            return temperature
            
        except Exception as e:
            logger.error(f"Temperature calculation failed: {e}")
            return 0.0
    
    def maxwell_boltzmann_pdf_2d(self, speeds: np.ndarray, temperature: float, mass: float) -> np.ndarray:
        """
        2D Maxwell-Boltzmann distribution: f(v) = (m/(2πkT)) * v * exp(-mv²/(2kT))
        """
        try:
            if temperature <= 0:
                return np.zeros_like(speeds)
            
            # Normalization constant
            norm = mass / (2 * np.pi * temperature)
            
            # Exponential factor
            exp_factor = np.exp(-mass * speeds**2 / (2 * temperature))
            
            # PDF
            pdf = norm * speeds * exp_factor
            
            return pdf
            
        except Exception as e:
            logger.error(f"M-B PDF calculation failed: {e}")
            return np.zeros_like(speeds)
    
    def validate_all(self, state_before, state_after, energy_injected: float, 
                    energy_dissipated: float, step: int) -> Dict[str, ValidationResult]:
        """
        Run all thermodynamic validations
        
        Returns:
            Dictionary with all validation results
        """
        results = {}
        
        # Energy conservation
        results['energy'] = self.validate_energy_conservation(
            state_before, state_after, energy_injected, energy_dissipated, step
        )
        
        # Momentum conservation
        results['momentum'] = self.validate_momentum_conservation(
            state_before, state_after, step
        )
        
        # Maxwell-Boltzmann distribution
        if hasattr(state_after, 'velocities'):
            velocities = state_after.velocities.to_numpy()
            active = state_after.active.to_numpy()
            valid_velocities = velocities[active == 1]
            
            if len(valid_velocities) > 0:
                temperature = self.compute_temperature(valid_velocities)
                results['maxwell_boltzmann'] = self.validate_maxwell_boltzmann(
                    valid_velocities, temperature, step
                )
            else:
                results['maxwell_boltzmann'] = ValidationResult(
                    passed=True, error=0.0, details={'note': 'No active particles'},
                    timestamp=time.time(), step=step
                )
        
        # Second law (entropy)
        results['second_law'] = self.validate_second_law(
            state_before, state_after, step
        )
        
        # Overall result
        all_passed = all(r.passed for r in results.values())
        results['all_passed'] = ValidationResult(
            passed=all_passed,
            error=0.0,
            details={'individual_results': {k: r.passed for k, r in results.items()}},
            timestamp=time.time(),
            step=step
        )
        
        # Log results
        if not all_passed:
            failed_tests = [k for k, r in results.items() if not r.passed and k != 'all_passed']
            logger.warning(f"Thermodynamic validation failed at step {step}: {failed_tests}")
        
        return results
    
    def log_validation_results(self, results: Dict[str, ValidationResult]):
        """Log validation results to file"""
        log_entry = {
            'timestamp': time.time(),
            'step': results['all_passed'].step,
            'all_passed': results['all_passed'].passed,
            'results': {}
        }
        
        for test_name, result in results.items():
            if test_name != 'all_passed':
                log_entry['results'][test_name] = {
                    'passed': result.passed,
                    'error': result.error,
                    'details': result.details
                }
        
        self.validation_log.append(log_entry)
    
    def get_validation_summary(self) -> Dict[str, Any]:
        """Get summary statistics of validation results"""
        if not self.validation_log:
            return {'message': 'No validation data available'}
        
        total_tests = len(self.validation_log)
        passed_tests = sum(1 for entry in self.validation_log if entry['all_passed'])
        
        # Per-test statistics
        test_stats = {}
        for test_name in ['energy', 'momentum', 'maxwell_boltzmann', 'second_law']:
            test_results = [entry['results'].get(test_name, {}).get('passed', True) 
                          for entry in self.validation_log]
            test_stats[test_name] = {
                'passed': sum(test_results),
                'total': len(test_results),
                'success_rate': sum(test_results) / len(test_results) if test_results else 1.0
            }
        
        return {
            'total_tests': total_tests,
            'overall_success_rate': passed_tests / total_tests,
            'test_statistics': test_stats,
            'latest_entry': self.validation_log[-1] if self.validation_log else None
        }
