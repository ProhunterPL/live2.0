"""
Reaction Kinetics Analysis Module
==================================

Analyzes reaction rates and equilibrium constants from simulation trajectories.

Provides:
- Rate constant calculation (k)
- Equilibrium constant calculation (K_eq)
- Activation energy estimation (E_a)
- Reaction order determination
"""

import numpy as np
import logging
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from scipy import optimize, stats

logger = logging.getLogger(__name__)


@dataclass
class KineticData:
    """Stores kinetic data for a reaction"""
    times: np.ndarray
    concentrations: Dict[str, np.ndarray]
    temperature: float


@dataclass
class RateConstant:
    """Rate constant with uncertainty"""
    value: float
    error: float
    units: str
    order: int
    r_squared: float
    
    def __repr__(self):
        return f"k = {self.value:.3e} +/- {self.error:.3e} {self.units} (order {self.order}, R^2={self.r_squared:.3f})"


@dataclass
class EquilibriumConstant:
    """Equilibrium constant with thermodynamic properties"""
    K_eq: float
    delta_G: float  # kJ/mol
    temperature: float  # K
    error: Optional[float] = None
    
    def __repr__(self):
        return f"K_eq = {self.K_eq:.3e}, DG = {self.delta_G:.2f} kJ/mol at {self.temperature}K"


class ReactionKineticsAnalyzer:
    """
    Analyzes reaction kinetics from concentration trajectories
    
    Usage:
        analyzer = ReactionKineticsAnalyzer()
        
        # Provide concentration vs time data
        kinetic_data = KineticData(
            times=times,
            concentrations={'A': [A], 'B': [B]},
            temperature=300
        )
        
        # Calculate rate constant
        rate_const = analyzer.fit_rate_constant(
            kinetic_data, reactant='A', order=1
        )
        
        # Calculate equilibrium constant
        K_eq = analyzer.calculate_equilibrium_constant(
            kinetic_data, reactants=['A'], products=['B']
        )
    """
    
    def __init__(self, R: float = 8.314):
        """
        Initialize kinetics analyzer
        
        Args:
            R: Gas constant (J/(mol·K))
        """
        self.R = R
    
    def fit_rate_constant(self,
                         data: KineticData,
                         reactant: str,
                         order: int = 1,
                         method: str = 'linear') -> RateConstant:
        """
        Fit rate constant from concentration trajectory
        
        For order n:
            -d[A]/dt = k[A]^n
        
        Integrated forms:
            n=0: [A] = [A]0 - kt
            n=1: ln([A]) = ln([A]0) - kt
            n=2: 1/[A] = 1/[A]0 + kt
        
        Args:
            data: KineticData with times and concentrations
            reactant: Which species to analyze
            order: Reaction order (0, 1, or 2)
            method: 'linear' or 'nonlinear' fitting
        
        Returns:
            RateConstant with fitted k value
        """
        times = data.times
        conc = data.concentrations[reactant]
        
        # Remove any invalid data (zeros, negatives)
        mask = conc > 0
        times = times[mask]
        conc = conc[mask]
        
        if len(times) < 3:
            raise ValueError(f"Insufficient data points for {reactant}")
        
        # Fit based on order
        if order == 0:
            # [A] = [A]0 - kt
            k, intercept, r_value, p_value, std_err = stats.linregress(times, conc)
            k = -k  # Rate is negative of slope
            units = "M/s"
            
        elif order == 1:
            # ln([A]) = ln([A]0) - kt
            ln_conc = np.log(conc)
            k, intercept, r_value, p_value, std_err = stats.linregress(times, ln_conc)
            k = -k
            units = "1/s"
            
        elif order == 2:
            # 1/[A] = 1/[A]0 + kt
            inv_conc = 1.0 / conc
            k, intercept, r_value, p_value, std_err = stats.linregress(times, inv_conc)
            units = "1/(M·s)"
            
        else:
            raise ValueError(f"Order {order} not supported (use 0, 1, or 2)")
        
        r_squared = r_value ** 2
        
        return RateConstant(
            value=k,
            error=std_err,
            units=units,
            order=order,
            r_squared=r_squared
        )
    
    def determine_reaction_order(self,
                                data: KineticData,
                                reactant: str,
                                max_order: int = 2) -> int:
        """
        Determine reaction order by trying different orders and comparing R^2
        
        Returns order (0, 1, or 2) with best fit
        """
        best_order = 1
        best_r_squared = 0.0
        
        for order in range(max_order + 1):
            try:
                rate_const = self.fit_rate_constant(data, reactant, order=order)
                if rate_const.r_squared > best_r_squared:
                    best_r_squared = rate_const.r_squared
                    best_order = order
            except Exception as e:
                logger.warning(f"Failed to fit order {order}: {e}")
        
        logger.info(f"Determined reaction order: {best_order} (R^2 = {best_r_squared:.3f})")
        return best_order
    
    def calculate_equilibrium_constant(self,
                                      data: KineticData,
                                      reactants: List[str],
                                      products: List[str],
                                      equilibrium_window: Tuple[float, float] = None) -> EquilibriumConstant:
        """
        Calculate equilibrium constant K_eq = [products] / [reactants]
        
        For reaction: aA + bB ⇌ cC + dD
            K_eq = [C]^c [D]^d / ([A]^a [B]^b)
        
        Args:
            data: KineticData with concentration trajectories
            reactants: List of reactant species names
            products: List of product species names
            equilibrium_window: (start_time, end_time) to average over
        
        Returns:
            EquilibriumConstant with K_eq and ΔG
        """
        # Determine equilibrium window (last 20% of simulation if not specified)
        if equilibrium_window is None:
            t_max = data.times[-1]
            equilibrium_window = (0.8 * t_max, t_max)
        
        # Select data in equilibrium window
        mask = (data.times >= equilibrium_window[0]) & (data.times <= equilibrium_window[1])
        
        if not np.any(mask):
            raise ValueError("No data in equilibrium window")
        
        # Calculate K_eq at each time point in window
        K_eq_values = []
        
        for i in np.where(mask)[0]:
            # Numerator: product of product concentrations
            numerator = 1.0
            for product in products:
                conc = data.concentrations[product][i]
                if conc > 0:
                    numerator *= conc
                else:
                    numerator = 0
                    break
            
            # Denominator: product of reactant concentrations
            denominator = 1.0
            for reactant in reactants:
                conc = data.concentrations[reactant][i]
                if conc > 0:
                    denominator *= conc
                else:
                    denominator = np.inf
                    break
            
            if numerator > 0 and np.isfinite(denominator) and denominator > 0:
                K_eq_values.append(numerator / denominator)
        
        if not K_eq_values:
            raise ValueError("Could not calculate K_eq (zero or invalid concentrations)")
        
        # Average K_eq over equilibrium window
        K_eq = np.mean(K_eq_values)
        K_eq_error = np.std(K_eq_values)
        
        # Calculate ΔG from K_eq
        # ΔG = -RT ln(K_eq)
        if K_eq > 0:
            delta_G = -self.R * data.temperature * np.log(K_eq) / 1000.0  # Convert to kJ/mol
        else:
            delta_G = np.inf
        
        return EquilibriumConstant(
            K_eq=K_eq,
            delta_G=delta_G,
            temperature=data.temperature,
            error=K_eq_error
        )
    
    def calculate_half_life(self, rate_constant: RateConstant, initial_conc: float = 1.0) -> float:
        """
        Calculate half-life from rate constant
        
        t_1/2:
            order 0: [A]0 / (2k)
            order 1: ln(2) / k
            order 2: 1 / (k[A]0)
        
        Args:
            rate_constant: RateConstant object
            initial_conc: Initial concentration (for order != 1)
        
        Returns:
            Half-life in same time units as rate constant
        """
        k = rate_constant.value
        order = rate_constant.order
        
        if order == 0:
            return initial_conc / (2 * k)
        elif order == 1:
            return np.log(2) / k
        elif order == 2:
            return 1.0 / (k * initial_conc)
        else:
            raise ValueError(f"Order {order} not supported")
    
    def estimate_activation_energy(self,
                                  rate_constants: List[Tuple[float, RateConstant]]) -> Tuple[float, float]:
        """
        Estimate activation energy from rate constants at different temperatures
        
        Arrhenius equation: k = A * exp(-E_a / RT)
        Linear form: ln(k) = ln(A) - E_a/(RT)
        
        Args:
            rate_constants: List of (temperature, RateConstant) tuples
        
        Returns:
            (E_a, A) - activation energy (kJ/mol) and pre-exponential factor
        """
        if len(rate_constants) < 2:
            raise ValueError("Need at least 2 temperatures to estimate E_a")
        
        temperatures = np.array([T for T, _ in rate_constants])
        k_values = np.array([rc.value for _, rc in rate_constants])
        
        # Arrhenius plot: ln(k) vs 1/T
        ln_k = np.log(k_values)
        inv_T = 1.0 / temperatures
        
        # Linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(inv_T, ln_k)
        
        # E_a = -R * slope
        E_a = -self.R * slope / 1000.0  # kJ/mol
        A = np.exp(intercept)
        
        logger.info(f"Activation energy: E_a = {E_a:.2f} kJ/mol (R^2 = {r_value**2:.3f})")
        
        return E_a, A
    
    def test_detailed_balance(self,
                             k_forward: RateConstant,
                             k_reverse: RateConstant,
                             K_eq: EquilibriumConstant) -> bool:
        """
        Test detailed balance: k_forward / k_reverse = K_eq
        
        Args:
            k_forward: Forward rate constant
            k_reverse: Reverse rate constant
            K_eq: Equilibrium constant
        
        Returns:
            True if detailed balance is satisfied (within tolerance)
        """
        ratio = k_forward.value / k_reverse.value
        
        # Allow 30% tolerance
        tolerance = 0.30
        lower = K_eq.K_eq * (1 - tolerance)
        upper = K_eq.K_eq * (1 + tolerance)
        
        is_valid = lower <= ratio <= upper
        
        if is_valid:
            logger.info(f"Detailed balance PASSED: k_f/k_r = {ratio:.3e}, K_eq = {K_eq.K_eq:.3e}")
        else:
            logger.warning(f"Detailed balance FAILED: k_f/k_r = {ratio:.3e}, K_eq = {K_eq.K_eq:.3e}")
        
        return is_valid


# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("=" * 70)
    print("REACTION KINETICS ANALYZER - EXAMPLE")
    print("=" * 70)
    
    analyzer = ReactionKineticsAnalyzer()
    
    # Example 1: First-order decay A -> B
    print("\nExample 1: First-order decay")
    print("-" * 70)
    
    times = np.linspace(0, 10, 100)
    k_true = 0.5  # 1/s
    A0 = 1.0
    A_conc = A0 * np.exp(-k_true * times)
    B_conc = A0 - A_conc
    
    data = KineticData(
        times=times,
        concentrations={'A': A_conc, 'B': B_conc},
        temperature=300
    )
    
    # Fit rate constant
    rate_const = analyzer.fit_rate_constant(data, 'A', order=1)
    print(f"Fitted: {rate_const}")
    print(f"True k = {k_true:.3f} 1/s")
    print(f"Error: {abs(rate_const.value - k_true) / k_true * 100:.1f}%")
    
    # Calculate half-life
    t_half = analyzer.calculate_half_life(rate_const)
    t_half_true = np.log(2) / k_true
    print(f"Half-life: {t_half:.3f} s (true: {t_half_true:.3f} s)")
    
    # Example 2: Equilibrium A ⇌ B
    print("\n\nExample 2: Equilibrium A <=> B")
    print("-" * 70)
    
    # Simulate approach to equilibrium
    times = np.linspace(0, 20, 200)
    K_eq_true = 2.0
    k_f = 0.3
    k_r = k_f / K_eq_true
    
    # Solve: d[A]/dt = -k_f[A] + k_r[B], [A] + [B] = A0
    A_eq = A0 / (1 + K_eq_true)
    B_eq = A0 * K_eq_true / (1 + K_eq_true)
    
    A_conc = A_eq + (A0 - A_eq) * np.exp(-(k_f + k_r) * times)
    B_conc = A0 - A_conc
    
    data_eq = KineticData(
        times=times,
        concentrations={'A': A_conc, 'B': B_conc},
        temperature=300
    )
    
    # Calculate equilibrium constant
    K_eq_calc = analyzer.calculate_equilibrium_constant(
        data_eq, reactants=['A'], products=['B']
    )
    print(f"Calculated: {K_eq_calc}")
    print(f"True K_eq = {K_eq_true:.3f}")
    print(f"Error: {abs(K_eq_calc.K_eq - K_eq_true) / K_eq_true * 100:.1f}%")
    
    # Test detailed balance
    print("\n\nExample 3: Detailed Balance Test")
    print("-" * 70)
    
    rate_f = analyzer.fit_rate_constant(data_eq, 'A', order=1)
    rate_r = RateConstant(value=k_r, error=0.0, units="1/s", order=1, r_squared=1.0)
    
    is_valid = analyzer.test_detailed_balance(rate_f, rate_r, K_eq_calc)
    
    print("\n" + "=" * 70)

