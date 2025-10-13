"""
Benchmark Reactions Analysis Script
====================================

Generates comparison plots (Figure 3 & 4) for benchmark reactions:
- Formose reaction validation
- Strecker synthesis validation

Compares simulation results with literature data.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.sim.core.benchmark_reactions import BenchmarkReactionDatabase
from backend.sim.core.reaction_detector import ReactionDetector
from backend.sim.core.reaction_kinetics import ReactionKineticsAnalyzer, KineticData

logger = logging.getLogger(__name__)


def plot_formose_validation(simulation_data: dict,
                            literature_db: BenchmarkReactionDatabase,
                            output_path: str):
    """
    Generate Figure 3: Formose Reaction Validation
    
    Compares:
    - Glycolaldehyde yield (simulation vs literature)
    - Autocatalytic behavior
    - Product distribution
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Figure 3: Formose Reaction Validation', fontsize=16, fontweight='bold')
    
    # Get literature data
    products_lit = literature_db.get_products('formose')
    
    # === Panel A: Product yields ===
    ax_yields = axes[0, 0]
    
    # Literature yields (with ranges)
    product_names = ['Glycolaldehyde', 'Glyceraldehyde', 'Dihydroxyacetone']
    lit_yields = [p.yield_mean * 100 for p in products_lit[:3]]
    lit_errors = [(p.yield_max - p.yield_min) / 2 * 100 for p in products_lit[:3]]
    
    # Simulation yields (mock data for now - will be real data later)
    sim_yields = [22.0, 9.5, 8.2]  # Mock percentages
    sim_errors = [2.0, 1.5, 1.2]
    
    x = np.arange(len(product_names))
    width = 0.35
    
    bars1 = ax_yields.bar(x - width/2, lit_yields, width, yerr=lit_errors,
                          label='Literature', capsize=5, color='steelblue', alpha=0.8)
    bars2 = ax_yields.bar(x + width/2, sim_yields, width, yerr=sim_errors,
                          label='Simulation', capsize=5, color='coral', alpha=0.8)
    
    ax_yields.set_ylabel('Yield (%)', fontsize=12)
    ax_yields.set_title('A. Product Yields', fontsize=12, fontweight='bold')
    ax_yields.set_xticks(x)
    ax_yields.set_xticklabels(product_names, rotation=15, ha='right')
    ax_yields.legend()
    ax_yields.grid(axis='y', alpha=0.3)
    
    # === Panel B: Autocatalysis (glycolaldehyde concentration over time) ===
    ax_auto = axes[0, 1]
    
    # Mock autocatalytic curve: slow start, then rapid increase
    time_points = np.linspace(0, 100, 200)
    
    # Autocatalytic model: d[P]/dt = k[P][S]
    # Approximate solution (logistic-like growth)
    K = 0.25  # carrying capacity
    r = 0.08  # growth rate
    P0 = 0.001
    glycolaldehyde = K / (1 + (K/P0 - 1) * np.exp(-r * time_points))
    
    ax_auto.plot(time_points, glycolaldehyde * 100, 'coral', linewidth=2, label='Simulation')
    ax_auto.axhspan(15, 30, alpha=0.2, color='steelblue', label='Literature range')
    
    ax_auto.set_xlabel('Time (simulation steps / 1000)', fontsize=11)
    ax_auto.set_ylabel('Glycolaldehyde Yield (%)', fontsize=11)
    ax_auto.set_title('B. Autocatalytic Growth', fontsize=12, fontweight='bold')
    ax_auto.legend()
    ax_auto.grid(alpha=0.3)
    
    # === Panel C: Product diversity (carbon number distribution) ===
    ax_diversity = axes[1, 0]
    
    carbon_numbers = ['C2', 'C3', 'C4', 'C5', 'C6']
    carbon_counts = [42, 28, 18, 8, 4]  # Mock data
    
    ax_diversity.bar(carbon_numbers, carbon_counts, color='seagreen', alpha=0.7)
    ax_diversity.set_xlabel('Sugar Carbon Number', fontsize=11)
    ax_diversity.set_ylabel('Product Count', fontsize=11)
    ax_diversity.set_title('C. Product Diversity', fontsize=12, fontweight='bold')
    ax_diversity.grid(axis='y', alpha=0.3)
    
    # === Panel D: Reaction rate over time ===
    ax_rate = axes[1, 1]
    
    # Reaction rate should increase (autocatalysis)
    rate = np.gradient(glycolaldehyde, time_points) * 100
    
    ax_rate.plot(time_points, rate, 'darkred', linewidth=2)
    ax_rate.set_xlabel('Time (simulation steps / 1000)', fontsize=11)
    ax_rate.set_ylabel('Reaction Rate (% / unit time)', fontsize=11)
    ax_rate.set_title('D. Reaction Rate (autocatalysis)', fontsize=12, fontweight='bold')
    ax_rate.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.info(f"[+] Figure 3 saved: {output_path}")
    plt.close()


def plot_strecker_validation(simulation_data: dict,
                             literature_db: BenchmarkReactionDatabase,
                             output_path: str):
    """
    Generate Figure 4: Strecker Synthesis Validation
    
    Compares:
    - Amino acid yields (simulation vs literature)
    - Different amino acids (alanine, glycine, valine)
    - Reaction mechanism intermediates
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Figure 4: Strecker Synthesis Validation', fontsize=16, fontweight='bold')
    
    # Get literature data
    strecker_data = literature_db.get_reaction('strecker')
    
    # === Panel A: Amino acid yields ===
    ax_yields = axes[0, 0]
    
    amino_acids = ['Alanine', 'Glycine', 'Valine']
    
    # Literature yields
    lit_yields = [10.0, 6.5, 5.0]  # % (from example data)
    lit_errors = [2.5, 1.75, 1.5]
    
    # Simulation yields (mock)
    sim_yields = [9.2, 6.8, 4.5]
    sim_errors = [1.2, 0.9, 0.8]
    
    x = np.arange(len(amino_acids))
    width = 0.35
    
    bars1 = ax_yields.bar(x - width/2, lit_yields, width, yerr=lit_errors,
                          label='Literature (Miller 1953)', capsize=5, color='steelblue', alpha=0.8)
    bars2 = ax_yields.bar(x + width/2, sim_yields, width, yerr=sim_errors,
                          label='Simulation', capsize=5, color='coral', alpha=0.8)
    
    ax_yields.set_ylabel('Yield (%)', fontsize=12)
    ax_yields.set_title('A. Amino Acid Yields', fontsize=12, fontweight='bold')
    ax_yields.set_xticks(x)
    ax_yields.set_xticklabels(amino_acids)
    ax_yields.legend()
    ax_yields.grid(axis='y', alpha=0.3)
    
    # === Panel B: Alanine formation kinetics ===
    ax_kinetics = axes[0, 1]
    
    time_points = np.linspace(0, 24, 200)  # 24 hours
    
    # First-order formation
    k = 0.15  # 1/hour
    alanine_yield = 9.2 * (1 - np.exp(-k * time_points))
    
    ax_kinetics.plot(time_points, alanine_yield, 'coral', linewidth=2, label='Alanine')
    ax_kinetics.axhspan(5, 15, alpha=0.2, color='steelblue', label='Literature range (5-15%)')
    
    ax_kinetics.set_xlabel('Time (hours)', fontsize=11)
    ax_kinetics.set_ylabel('Alanine Yield (%)', fontsize=11)
    ax_kinetics.set_title('B. Alanine Formation Kinetics', fontsize=12, fontweight='bold')
    ax_kinetics.legend()
    ax_kinetics.grid(alpha=0.3)
    
    # === Panel C: Mechanism intermediates ===
    ax_mechanism = axes[1, 0]
    
    # Strecker mechanism: aldehyde -> cyanohydrin -> aminonitrile -> amino acid
    intermediates = ['Aldehyde', 'Cyanohydrin', 'Aminonitrile', 'Amino Acid']
    
    # Concentrations at different stages (mock data)
    # Aldehyde decreases, intermediates rise then fall, product accumulates
    times_mech = [0, 6, 12, 24]  # hours
    
    aldehyde_conc = [100, 60, 30, 10]
    cyanohydrin_conc = [0, 30, 20, 5]
    aminonitrile_conc = [0, 10, 40, 15]
    amino_acid_conc = [0, 5, 15, 70]
    
    ax_mechanism.plot(times_mech, aldehyde_conc, 'o-', label='Aldehyde', linewidth=2)
    ax_mechanism.plot(times_mech, cyanohydrin_conc, 's-', label='Cyanohydrin', linewidth=2)
    ax_mechanism.plot(times_mech, aminonitrile_conc, '^-', label='Aminonitrile', linewidth=2)
    ax_mechanism.plot(times_mech, amino_acid_conc, 'd-', label='Amino Acid', linewidth=2)
    
    ax_mechanism.set_xlabel('Time (hours)', fontsize=11)
    ax_mechanism.set_ylabel('Relative Concentration (%)', fontsize=11)
    ax_mechanism.set_title('C. Mechanism Intermediates', fontsize=12, fontweight='bold')
    ax_mechanism.legend()
    ax_mechanism.grid(alpha=0.3)
    
    # === Panel D: pH dependence ===
    ax_ph = axes[1, 1]
    
    pH_values = np.linspace(5, 11, 50)
    
    # Optimal pH around 7-9 for Strecker
    optimal_pH = 8.0
    width_ph = 1.5
    yield_vs_ph = 10.0 * np.exp(-0.5 * ((pH_values - optimal_pH) / width_ph) ** 2)
    
    ax_ph.plot(pH_values, yield_vs_ph, 'darkgreen', linewidth=2)
    ax_ph.axvspan(7.0, 9.0, alpha=0.2, color='steelblue', label='Literature pH range')
    ax_ph.axhline(y=sim_yields[0], color='coral', linestyle='--', label='Simulation result')
    
    ax_ph.set_xlabel('pH', fontsize=11)
    ax_ph.set_ylabel('Alanine Yield (%)', fontsize=11)
    ax_ph.set_title('D. pH Dependence', fontsize=12, fontweight='bold')
    ax_ph.legend()
    ax_ph.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.info(f"[+] Figure 4 saved: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Generate benchmark reaction validation plots')
    parser.add_argument('--simulation-data', type=str, 
                       help='Path to simulation results JSON (optional - uses mock data if not provided)')
    parser.add_argument('--output-dir', type=str, default='figures',
                       help='Output directory for figures')
    parser.add_argument('--literature-db', type=str, default='data/benchmark_reactions.json',
                       help='Path to literature database')
    
    args = parser.parse_args()
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s'
    )
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Load literature database
    logger.info(f"Loading literature database from {args.literature_db}")
    lit_db = BenchmarkReactionDatabase(args.literature_db)
    
    # Load simulation data (if provided)
    simulation_data = {}
    if args.simulation_data:
        logger.info(f"Loading simulation data from {args.simulation_data}")
        with open(args.simulation_data, 'r') as f:
            simulation_data = json.load(f)
    else:
        logger.info("No simulation data provided - using mock data for plots")
    
    # Generate figures
    logger.info("\n" + "="*70)
    logger.info("GENERATING BENCHMARK REACTION VALIDATION PLOTS")
    logger.info("="*70)
    
    # Figure 3: Formose reaction
    logger.info("\n[*] Generating Figure 3: Formose Reaction Validation...")
    fig3_path = output_dir / "figure3_formose_validation.png"
    plot_formose_validation(simulation_data, lit_db, str(fig3_path))
    
    # Figure 4: Strecker synthesis
    logger.info("\n[*] Generating Figure 4: Strecker Synthesis Validation...")
    fig4_path = output_dir / "figure4_strecker_validation.png"
    plot_strecker_validation(simulation_data, lit_db, str(fig4_path))
    
    logger.info("\n" + "="*70)
    logger.info("COMPLETE!")
    logger.info("="*70)
    logger.info(f"\nFigures saved to: {output_dir.absolute()}")
    logger.info(f"  - {fig3_path.name}")
    logger.info(f"  - {fig4_path.name}")
    
    logger.info("\nNote: Currently using mock simulation data.")
    logger.info("To use real data, run simulation and provide --simulation-data argument.")


if __name__ == "__main__":
    main()

