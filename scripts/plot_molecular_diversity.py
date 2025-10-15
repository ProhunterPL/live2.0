#!/usr/bin/env python3
"""
Figure 3: Molecular Diversity Across Scenarios
==============================================

Generates Figure 3 for the publication:
- Panel A: Species accumulation curves (3 scenarios)
- Panel B: Size distribution (violin plots)
- Panel C: Shannon entropy evolution
- Panel D: Venn diagram (shared vs unique molecules)

Usage:
    python scripts/plot_molecular_diversity.py --input results/phase2_full --output figures/
"""

import sys
import argparse
import logging
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict
import pandas as pd
from matplotlib_venn import venn3

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_scenario_data(results_dir: Path) -> Dict[str, Dict]:
    """Load all simulation results for each scenario"""
    scenario_data = defaultdict(list)
    
    for scenario_dir in results_dir.iterdir():
        if not scenario_dir.is_dir():
            continue
            
        scenario_name = scenario_dir.name
        logger.info(f"Loading {scenario_name} data...")
        
        for run_dir in scenario_dir.iterdir():
            if not run_dir.is_dir() or not run_dir.name.startswith('run_'):
                continue
                
            results_file = run_dir / "results.json"
            if not results_file.exists():
                logger.warning(f"No results.json in {run_dir}")
                continue
                
            try:
                with open(results_file, 'r') as f:
                    data = json.load(f)
                
                molecules = data.get('molecules_detected', [])
                scenario_data[scenario_name].append({
                    'run': run_dir.name,
                    'molecules': molecules,
                    'total_instances': sum(m.get('count', 0) for m in molecules),
                    'unique_molecules': len(molecules)
                })
                
            except Exception as e:
                logger.error(f"Failed to load {results_file}: {e}")
    
    return dict(scenario_data)

def plot_species_accumulation(scenario_data: Dict[str, Dict], output_path: Path):
    """Panel A: Species accumulation curves"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = {'miller_urey': '#1f77b4', 'hydrothermal': '#ff7f0e', 'formamide': '#2ca02c'}
    
    for scenario, runs in scenario_data.items():
        if not runs:
            continue
            
        # Calculate accumulation curve
        all_molecules = set()
        accumulation = []
        
        for run in runs:
            run_molecules = set(m['formula'] for m in run['molecules'])
            all_molecules.update(run_molecules)
            accumulation.append(len(all_molecules))
        
        # Plot
        x = range(1, len(accumulation) + 1)
        ax.plot(x, accumulation, 'o-', color=colors.get(scenario, 'gray'), 
                label=scenario.replace('_', ' ').title(), linewidth=2, markersize=6)
    
    ax.set_xlabel('Number of Runs')
    ax.set_ylabel('Cumulative Unique Molecules')
    ax.set_title('Species Accumulation Curves')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure3_panelA_species_accumulation.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel A: {output_path / 'figure3_panelA_species_accumulation.png'}")

def plot_size_distribution(scenario_data: Dict[str, Dict], output_path: Path):
    """Panel B: Size distribution (violin plots)"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Collect data
    data_for_plot = []
    for scenario, runs in scenario_data.items():
        for run in runs:
            for molecule in run['molecules']:
                data_for_plot.append({
                    'scenario': scenario.replace('_', ' ').title(),
                    'mass': molecule.get('mass', 0),
                    'complexity': molecule.get('complexity', 0)
                })
    
    if not data_for_plot:
        logger.warning("No data for size distribution plot")
        return
    
    df = pd.DataFrame(data_for_plot)
    
    # Create violin plot
    sns.violinplot(data=df, x='scenario', y='mass', ax=ax)
    ax.set_xlabel('Scenario')
    ax.set_ylabel('Molecular Mass (amu)')
    ax.set_title('Molecular Size Distribution')
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure3_panelB_size_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel B: {output_path / 'figure3_panelB_size_distribution.png'}")

def plot_shannon_entropy(scenario_data: Dict[str, Dict], output_path: Path):
    """Panel C: Shannon entropy evolution"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = {'miller_urey': '#1f77b4', 'hydrothermal': '#ff7f0e', 'formamide': '#2ca02c'}
    
    for scenario, runs in scenario_data.items():
        if not runs:
            continue
            
        # Calculate Shannon entropy for each run
        entropies = []
        for run in runs:
            molecules = run['molecules']
            if not molecules:
                entropies.append(0)
                continue
                
            # Calculate relative abundances
            total_count = sum(m.get('count', 0) for m in molecules)
            if total_count == 0:
                entropies.append(0)
                continue
                
            abundances = [m.get('count', 0) / total_count for m in molecules]
            abundances = [p for p in abundances if p > 0]  # Remove zeros
            
            # Shannon entropy: H = -Σ(p * log(p))
            entropy = -sum(p * np.log(p) for p in abundances)
            entropies.append(entropy)
        
        # Plot
        x = range(1, len(entropies) + 1)
        ax.plot(x, entropies, 'o-', color=colors.get(scenario, 'gray'),
                label=scenario.replace('_', ' ').title(), linewidth=2, markersize=6)
    
    ax.set_xlabel('Run Number')
    ax.set_ylabel('Shannon Entropy')
    ax.set_title('Molecular Diversity (Shannon Entropy)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure3_panelC_shannon_entropy.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel C: {output_path / 'figure3_panelC_shannon_entropy.png'}")

def plot_venn_diagram(scenario_data: Dict[str, Dict], output_path: Path):
    """Panel D: Venn diagram (shared vs unique molecules)"""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Collect molecules for each scenario
    scenario_molecules = {}
    for scenario, runs in scenario_data.items():
        all_molecules = set()
        for run in runs:
            run_molecules = set(m['formula'] for m in run['molecules'])
            all_molecules.update(run_molecules)
        scenario_molecules[scenario] = all_molecules
    
    # Create Venn diagram
    if len(scenario_molecules) == 3:
        miller = scenario_molecules.get('miller_urey', set())
        hydro = scenario_molecules.get('hydrothermal', set())
        forma = scenario_molecules.get('formamide', set())
        
        venn3([miller, hydro, forma], 
              ['Miller-Urey', 'Hydrothermal', 'Formamide'],
              ax=ax)
        
        ax.set_title('Shared vs Unique Molecules')
        
    else:
        ax.text(0.5, 0.5, 'Need 3 scenarios for Venn diagram', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Venn Diagram (Insufficient Data)')
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure3_panelD_venn_diagram.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel D: {output_path / 'figure3_panelD_venn_diagram.png'}")

def create_combined_figure(output_path: Path):
    """Create combined Figure 3 with all panels"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Figure 3: Molecular Diversity Across Scenarios', fontsize=16, fontweight='bold')
    
    # Load individual panels
    panels = [
        ('figure3_panelA_species_accumulation.png', 'A', 'Species Accumulation'),
        ('figure3_panelB_size_distribution.png', 'B', 'Size Distribution'),
        ('figure3_panelC_shannon_entropy.png', 'C', 'Shannon Entropy'),
        ('figure3_panelD_venn_diagram.png', 'D', 'Venn Diagram')
    ]
    
    for i, (filename, label, title) in enumerate(panels):
        panel_path = output_path / filename
        if panel_path.exists():
            # Load and display the panel
            img = plt.imread(panel_path)
            row, col = i // 2, i % 2
            axes[row, col].imshow(img)
            axes[row, col].set_title(f'{label}) {title}')
            axes[row, col].axis('off')
        else:
            axes[i // 2, i % 2].text(0.5, 0.5, f'Panel {label} not found', 
                                    ha='center', va='center', transform=axes[i // 2, i % 2].transAxes)
            axes[i // 2, i % 2].set_title(f'{label}) {title}')
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure3_combined.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved combined figure: {output_path / 'figure3_combined.png'}")

def main():
    parser = argparse.ArgumentParser(description='Generate Figure 3: Molecular Diversity')
    parser.add_argument('--input', type=Path, required=True, help='Results directory')
    parser.add_argument('--output', type=Path, required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("FIGURE 3: MOLECULAR DIVERSITY ACROSS SCENARIOS")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    # Load data
    scenario_data = load_scenario_data(args.input)
    
    if not scenario_data:
        logger.error("No scenario data found!")
        return 1
    
    logger.info(f"Found data for scenarios: {list(scenario_data.keys())}")
    
    # Generate panels
    plot_species_accumulation(scenario_data, args.output)
    plot_size_distribution(scenario_data, args.output)
    plot_shannon_entropy(scenario_data, args.output)
    plot_venn_diagram(scenario_data, args.output)
    
    # Create combined figure
    create_combined_figure(args.output)
    
    logger.info("=" * 70)
    logger.info("FIGURE 3 GENERATION COMPLETE!")
    logger.info("=" * 70)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
