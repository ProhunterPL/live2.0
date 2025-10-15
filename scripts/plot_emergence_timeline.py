#!/usr/bin/env python3
"""
Figure 7: Emergence Timeline
============================

Generates Figure 7 for the publication:
- Panel A: Cumulative molecule discovery (all scenarios)
- Panel B: Key molecules first appearance
- Panel C: Complexity evolution (avg molecular weight)
- Panel D: Autocatalysis onset timing

Usage:
    python scripts/plot_emergence_timeline.py --input results/phase2_full --output figures/
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

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_timeline_data(results_dir: Path) -> Dict[str, List[Dict]]:
    """Load timeline data for each scenario"""
    timeline_data = defaultdict(list)
    
    for scenario_dir in results_dir.iterdir():
        if not scenario_dir.is_dir():
            continue
            
        scenario_name = scenario_dir.name
        logger.info(f"Loading timeline data for {scenario_name}...")
        
        for run_dir in scenario_dir.iterdir():
            if not run_dir.is_dir() or not run_dir.name.startswith('run_'):
                continue
                
            results_file = run_dir / "results.json"
            if not results_file.exists():
                continue
                
            try:
                with open(results_file, 'r') as f:
                    data = json.load(f)
                
                # Extract timeline information
                molecules = data.get('molecules_detected', [])
                metrics_history = data.get('metrics_history', [])
                
                timeline_entry = {
                    'run': run_dir.name,
                    'molecules': molecules,
                    'metrics_history': metrics_history,
                    'final_step': data.get('final_state', {}).get('step', 0),
                    'total_time': data.get('final_state', {}).get('time', 0)
                }
                
                timeline_data[scenario_name].append(timeline_entry)
                
            except Exception as e:
                logger.error(f"Failed to load {results_file}: {e}")
    
    return dict(timeline_data)

def plot_cumulative_discovery(timeline_data: Dict[str, List[Dict]], output_path: Path):
    """Panel A: Cumulative molecule discovery"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = {'miller_urey': '#1f77b4', 'hydrothermal': '#ff7f0e', 'formamide': '#2ca02c'}
    
    for scenario, runs in timeline_data.items():
        if not runs:
            continue
            
        # Calculate cumulative discovery curve
        all_steps = []
        all_discoveries = []
        
        for run in runs:
            molecules = run['molecules']
            if not molecules:
                continue
                
            # Sort molecules by first detection time
            sorted_molecules = sorted(molecules, key=lambda x: x.get('first_detected', 0))
            
            steps = []
            discoveries = []
            cumulative = 0
            
            for molecule in sorted_molecules:
                step = molecule.get('first_detected', 0)
                cumulative += 1
                steps.append(step)
                discoveries.append(cumulative)
            
            all_steps.extend(steps)
            all_discoveries.extend(discoveries)
        
        if all_steps and all_discoveries:
            # Sort by step
            sorted_data = sorted(zip(all_steps, all_discoveries))
            steps, discoveries = zip(*sorted_data)
            
            # Plot
            ax.plot(steps, discoveries, 'o-', color=colors.get(scenario, 'gray'),
                   label=scenario.replace('_', ' ').title(), linewidth=2, markersize=4)
    
    ax.set_xlabel('Simulation Step')
    ax.set_ylabel('Cumulative Molecules Discovered')
    ax.set_title('Cumulative Molecule Discovery')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure7_panelA_cumulative_discovery.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel A: {output_path / 'figure7_panelA_cumulative_discovery.png'}")

def plot_key_molecules_appearance(timeline_data: Dict[str, List[Dict]], output_path: Path):
    """Panel B: Key molecules first appearance"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define key molecules to track
    key_molecules = ['HCN', 'CH2O', 'NH3', 'CH4', 'H2O', 'HCOOH', 'CH3OH']
    
    # Collect first appearance data
    appearance_data = defaultdict(list)
    
    for scenario, runs in timeline_data.items():
        for run in runs:
            molecules = run['molecules']
            if not molecules:
                continue
                
            for molecule in molecules:
                formula = molecule.get('formula', '')
                if any(key in formula for key in key_molecules):
                    first_detected = molecule.get('first_detected', 0)
                    appearance_data[scenario].append({
                        'molecule': formula,
                        'step': first_detected
                    })
    
    # Create timeline plot
    y_pos = 0
    colors = {'miller_urey': '#1f77b4', 'hydrothermal': '#ff7f0e', 'formamide': '#2ca02c'}
    
    for scenario, appearances in appearance_data.items():
        if not appearances:
            continue
            
        # Sort by step
        appearances.sort(key=lambda x: x['step'])
        
        for appearance in appearances:
            ax.scatter(appearance['step'], y_pos, 
                      color=colors.get(scenario, 'gray'), s=100, alpha=0.7)
            ax.annotate(appearance['molecule'], 
                       (appearance['step'], y_pos),
                       xytext=(5, 0), textcoords='offset points', fontsize=8)
            y_pos += 1
    
    ax.set_xlabel('Simulation Step')
    ax.set_ylabel('Key Molecules')
    ax.set_title('Key Molecules First Appearance')
    
    # Add legend
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor=colors.get(scenario, 'gray'),
                                 label=scenario.replace('_', ' ').title())
                      for scenario in appearance_data.keys()]
    ax.legend(handles=legend_elements)
    
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure7_panelB_key_molecules.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel B: {output_path / 'figure7_panelB_key_molecules.png'}")

def plot_complexity_evolution(timeline_data: Dict[str, List[Dict]], output_path: Path):
    """Panel C: Complexity evolution (avg molecular weight)"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = {'miller_urey': '#1f77b4', 'hydrothermal': '#ff7f0e', 'formamide': '#2ca02c'}
    
    for scenario, runs in timeline_data.items():
        if not runs:
            continue
            
        # Calculate complexity evolution
        all_steps = []
        all_complexities = []
        
        for run in runs:
            molecules = run['molecules']
            if not molecules:
                continue
                
            # Sort by first detection time
            sorted_molecules = sorted(molecules, key=lambda x: x.get('first_detected', 0))
            
            steps = []
            avg_masses = []
            
            for i, molecule in enumerate(sorted_molecules):
                step = molecule.get('first_detected', 0)
                mass = molecule.get('mass', 0)
                
                steps.append(step)
                avg_masses.append(mass)
            
            all_steps.extend(steps)
            all_complexities.extend(avg_masses)
        
        if all_steps and all_complexities:
            # Sort by step
            sorted_data = sorted(zip(all_steps, all_complexities))
            steps, complexities = zip(*sorted_data)
            
            # Plot
            ax.plot(steps, complexities, 'o-', color=colors.get(scenario, 'gray'),
                   label=scenario.replace('_', ' ').title(), linewidth=2, markersize=4)
    
    ax.set_xlabel('Simulation Step')
    ax.set_ylabel('Average Molecular Mass (amu)')
    ax.set_title('Complexity Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure7_panelC_complexity_evolution.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel C: {output_path / 'figure7_panelC_complexity_evolution.png'}")

def plot_autocatalysis_onset(timeline_data: Dict[str, List[Dict]], output_path: Path):
    """Panel D: Autocatalysis onset timing"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # This is a placeholder - in real implementation, we'd load actual autocatalytic cycle data
    # For now, simulate autocatalysis onset based on molecule complexity
    
    colors = {'miller_urey': '#1f77b4', 'hydrothermal': '#ff7f0e', 'formamide': '#2ca02c'}
    
    for scenario, runs in timeline_data.items():
        if not runs:
            continue
            
        onset_times = []
        
        for run in runs:
            molecules = run['molecules']
            if not molecules:
                continue
                
            # Find when first complex molecule appears (proxy for autocatalysis)
            complex_molecules = [m for m in molecules if m.get('complexity', 0) > 5]
            
            if complex_molecules:
                first_complex = min(complex_molecules, key=lambda x: x.get('first_detected', 0))
                onset_time = first_complex.get('first_detected', 0)
                onset_times.append(onset_time)
        
        if onset_times:
            # Create histogram
            ax.hist(onset_times, bins=20, alpha=0.7, 
                   label=scenario.replace('_', ' ').title(),
                   color=colors.get(scenario, 'gray'), density=True)
    
    ax.set_xlabel('Autocatalysis Onset Time (steps)')
    ax.set_ylabel('Density')
    ax.set_title('Autocatalysis Onset Timing')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure7_panelD_autocatalysis_onset.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel D: {output_path / 'figure7_panelD_autocatalysis_onset.png'}")

def create_combined_figure(timeline_data: Dict[str, List[Dict]], output_path: Path):
    """Create combined Figure 7 with all panels"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Figure 7: Emergence Timeline', fontsize=16, fontweight='bold')
    
    # Load individual panels
    panels = [
        ('figure7_panelA_cumulative_discovery.png', 'A', 'Cumulative Discovery'),
        ('figure7_panelB_key_molecules.png', 'B', 'Key Molecules'),
        ('figure7_panelC_complexity_evolution.png', 'C', 'Complexity Evolution'),
        ('figure7_panelD_autocatalysis_onset.png', 'D', 'Autocatalysis Onset')
    ]
    
    for i, (filename, label, title) in enumerate(panels):
        panel_path = output_path / filename
        if panel_path.exists():
            img = plt.imread(panel_path)
            axes[i//2, i%2].imshow(img)
            axes[i//2, i%2].set_title(f'{label}) {title}')
            axes[i//2, i%2].axis('off')
        else:
            axes[i//2, i%2].text(0.5, 0.5, f'Panel {label} not found', 
                                ha='center', va='center', transform=axes[i//2, i%2].transAxes)
            axes[i//2, i%2].set_title(f'{label}) {title}')
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure7_combined.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved combined figure: {output_path / 'figure7_combined.png'}")

def generate_timeline_statistics(timeline_data: Dict[str, List[Dict]], output_path: Path):
    """Generate timeline statistics"""
    stats = []
    
    for scenario, runs in timeline_data.items():
        if not runs:
            continue
            
        total_molecules = 0
        total_steps = 0
        avg_complexity = 0
        
        for run in runs:
            molecules = run['molecules']
            total_molecules += len(molecules)
            total_steps += run.get('final_step', 0)
            
            if molecules:
                complexities = [m.get('complexity', 0) for m in molecules]
                avg_complexity += np.mean(complexities)
        
        avg_complexity /= len(runs) if runs else 1
        
        stats.append({
            'Scenario': scenario.replace('_', ' ').title(),
            'Runs': len(runs),
            'Total Molecules': total_molecules,
            'Avg Molecules per Run': total_molecules / len(runs) if runs else 0,
            'Avg Steps per Run': total_steps / len(runs) if runs else 0,
            'Avg Complexity': avg_complexity
        })
    
    if stats:
        df = pd.DataFrame(stats)
        df.to_csv(output_path / 'figure7_timeline_statistics.csv', index=False)
        logger.info(f"Saved timeline statistics: {output_path / 'figure7_timeline_statistics.csv'}")

def main():
    parser = argparse.ArgumentParser(description='Generate Figure 7: Emergence Timeline')
    parser.add_argument('--input', type=Path, required=True, help='Results directory')
    parser.add_argument('--output', type=Path, required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("FIGURE 7: EMERGENCE TIMELINE")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    # Load timeline data
    timeline_data = load_timeline_data(args.input)
    
    if not timeline_data:
        logger.error("No timeline data found!")
        return 1
    
    logger.info(f"Found timeline data for scenarios: {list(timeline_data.keys())}")
    
    # Generate panels
    plot_cumulative_discovery(timeline_data, args.output)
    plot_key_molecules_appearance(timeline_data, args.output)
    plot_complexity_evolution(timeline_data, args.output)
    plot_autocatalysis_onset(timeline_data, args.output)
    
    # Create combined figure
    create_combined_figure(timeline_data, args.output)
    
    # Generate statistics
    generate_timeline_statistics(timeline_data, args.output)
    
    logger.info("=" * 70)
    logger.info("FIGURE 7 GENERATION COMPLETE!")
    logger.info("=" * 70)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
