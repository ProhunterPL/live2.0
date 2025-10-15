#!/usr/bin/env python3
"""
Figure 5: Autocatalytic Cycles
==============================

Generates Figure 5 for the publication:
- Panel A: Direct autocatalysis examples (3 best)
- Panel B: Indirect cycle example (Miller-Urey)
- Panel C: Hypercycle (if detected)
- Panel D: Cycle frequency by scenario

Usage:
    python scripts/plot_autocatalytic_cycles.py --input analysis/autocatalytic_cycles --output figures/
"""

import sys
import argparse
import logging
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict
import pandas as pd

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_cycle_data(cycles_dir: Path) -> Dict[str, Dict]:
    """Load autocatalytic cycle data for each scenario"""
    cycle_data = {}
    
    for scenario_dir in cycles_dir.iterdir():
        if not scenario_dir.is_dir():
            continue
            
        scenario_name = scenario_dir.name
        logger.info(f"Loading cycles for {scenario_name}...")
        
        cycles_file = scenario_dir / "autocatalytic_cycles.json"
        if not cycles_file.exists():
            logger.warning(f"No cycles file found for {scenario_name}")
            continue
            
        try:
            with open(cycles_file, 'r') as f:
                data = json.load(f)
            
            cycle_data[scenario_name] = {
                'direct_cycles': data.get('direct_autocatalysis', []),
                'indirect_cycles': data.get('indirect_cycles', []),
                'hypercycles': data.get('hypercycles', []),
                'raf_sets': data.get('raf_sets', [])
            }
            
            logger.info(f"  {scenario_name}: {len(data.get('direct_autocatalysis', []))} direct, "
                       f"{len(data.get('indirect_cycles', []))} indirect, "
                       f"{len(data.get('hypercycles', []))} hypercycles")
            
        except Exception as e:
            logger.error(f"Failed to load {cycles_file}: {e}")
    
    return cycle_data

def plot_direct_autocatalysis(cycle_data: Dict[str, Dict], output_path: Path):
    """Panel A: Direct autocatalysis examples"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Direct Autocatalysis Examples', fontsize=14, fontweight='bold')
    
    # Find best direct cycles across all scenarios
    all_direct_cycles = []
    for scenario, data in cycle_data.items():
        for cycle in data.get('direct_autocatalysis', []):
            cycle['scenario'] = scenario
            all_direct_cycles.append(cycle)
    
    # Sort by confidence/strength
    all_direct_cycles.sort(key=lambda x: x.get('confidence', 0), reverse=True)
    
    # Plot top 3
    for i, cycle in enumerate(all_direct_cycles[:3]):
        if i >= 3:
            break
            
        ax = axes[i]
        
        # Create simple cycle visualization
        G = nx.DiGraph()
        
        # Add nodes (molecules)
        reactants = cycle.get('reactants', [])
        products = cycle.get('products', [])
        
        for mol in reactants + products:
            G.add_node(mol, type='reactant' if mol in reactants else 'product')
        
        # Add edges (reactions)
        for reactant in reactants:
            for product in products:
                G.add_edge(reactant, product)
        
        if G.number_of_nodes() > 0:
            pos = nx.spring_layout(G, k=1, iterations=50)
            
            # Color nodes by type
            node_colors = ['lightcoral' if G.nodes[node]['type'] == 'reactant' else 'lightblue' 
                          for node in G.nodes()]
            
            nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                                 node_size=500, alpha=0.8, ax=ax)
            nx.draw_networkx_edges(G, pos, edge_color='gray', 
                                 arrows=True, arrowsize=20, ax=ax)
            nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
        
        ax.set_title(f"{cycle.get('scenario', 'Unknown').replace('_', ' ').title()}\n"
                    f"Confidence: {cycle.get('confidence', 0):.2f}")
        ax.axis('off')
    
    # Hide unused subplots
    for i in range(len(all_direct_cycles), 3):
        axes[i].text(0.5, 0.5, 'No data', ha='center', va='center', 
                    transform=axes[i].transAxes)
        axes[i].set_title('No Cycle')
        axes[i].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure5_panelA_direct_autocatalysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel A: {output_path / 'figure5_panelA_direct_autocatalysis.png'}")

def plot_indirect_cycle(cycle_data: Dict[str, Dict], output_path: Path):
    """Panel B: Indirect cycle example (Miller-Urey)"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Find best indirect cycle from Miller-Urey
    miller_cycles = cycle_data.get('miller_urey', {}).get('indirect_cycles', [])
    
    if not miller_cycles:
        ax.text(0.5, 0.5, 'No indirect cycles found in Miller-Urey', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Indirect Cycle Example (Miller-Urey)')
        plt.tight_layout()
        plt.savefig(output_path / 'figure5_panelB_indirect_cycle.png', dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Use the longest cycle
    best_cycle = max(miller_cycles, key=lambda x: len(x.get('cycle', [])))
    
    # Create cycle visualization
    G = nx.DiGraph()
    cycle_molecules = best_cycle.get('cycle', [])
    
    for i, mol in enumerate(cycle_molecules):
        G.add_node(mol, step=i)
    
    # Add edges in cycle order
    for i in range(len(cycle_molecules)):
        next_i = (i + 1) % len(cycle_molecules)
        G.add_edge(cycle_molecules[i], cycle_molecules[next_i])
    
    if G.number_of_nodes() > 0:
        # Circular layout
        pos = nx.circular_layout(G)
        
        nx.draw_networkx_nodes(G, pos, node_color='lightgreen', 
                              node_size=800, alpha=0.8, ax=ax)
        nx.draw_networkx_edges(G, pos, edge_color='darkgreen', 
                              arrows=True, arrowsize=20, width=2, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=10, ax=ax)
        
        # Add cycle information
        ax.text(0.02, 0.98, f"Cycle Length: {len(cycle_molecules)}\n"
                           f"Amplification: {best_cycle.get('amplification_factor', 'N/A')}", 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.set_title('Indirect Cycle Example (Miller-Urey)')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure5_panelB_indirect_cycle.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel B: {output_path / 'figure5_panelB_indirect_cycle.png'}")

def plot_hypercycle(cycle_data: Dict[str, Dict], output_path: Path):
    """Panel C: Hypercycle (if detected)"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Find any hypercycle
    hypercycle = None
    scenario_name = None
    
    for scenario, data in cycle_data.items():
        hypercycles = data.get('hypercycles', [])
        if hypercycles:
            hypercycle = hypercycles[0]  # Use first one
            scenario_name = scenario
            break
    
    if not hypercycle:
        ax.text(0.5, 0.5, 'No hypercycles detected', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Hypercycle Detection')
        plt.tight_layout()
        plt.savefig(output_path / 'figure5_panelC_hypercycle.png', dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Create hypercycle visualization
    G = nx.DiGraph()
    
    # Add nodes (catalytic species)
    species = hypercycle.get('species', [])
    for species_name in species:
        G.add_node(species_name)
    
    # Add edges (catalytic relationships)
    relationships = hypercycle.get('relationships', [])
    for rel in relationships:
        if len(rel) >= 2:
            G.add_edge(rel[0], rel[1])
    
    if G.number_of_nodes() > 0:
        pos = nx.spring_layout(G, k=2, iterations=100)
        
        nx.draw_networkx_nodes(G, pos, node_color='gold', 
                              node_size=1000, alpha=0.8, ax=ax)
        nx.draw_networkx_edges(G, pos, edge_color='orange', 
                              arrows=True, arrowsize=20, width=2, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=10, ax=ax)
        
        # Add hypercycle info
        ax.text(0.02, 0.98, f"Scenario: {scenario_name.replace('_', ' ').title()}\n"
                           f"Species: {len(species)}\n"
                           f"Relationships: {len(relationships)}", 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    ax.set_title('Hypercycle Detection')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure5_panelC_hypercycle.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel C: {output_path / 'figure5_panelC_hypercycle.png'}")

def plot_cycle_frequency(cycle_data: Dict[str, Dict], output_path: Path):
    """Panel D: Cycle frequency by scenario"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Collect data
    scenarios = []
    direct_counts = []
    indirect_counts = []
    hypercycle_counts = []
    
    for scenario, data in cycle_data.items():
        scenarios.append(scenario.replace('_', ' ').title())
        direct_counts.append(len(data.get('direct_autocatalysis', [])))
        indirect_counts.append(len(data.get('indirect_cycles', [])))
        hypercycle_counts.append(len(data.get('hypercycles', [])))
    
    if not scenarios:
        ax.text(0.5, 0.5, 'No cycle data available', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Cycle Frequency by Scenario')
        plt.tight_layout()
        plt.savefig(output_path / 'figure5_panelD_cycle_frequency.png', dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Create grouped bar chart
    x = np.arange(len(scenarios))
    width = 0.25
    
    ax.bar(x - width, direct_counts, width, label='Direct Autocatalysis', color='lightcoral')
    ax.bar(x, indirect_counts, width, label='Indirect Cycles', color='lightblue')
    ax.bar(x + width, hypercycle_counts, width, label='Hypercycles', color='gold')
    
    ax.set_xlabel('Scenario')
    ax.set_ylabel('Number of Cycles')
    ax.set_title('Cycle Frequency by Scenario')
    ax.set_xticks(x)
    ax.set_xticklabels(scenarios)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'figure5_panelD_cycle_frequency.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Panel D: {output_path / 'figure5_panelD_cycle_frequency.png'}")

def create_combined_figure(cycle_data: Dict[str, Dict], output_path: Path):
    """Create combined Figure 5 with all panels"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Figure 5: Autocatalytic Cycles', fontsize=16, fontweight='bold')
    
    # Load individual panels
    panels = [
        ('figure5_panelA_direct_autocatalysis.png', 'A', 'Direct Autocatalysis'),
        ('figure5_panelB_indirect_cycle.png', 'B', 'Indirect Cycle'),
        ('figure5_panelC_hypercycle.png', 'C', 'Hypercycle'),
        ('figure5_panelD_cycle_frequency.png', 'D', 'Cycle Frequency')
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
    plt.savefig(output_path / 'figure5_combined.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved combined figure: {output_path / 'figure5_combined.png'}")

def main():
    parser = argparse.ArgumentParser(description='Generate Figure 5: Autocatalytic Cycles')
    parser.add_argument('--input', type=Path, required=True, help='Autocatalytic cycles directory')
    parser.add_argument('--output', type=Path, required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("FIGURE 5: AUTOCATALYTIC CYCLES")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    # Load cycle data
    cycle_data = load_cycle_data(args.input)
    
    if not cycle_data:
        logger.error("No cycle data found!")
        return 1
    
    logger.info(f"Found cycle data for scenarios: {list(cycle_data.keys())}")
    
    # Generate panels
    plot_direct_autocatalysis(cycle_data, args.output)
    plot_indirect_cycle(cycle_data, args.output)
    plot_hypercycle(cycle_data, args.output)
    plot_cycle_frequency(cycle_data, args.output)
    
    # Create combined figure
    create_combined_figure(cycle_data, args.output)
    
    logger.info("=" * 70)
    logger.info("FIGURE 5 GENERATION COMPLETE!")
    logger.info("=" * 70)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
