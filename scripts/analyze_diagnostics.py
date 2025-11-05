#!/usr/bin/env python3
"""
Diagnostics Analysis Script for Live 2.0

Analyzes CSV diagnostics files and generates plots showing:
- Bond dynamics
- Cluster evolution
- Event rates
- Energy conservation
- Phase transitions

Usage:
    python analyze_diagnostics.py diagnostics/
    python analyze_diagnostics.py diagnostics/metrics_20250101_120000.csv
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter

def load_diagnostics(path: Path):
    """Load diagnostics CSV files"""
    if path.is_file():
        # Single file provided
        metrics_file = path
        base_dir = path.parent
        timestamp = path.stem.split('_', 1)[1]
    else:
        # Directory provided, find most recent
        metrics_files = list(path.glob('metrics_*.csv'))
        if not metrics_files:
            raise FileNotFoundError(f"No metrics_*.csv files found in {path}")
        metrics_file = sorted(metrics_files)[-1]
        base_dir = path
        timestamp = metrics_file.stem.split('_', 1)[1]
    
    print(f"Loading metrics from: {metrics_file}")
    
    data = {
        'metrics': pd.read_csv(metrics_file),
        'bond_types': None,
        'cluster_dist': None,
        'lifetimes': None
    }
    
    # Try to load other files
    bond_types_file = base_dir / f'bond_types_{timestamp}.csv'
    if bond_types_file.exists():
        data['bond_types'] = pd.read_csv(bond_types_file)
        print(f"Loaded bond types from: {bond_types_file}")
    
    cluster_dist_file = base_dir / f'cluster_dist_{timestamp}.csv'
    if cluster_dist_file.exists():
        data['cluster_dist'] = pd.read_csv(cluster_dist_file)
        print(f"Loaded cluster distribution from: {cluster_dist_file}")
    
    lifetimes_file = base_dir / f'bond_lifetimes_{timestamp}.csv'
    if lifetimes_file.exists():
        data['lifetimes'] = pd.read_csv(lifetimes_file)
        print(f"Loaded bond lifetimes from: {lifetimes_file}")
    
    return data

def plot_bond_dynamics(df: pd.DataFrame, output_dir: Path):
    """Plot bond formation and breaking dynamics"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Bond Dynamics Analysis', fontsize=16)
    
    # Total bonds over time
    ax = axes[0, 0]
    ax.plot(df['sim_time'], df['num_bonds_total'], 'b-', linewidth=1)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Total Bonds')
    ax.set_title('Bond Count Evolution')
    ax.grid(True, alpha=0.3)
    
    # Bond length and tension
    ax = axes[0, 1]
    ax.plot(df['sim_time'], df['avg_bond_length'], 'g-', label='Length', linewidth=1)
    ax_twin = ax.twinx()
    ax_twin.plot(df['sim_time'], df['avg_bond_tension'], 'r-', label='Tension', linewidth=1)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Average Bond Length', color='g')
    ax_twin.set_ylabel('Average Bond Tension', color='r')
    ax.set_title('Bond Length and Tension')
    ax.grid(True, alpha=0.3)
    
    # Formation vs Breaking events
    ax = axes[1, 0]
    window = min(20, len(df) // 10)
    if window > 0:
        formed_smooth = savgol_filter(df['events_formed'], window, 3) if len(df) > window else df['events_formed']
        broken_smooth = savgol_filter(df['events_broken'], window, 3) if len(df) > window else df['events_broken']
        ax.plot(df['sim_time'], formed_smooth, 'g-', label='Formed', linewidth=1.5)
        ax.plot(df['sim_time'], broken_smooth, 'r-', label='Broken', linewidth=1.5)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Events per Step')
    ax.set_title('Bond Formation/Breaking Events')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Cumulative events
    ax = axes[1, 1]
    ax.plot(df['sim_time'], df['events_formed'].cumsum(), 'g-', label='Formed', linewidth=1.5)
    ax.plot(df['sim_time'], df['events_broken'].cumsum(), 'r-', label='Broken', linewidth=1.5)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Cumulative Events')
    ax.set_title('Cumulative Bond Events')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = output_dir / 'bond_dynamics.png'
    plt.savefig(output_file, dpi=150)
    print(f"Saved bond dynamics plot: {output_file}")
    plt.close()

def plot_cluster_dynamics(df: pd.DataFrame, output_dir: Path):
    """Plot cluster formation and evolution"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Cluster Dynamics Analysis', fontsize=16)
    
    # Number of clusters
    ax = axes[0, 0]
    ax.plot(df['sim_time'], df['num_clusters'], 'b-', linewidth=1)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Number of Clusters')
    ax.set_title('Cluster Count Evolution')
    ax.grid(True, alpha=0.3)
    
    # Largest cluster size
    ax = axes[0, 1]
    ax.plot(df['sim_time'], df['largest_cluster_size'], 'r-', linewidth=1)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Largest Cluster Size')
    ax.set_title('Largest Cluster Evolution')
    ax.grid(True, alpha=0.3)
    
    # Cluster energy statistics
    ax = axes[1, 0]
    ax.plot(df['sim_time'], df['cluster_energy_mean'], 'g-', label='Mean', linewidth=1)
    ax.fill_between(df['sim_time'], 
                     df['cluster_energy_mean'] - np.sqrt(df['cluster_energy_var']),
                     df['cluster_energy_mean'] + np.sqrt(df['cluster_energy_var']),
                     alpha=0.3, color='g')
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Cluster Energy')
    ax.set_title('Cluster Energy Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Radius of gyration
    ax = axes[1, 1]
    ax.plot(df['sim_time'], df['R_g_mean'], 'purple', linewidth=1)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Mean Radius of Gyration')
    ax.set_title('Cluster Spatial Extent')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = output_dir / 'cluster_dynamics.png'
    plt.savefig(output_file, dpi=150)
    print(f"Saved cluster dynamics plot: {output_file}")
    plt.close()

def plot_phase_analysis(df: pd.DataFrame, output_dir: Path):
    """Detect and plot phase transitions"""
    fig, axes = plt.subplots(2, 1, figsize=(15, 10))
    fig.suptitle('Phase Transition Analysis', fontsize=16)
    
    # Compute order parameters
    df['bond_density'] = df['num_bonds_total'] / df['num_particles']
    df['cluster_fraction'] = df['largest_cluster_size'] / df['num_particles']
    
    # Order parameters
    ax = axes[0]
    ax.plot(df['sim_time'], df['bond_density'], 'b-', label='Bond Density', linewidth=1.5)
    ax.plot(df['sim_time'], df['cluster_fraction'], 'r-', label='Largest Cluster Fraction', linewidth=1.5)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Order Parameter')
    ax.set_title('Order Parameters (Phase Indicators)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Event rate (activity)
    ax = axes[1]
    df['total_events'] = df['events_formed'] + df['events_broken'] + df['events_merged'] + df['events_split']
    window = min(20, len(df) // 10)
    if window > 0:
        activity = savgol_filter(df['total_events'], window, 3) if len(df) > window else df['total_events']
        ax.plot(df['sim_time'], activity, 'g-', linewidth=1.5)
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Total Events per Step')
    ax.set_title('System Activity Level')
    ax.grid(True, alpha=0.3)
    
    # Detect peaks (high activity phases)
    if len(df) > 10:
        peaks, properties = find_peaks(df['total_events'], prominence=df['total_events'].std())
        if len(peaks) > 0:
            ax.scatter(df['sim_time'].iloc[peaks], df['total_events'].iloc[peaks], 
                      color='red', s=100, marker='*', label='High Activity Peaks', zorder=5)
            ax.legend()
    
    plt.tight_layout()
    output_file = output_dir / 'phase_analysis.png'
    plt.savefig(output_file, dpi=150)
    print(f"Saved phase analysis plot: {output_file}")
    plt.close()

def print_summary_statistics(df: pd.DataFrame):
    """Print summary statistics"""
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    print(f"\nSimulation Duration:")
    print(f"  Steps: {df['step'].max()}")
    print(f"  Time: {df['sim_time'].max():.2f}")
    print(f"  Wall Time: {df['wall_time'].max():.2f} seconds")
    
    print(f"\nBond Statistics:")
    print(f"  Max Bonds: {df['num_bonds_total'].max()}")
    print(f"  Mean Bonds: {df['num_bonds_total'].mean():.1f}")
    print(f"  Total Formed: {df['events_formed'].sum()}")
    print(f"  Total Broken: {df['events_broken'].sum()}")
    print(f"  Net Change: {df['events_formed'].sum() - df['events_broken'].sum()}")
    
    print(f"\nCluster Statistics:")
    print(f"  Max Clusters: {df['num_clusters'].max()}")
    print(f"  Mean Clusters: {df['num_clusters'].mean():.1f}")
    print(f"  Max Cluster Size: {df['largest_cluster_size'].max()}")
    print(f"  Total Merges: {df['events_merged'].sum()}")
    print(f"  Total Splits: {df['events_split'].sum()}")
    
    print(f"\nSystem Activity:")
    total_events = df['events_formed'] + df['events_broken'] + df['events_merged'] + df['events_split']
    print(f"  Total Events: {total_events.sum()}")
    print(f"  Events/Step (mean): {total_events.mean():.2f}")
    print(f"  Events/Step (max): {total_events.max()}")
    
    print("\n" + "="*60 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Analyze Live 2.0 diagnostics data')
    parser.add_argument('path', type=str, help='Path to diagnostics directory or metrics CSV file')
    parser.add_argument('-o', '--output', type=str, default=None, 
                       help='Output directory for plots (default: same as input)')
    args = parser.parse_args()
    
    # Load data
    path = Path(args.path)
    if not path.exists():
        print(f"Error: Path not found: {path}")
        sys.exit(1)
    
    data = load_diagnostics(path)
    df = data['metrics']
    
    # Output directory
    if args.output:
        output_dir = Path(args.output)
    elif path.is_file():
        output_dir = path.parent
    else:
        output_dir = path
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Print summary
    print_summary_statistics(df)
    
    # Generate plots
    print("Generating plots...")
    plot_bond_dynamics(df, output_dir)
    plot_cluster_dynamics(df, output_dir)
    plot_phase_analysis(df, output_dir)
    
    print(f"\nAnalysis complete! Plots saved to: {output_dir}")

if __name__ == '__main__':
    main()

