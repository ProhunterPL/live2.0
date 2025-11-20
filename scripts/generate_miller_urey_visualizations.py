#!/usr/bin/env python3
"""
Generate visualizations for Miller-Urey Phase 2B results
"""
import json
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter

# Set style
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'

def load_data(analysis_file):
    """Load batch analysis"""
    with open(analysis_file, 'r') as f:
        return json.load(f)

def plot_molecules_per_run(data, output_dir):
    """Plot number of molecules per run"""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    runs = []
    molecules = []
    instances = []
    
    for run in sorted(data['runs'], key=lambda x: int(x['run_name'].split('_')[1])):
        runs.append(run['run_name'])
        mols = run.get('molecules', [])
        molecules.append(len(mols))
        instances.append(sum(m.get('count', 1) for m in mols))
    
    x = np.arange(len(runs))
    width = 0.35
    
    ax.bar(x - width/2, molecules, width, label='Unique Molecules', color='steelblue')
    ax.bar(x + width/2, instances, width, label='Total Instances', color='coral', alpha=0.7)
    
    ax.set_xlabel('Run Number', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Miller-Urey Phase 2B: Molecules Per Run', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([r.split('_')[1] for r in runs], rotation=45)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'figure1_molecules_per_run.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated: figure1_molecules_per_run.png")

def plot_top_molecules(data, output_dir):
    """Plot top 20 molecules"""
    # Count molecules
    molecule_counts = Counter()
    for run in data['runs']:
        for mol in run.get('molecules', []):
            formula = mol.get('formula', 'Unknown')
            count = mol.get('count', 1)
            molecule_counts[formula] += count
    
    # Get top 20
    top_20 = molecule_counts.most_common(20)
    formulas = [f[:12] + '...' for f, _ in top_20]  # Shorten hash
    counts = [c for _, c in top_20]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    y_pos = np.arange(len(formulas))
    ax.barh(y_pos, counts, color='mediumseagreen')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(formulas, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel('Total Occurrences', fontsize=12)
    ax.set_title('Top 20 Molecules by Occurrence', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add values
    for i, v in enumerate(counts):
        ax.text(v + 20, i, str(v), va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'figure2_top_molecules.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated: figure2_top_molecules.png")

def plot_complexity_distribution(data, output_dir):
    """Plot molecule complexity distribution"""
    
    size_dist = Counter()
    bond_dist = Counter()
    
    for run in data['runs']:
        for mol in run.get('molecules', []):
            size = mol.get('size', 0)
            bonds = len(mol.get('bonds', []))
            count = mol.get('count', 1)
            
            size_dist[size] += count
            bond_dist[bonds] += count
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Size distribution
    sizes = sorted([k for k in size_dist.keys() if k > 0 and k <= 30])
    size_counts = [size_dist[s] for s in sizes]
    
    ax1.bar(sizes, size_counts, color='royalblue', alpha=0.7)
    ax1.set_xlabel('Number of Atoms', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Molecule Size Distribution', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Bond distribution
    bonds = sorted([k for k in bond_dist.keys() if k >= 0 and k <= 20])
    bond_counts = [bond_dist[b] for b in bonds]
    
    ax2.bar(bonds, bond_counts, color='darkorange', alpha=0.7)
    ax2.set_xlabel('Number of Bonds', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Bond Complexity Distribution', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('Miller-Urey: Molecular Complexity', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'figure3_complexity_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated: figure3_complexity_distribution.png")

def plot_summary_statistics(data, output_dir):
    """Plot summary statistics"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Unique vs Total molecules
    runs_data = sorted(data['runs'], key=lambda x: int(x['run_name'].split('_')[1]))
    run_nums = [int(r['run_name'].split('_')[1]) for r in runs_data]
    unique_mols = [len(r.get('molecules', [])) for r in runs_data]
    total_insts = [sum(m.get('count', 1) for m in r.get('molecules', [])) for r in runs_data]
    
    ax1.scatter(run_nums, unique_mols, s=100, alpha=0.6, c='steelblue', label='Unique')
    ax1.axhline(np.mean(unique_mols), color='red', linestyle='--', alpha=0.5, label=f'Mean: {np.mean(unique_mols):.1f}')
    ax1.set_xlabel('Run Number', fontsize=11)
    ax1.set_ylabel('Unique Molecules', fontsize=11)
    ax1.set_title('Unique Molecules per Run', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Total instances
    ax2.scatter(run_nums, total_insts, s=100, alpha=0.6, c='coral')
    ax2.axhline(np.mean(total_insts), color='red', linestyle='--', alpha=0.5, label=f'Mean: {np.mean(total_insts):.1f}')
    ax2.set_xlabel('Run Number', fontsize=11)
    ax2.set_ylabel('Total Instances', fontsize=11)
    ax2.set_title('Total Molecule Instances per Run', fontsize=12, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Average molecule size per run
    avg_sizes = []
    for run in runs_data:
        mols = run.get('molecules', [])
        if mols:
            total_size = sum(m.get('size', 0) * m.get('count', 1) for m in mols)
            total_count = sum(m.get('count', 1) for m in mols)
            avg_sizes.append(total_size / total_count if total_count > 0 else 0)
        else:
            avg_sizes.append(0)
    
    ax3.bar(run_nums, avg_sizes, color='mediumseagreen', alpha=0.7)
    ax3.axhline(np.mean(avg_sizes), color='red', linestyle='--', alpha=0.5, label=f'Mean: {np.mean(avg_sizes):.2f}')
    ax3.set_xlabel('Run Number', fontsize=11)
    ax3.set_ylabel('Average Size (atoms)', fontsize=11)
    ax3.set_title('Average Molecule Size per Run', fontsize=12, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Distribution histogram
    ax4.hist(unique_mols, bins=15, color='purple', alpha=0.6, edgecolor='black')
    ax4.axvline(np.mean(unique_mols), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(unique_mols):.1f}')
    ax4.axvline(np.median(unique_mols), color='orange', linestyle='--', linewidth=2, label=f'Median: {np.median(unique_mols):.1f}')
    ax4.set_xlabel('Unique Molecules', fontsize=11)
    ax4.set_ylabel('Frequency', fontsize=11)
    ax4.set_title('Distribution of Unique Molecules', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('Miller-Urey Phase 2B: Summary Statistics', fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(output_dir / 'figure4_summary_statistics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated: figure4_summary_statistics.png")

def main():
    analysis_file = Path("analysis/phase2b_miller_urey/batch_analysis.json")
    output_dir = Path("analysis/phase2b_miller_urey/figures")
    
    if len(sys.argv) > 1:
        analysis_file = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading data...")
    data = load_data(analysis_file)
    
    print("\nGenerating visualizations...")
    plot_molecules_per_run(data, output_dir)
    plot_top_molecules(data, output_dir)
    plot_complexity_distribution(data, output_dir)
    plot_summary_statistics(data, output_dir)
    
    print(f"\nAll visualizations saved to: {output_dir}/")
    print("\nGenerated figures:")
    print("  1. figure1_molecules_per_run.png")
    print("  2. figure2_top_molecules.png")
    print("  3. figure3_complexity_distribution.png")
    print("  4. figure4_summary_statistics.png")

if __name__ == "__main__":
    main()

