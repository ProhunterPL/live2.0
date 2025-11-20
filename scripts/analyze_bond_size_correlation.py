#!/usr/bin/env python3
"""
Analyze bond-size correlation to check if results are realistic
"""
import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

def analyze_bond_size():
    # Load data
    with open("analysis/phase2b_miller_urey/batch_analysis.json", 'r') as f:
        data = json.load(f)
    
    # Collect size-bond pairs
    sizes = []
    bonds = []
    counts = []
    
    for run in data['runs']:
        for mol in run.get('molecules', []):
            size = mol.get('size', 0)
            n_bonds = len(mol.get('bonds', []))
            count = mol.get('count', 1)
            
            sizes.append(size)
            bonds.append(n_bonds)
            counts.append(count)
    
    sizes = np.array(sizes)
    bonds = np.array(bonds)
    counts = np.array(counts)
    
    # Print statistics
    print("=" * 80)
    print("BOND-SIZE CORRELATION ANALYSIS")
    print("=" * 80)
    print()
    
    print("Overall Statistics:")
    print(f"  Total molecules: {len(sizes)}")
    print(f"  Size range: {sizes.min()} - {sizes.max()} atoms")
    print(f"  Bond range: {bonds.min()} - {bonds.max()} bonds")
    print()
    
    # Expected bonds for linear molecule: size - 1
    # Expected bonds for fully connected: size * (size-1) / 2
    print("Theoretical Expectations:")
    print("  Linear molecule (n atoms):     n-1 bonds")
    print("  Branched/cyclic (n atoms):     n to 3n bonds")
    print("  Cluster/aggregate (n atoms):   << n bonds (weak connections)")
    print()
    
    # Analyze problematic cases
    problematic = []
    for i, (s, b) in enumerate(zip(sizes, bonds)):
        if s > 10 and b < s * 0.3:  # Large molecule with very few bonds
            problematic.append((s, b, counts[i]))
    
    print(f"Potentially Problematic Cases (size > 10, bonds < 30% of size):")
    print(f"  Found: {len(problematic)} cases")
    print()
    
    if problematic:
        print("Examples (size, bonds, count):")
        for s, b, c in sorted(problematic, key=lambda x: x[0], reverse=True)[:20]:
            ratio = b / s if s > 0 else 0
            print(f"    {s:3d} atoms, {b:3d} bonds ({ratio:.2f} ratio), count={c}")
        print()
    
    # Categorize molecules
    print("Molecular Categories:")
    
    monomers = sum((s == 1) for s in sizes)
    dimers = sum((s == 2) for s in sizes)
    small_linear = sum((s >= 3) and (s <= 5) and (b == s - 1) for s, b in zip(sizes, bonds))
    small_branched = sum((s >= 3) and (s <= 5) and (b >= s) for s, b in zip(sizes, bonds))
    medium = sum((s >= 6) and (s <= 15) for s in sizes)
    large = sum(s > 15 for s in sizes)
    
    print(f"  Monomers (1 atom):              {monomers:6d} ({100*monomers/len(sizes):.1f}%)")
    print(f"  Dimers (2 atoms):               {dimers:6d} ({100*dimers/len(sizes):.1f}%)")
    print(f"  Small linear (3-5 atoms, linear): {small_linear:6d} ({100*small_linear/len(sizes):.1f}%)")
    print(f"  Small branched (3-5 atoms, >n-1): {small_branched:6d} ({100*small_branched/len(sizes):.1f}%)")
    print(f"  Medium (6-15 atoms):            {medium:6d} ({100*medium/len(sizes):.1f}%)")
    print(f"  Large (>15 atoms):              {large:6d} ({100*large/len(sizes):.1f}%)")
    print()
    
    # Check bond/size ratio
    print("Bond/Size Ratio Distribution:")
    for threshold in [0.1, 0.3, 0.5, 0.7, 0.9]:
        mask = (sizes > 5)  # Only molecules larger than 5 atoms
        if mask.sum() > 0:
            ratios = bonds[mask] / sizes[mask]
            below = sum(ratios < threshold)
            print(f"  Molecules >5 atoms with ratio < {threshold}: {below}/{mask.sum()} ({100*below/mask.sum():.1f}%)")
    print()
    
    # Create visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Scatter plot: bonds vs size
    ax1.scatter(sizes, bonds, alpha=0.3, s=20, c='steelblue')
    ax1.plot([0, 30], [0, 30], 'r--', label='Linear (bonds = atoms)', alpha=0.5)
    ax1.plot([0, 30], [0, 60], 'orange', linestyle='--', label='Branched (bonds = 2*atoms)', alpha=0.5)
    ax1.set_xlabel('Number of Atoms', fontsize=11)
    ax1.set_ylabel('Number of Bonds', fontsize=11)
    ax1.set_title('Bonds vs Size Correlation', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 35)
    ax1.set_ylim(0, 70)
    
    # 2. Histogram of bond/size ratio (for molecules > 5 atoms)
    mask = sizes > 5
    if mask.sum() > 0:
        ratios = bonds[mask] / sizes[mask]
        ax2.hist(ratios, bins=30, color='coral', alpha=0.7, edgecolor='black')
        ax2.axvline(1.0, color='red', linestyle='--', label='Linear (ratio=1)')
        ax2.set_xlabel('Bonds/Size Ratio', fontsize=11)
        ax2.set_ylabel('Frequency', fontsize=11)
        ax2.set_title('Bond/Size Ratio (molecules > 5 atoms)', fontsize=12, fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
    
    # 3. Box plot: bonds by size category
    size_categories = []
    bond_lists = {f'{i}-{i+4}': [] for i in range(2, 30, 5)}
    bond_lists['30+'] = []
    
    for s, b in zip(sizes, bonds):
        if s >= 30:
            bond_lists['30+'].append(b)
        else:
            for i in range(2, 30, 5):
                if i <= s < i+5:
                    bond_lists[f'{i}-{i+4}'].append(b)
                    break
    
    categories = [k for k, v in bond_lists.items() if len(v) > 0]
    bond_data = [bond_lists[k] for k in categories]
    
    ax3.boxplot(bond_data, labels=categories)
    ax3.set_xlabel('Size Category (atoms)', fontsize=11)
    ax3.set_ylabel('Number of Bonds', fontsize=11)
    ax3.set_title('Bond Distribution by Size Category', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.tick_params(axis='x', rotation=45)
    
    # 4. Hexbin plot for density
    ax4.hexbin(sizes, bonds, gridsize=20, cmap='YlOrRd', mincnt=1)
    ax4.plot([0, 30], [0, 30], 'b--', label='Linear', alpha=0.7, linewidth=2)
    ax4.set_xlabel('Number of Atoms', fontsize=11)
    ax4.set_ylabel('Number of Bonds', fontsize=11)
    ax4.set_title('Density Plot: Bonds vs Size', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.set_xlim(0, 35)
    ax4.set_ylim(0, 70)
    
    plt.suptitle('Miller-Urey: Bond-Size Correlation Analysis', fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig('analysis/phase2b_miller_urey/figures/figure5_bond_size_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Generated: figure5_bond_size_correlation.png")
    print()
    print("=" * 80)
    print("INTERPRETATION")
    print("=" * 80)
    print()
    print("If you see many large molecules with very few bonds:")
    print("  -> May indicate CLUSTERS (weakly bound aggregates) rather than molecules")
    print("  -> This is EXPECTED in MD simulations (transient associations)")
    print("  -> Not necessarily chemical bonds, but spatial proximity")
    print()
    print("Realistic molecules should have:")
    print("  - Linear:   bonds ~ atoms - 1")
    print("  - Branched: bonds ~ atoms to 1.5*atoms")
    print("  - Cyclic:   bonds ~ atoms")
    print()
    print("Large clusters with few bonds are PHYSICAL, not CHEMICAL species.")
    print("=" * 80)
    print()
    print("RECOMMENDATION:")
    print("  Filter molecules to keep only those with reasonable bond/size ratios")
    print("  Example: Keep molecules where bonds >= 0.5 * (size - 1)")
    print("  This will give TRUE MOLECULES, not spatial clusters")
    print("=" * 80)

if __name__ == "__main__":
    analyze_bond_size()

