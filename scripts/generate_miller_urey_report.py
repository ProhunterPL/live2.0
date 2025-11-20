#!/usr/bin/env python3
"""
Generate detailed Miller-Urey analysis report
"""
import json
import sys
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime

def load_batch_analysis(analysis_file):
    """Load batch analysis JSON"""
    with open(analysis_file, 'r') as f:
        return json.load(f)

def analyze_top_molecules(data):
    """Analyze top molecules across all runs"""
    
    # Collect all molecules
    all_molecules = []
    molecule_formulas = Counter()
    molecule_runs = defaultdict(set)
    
    for run in data['runs']:
        run_name = run['run_name']
        for mol in run.get('molecules', []):
            formula = mol.get('formula', 'Unknown')
            count = mol.get('count', 1)
            
            all_molecules.append({
                'formula': formula,
                'count': count,
                'run': run_name,
                'size': mol.get('size', 0),
                'bonds': len(mol.get('bonds', []))
            })
            
            molecule_formulas[formula] += count
            molecule_runs[formula].add(run_name)
    
    return all_molecules, molecule_formulas, molecule_runs

def generate_report(analysis_file, output_file):
    """Generate detailed report"""
    
    print("Loading analysis data...")
    data = load_batch_analysis(analysis_file)
    
    print("Analyzing molecules...")
    all_molecules, molecule_formulas, molecule_runs = analyze_top_molecules(data)
    
    # Generate report
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("MILLER-UREY PHASE 2B - DETAILED ANALYSIS REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"Analysis file: {analysis_file}")
    report_lines.append("")
    
    # Overall summary
    report_lines.append("=" * 80)
    report_lines.append("OVERALL SUMMARY")
    report_lines.append("=" * 80)
    report_lines.append(f"Total runs: {data['summary']['total_runs']}")
    report_lines.append(f"Successful runs: {data['summary']['successful_runs']}")
    report_lines.append(f"Unique molecules (all runs): {len(molecule_formulas)}")
    report_lines.append(f"Total molecule instances: {sum(molecule_formulas.values())}")
    report_lines.append("")
    
    # Per-run statistics
    report_lines.append("=" * 80)
    report_lines.append("PER-RUN STATISTICS")
    report_lines.append("=" * 80)
    report_lines.append(f"{'Run':<10} {'Molecules':<12} {'Instances':<12} {'Avg Size':<12}")
    report_lines.append("-" * 80)
    
    for run in sorted(data['runs'], key=lambda x: x['run_name']):
        run_name = run['run_name']
        molecules = run.get('molecules', [])
        n_mol = len(molecules)
        n_inst = sum(m.get('count', 1) for m in molecules)
        avg_size = sum(m.get('size', 0) * m.get('count', 1) for m in molecules) / n_inst if n_inst > 0 else 0
        
        report_lines.append(f"{run_name:<10} {n_mol:<12} {n_inst:<12} {avg_size:<12.2f}")
    
    report_lines.append("")
    
    # Top 50 molecules by occurrence
    report_lines.append("=" * 80)
    report_lines.append("TOP 50 MOLECULES (by total occurrences across all runs)")
    report_lines.append("=" * 80)
    report_lines.append(f"{'Rank':<6} {'Formula Hash':<35} {'Count':<10} {'Runs':<8}")
    report_lines.append("-" * 80)
    
    for i, (formula, count) in enumerate(molecule_formulas.most_common(50), 1):
        n_runs = len(molecule_runs[formula])
        report_lines.append(f"{i:<6} {formula[:33]:<35} {count:<10} {n_runs:<8}")
    
    report_lines.append("")
    
    # Molecules by run frequency
    report_lines.append("=" * 80)
    report_lines.append("MOST FREQUENT MOLECULES (appearing in most runs)")
    report_lines.append("=" * 80)
    report_lines.append(f"{'Rank':<6} {'Formula Hash':<35} {'Runs':<8} {'Total':<10}")
    report_lines.append("-" * 80)
    
    by_frequency = sorted(molecule_runs.items(), key=lambda x: (len(x[1]), molecule_formulas[x[0]]), reverse=True)
    for i, (formula, runs) in enumerate(by_frequency[:50], 1):
        n_runs = len(runs)
        total_count = molecule_formulas[formula]
        report_lines.append(f"{i:<6} {formula[:33]:<35} {n_runs:<8} {total_count:<10}")
    
    report_lines.append("")
    
    # Complexity distribution
    report_lines.append("=" * 80)
    report_lines.append("MOLECULE COMPLEXITY DISTRIBUTION")
    report_lines.append("=" * 80)
    
    size_dist = Counter()
    bond_dist = Counter()
    
    for mol in all_molecules:
        size = mol['size']
        bonds = mol['bonds']
        count = mol['count']
        
        size_dist[size] += count
        bond_dist[bonds] += count
    
    report_lines.append("\nBy Size (atoms):")
    report_lines.append(f"{'Atoms':<10} {'Count':<15} {'Percent':<10}")
    report_lines.append("-" * 40)
    total = sum(size_dist.values())
    for size in sorted(size_dist.keys())[:20]:
        count = size_dist[size]
        pct = 100 * count / total
        report_lines.append(f"{size:<10} {count:<15} {pct:<10.2f}%")
    
    report_lines.append("\nBy Bonds:")
    report_lines.append(f"{'Bonds':<10} {'Count':<15} {'Percent':<10}")
    report_lines.append("-" * 40)
    for bonds in sorted(bond_dist.keys())[:20]:
        count = bond_dist[bonds]
        pct = 100 * count / total
        report_lines.append(f"{bonds:<10} {count:<15} {pct:<10.2f}%")
    
    report_lines.append("")
    report_lines.append("=" * 80)
    report_lines.append("END OF REPORT")
    report_lines.append("=" * 80)
    
    # Write report
    report_text = "\n".join(report_lines)
    
    with open(output_file, 'w') as f:
        f.write(report_text)
    
    print(f"\nReport generated: {output_file}")
    print(report_text)
    
    return report_text

if __name__ == "__main__":
    analysis_file = "analysis/phase2b_miller_urey/batch_analysis.json"
    output_file = "analysis/phase2b_miller_urey/detailed_report.txt"
    
    if len(sys.argv) > 1:
        analysis_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    generate_report(analysis_file, output_file)

