"""
Scenario Comparison Tool
=========================

Compare results across different Phase 2 scenarios.

Usage:
    # Compare all scenarios
    python scripts/compare_scenarios.py results/phase2
    
    # Compare specific scenarios
    python scripts/compare_scenarios.py \
        results/phase2/miller_urey \
        results/phase2/hydrothermal \
        results/phase2/formamide
"""

import sys
import argparse
import logging
import json
from pathlib import Path
from collections import defaultdict
from typing import List, Dict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ScenarioComparator:
    """Compare results across different scenarios"""
    
    def __init__(self, scenario_dirs: List[Path], output_dir: Path):
        self.scenario_dirs = scenario_dirs
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.results = {}
        
    def load_results(self):
        """Load results from all scenario directories"""
        logger.info("Loading results from all scenarios...")
        
        for scenario_dir in self.scenario_dirs:
            scenario_name = scenario_dir.name
            logger.info(f"  Loading: {scenario_name}")
            
            # Find all result directories
            result_dirs = [d for d in scenario_dir.rglob('*') 
                          if d.is_dir() and (d / 'results.json').exists()]
            
            scenario_results = {
                'name': scenario_name,
                'num_runs': len(result_dirs),
                'molecules': [],
                'metrics': []
            }
            
            for result_dir in result_dirs:
                try:
                    # Load results.json
                    with open(result_dir / 'results.json') as f:
                        data = json.load(f)
                        
                    # Extract molecules
                    if 'molecules' in data:
                        scenario_results['molecules'].extend(data['molecules'])
                    
                    # Extract metrics
                    if 'metrics' in data:
                        scenario_results['metrics'].append(data['metrics'])
                        
                except Exception as e:
                    logger.debug(f"Failed to load {result_dir}: {e}")
            
            self.results[scenario_name] = scenario_results
            logger.info(f"    {len(scenario_results['molecules'])} molecules, "
                       f"{len(scenario_results['metrics'])} runs")
    
    def compare_molecule_diversity(self) -> Dict:
        """Compare molecular diversity across scenarios"""
        logger.info("\nComparing molecular diversity...")
        
        comparison = {}
        
        for name, data in self.results.items():
            # Unique molecules by formula
            formulas = set()
            for mol in data['molecules']:
                if 'formula' in mol:
                    formulas.add(mol['formula'])
            
            # Molecule size distribution
            sizes = []
            for mol in data['molecules']:
                if 'atoms' in mol:
                    sizes.append(len(mol['atoms']))
            
            comparison[name] = {
                'unique_molecules': len(formulas),
                'total_molecules': len(data['molecules']),
                'avg_molecule_size': sum(sizes) / len(sizes) if sizes else 0,
                'max_molecule_size': max(sizes) if sizes else 0,
                'min_molecule_size': min(sizes) if sizes else 0
            }
        
        return comparison
    
    def compare_reaction_rates(self) -> Dict:
        """Compare reaction rates across scenarios"""
        logger.info("Comparing reaction rates...")
        
        comparison = {}
        
        for name, data in self.results.items():
            if not data['metrics']:
                continue
            
            # Average metrics across runs
            total_reactions = 0
            total_steps = 0
            
            for metrics in data['metrics']:
                if isinstance(metrics, dict):
                    total_reactions += metrics.get('total_reactions', 0)
                    total_steps += metrics.get('total_steps', 0)
            
            comparison[name] = {
                'total_reactions': total_reactions,
                'avg_reactions_per_run': total_reactions / len(data['metrics']) if data['metrics'] else 0,
                'reactions_per_step': total_reactions / total_steps if total_steps > 0 else 0
            }
        
        return comparison
    
    def generate_report(self):
        """Generate comparison report"""
        logger.info("\nGenerating comparison report...")
        
        # Compare metrics
        diversity = self.compare_molecule_diversity()
        reactions = self.compare_reaction_rates()
        
        # Write report
        report_file = self.output_dir / "scenario_comparison.txt"
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("PHASE 2 SCENARIO COMPARISON\n")
            f.write("=" * 70 + "\n\n")
            
            # Overview
            f.write("Overview:\n")
            f.write("-" * 70 + "\n")
            for name, data in self.results.items():
                f.write(f"  {name}:\n")
                f.write(f"    Runs: {data['num_runs']}\n")
                f.write(f"    Total molecules: {len(data['molecules'])}\n")
                f.write("\n")
            
            # Molecular diversity
            f.write("\nMolecular Diversity:\n")
            f.write("-" * 70 + "\n")
            for name, metrics in diversity.items():
                f.write(f"  {name}:\n")
                f.write(f"    Unique molecules: {metrics['unique_molecules']}\n")
                f.write(f"    Avg size: {metrics['avg_molecule_size']:.1f} atoms\n")
                f.write(f"    Size range: {metrics['min_molecule_size']}-{metrics['max_molecule_size']} atoms\n")
                f.write("\n")
            
            # Reaction rates
            f.write("\nReaction Rates:\n")
            f.write("-" * 70 + "\n")
            for name, metrics in reactions.items():
                f.write(f"  {name}:\n")
                f.write(f"    Total reactions: {metrics['total_reactions']}\n")
                f.write(f"    Avg per run: {metrics['avg_reactions_per_run']:.1f}\n")
                f.write(f"    Per step: {metrics['reactions_per_step']:.6f}\n")
                f.write("\n")
            
            # Summary
            f.write("\nSummary:\n")
            f.write("-" * 70 + "\n")
            
            # Most productive scenario
            if diversity:
                most_diverse = max(diversity.items(), key=lambda x: x[1]['unique_molecules'])
                f.write(f"  Most diverse: {most_diverse[0]} ({most_diverse[1]['unique_molecules']} unique molecules)\n")
            
            if reactions:
                most_reactive = max(reactions.items(), key=lambda x: x[1]['total_reactions'])
                f.write(f"  Most reactive: {most_reactive[0]} ({most_reactive[1]['total_reactions']} total reactions)\n")
            
            f.write("\n" + "=" * 70 + "\n")
        
        logger.info(f"Report saved to: {report_file}")
        
        # Save JSON for further analysis
        json_file = self.output_dir / "comparison_data.json"
        with open(json_file, 'w') as f:
            json.dump({
                'diversity': diversity,
                'reactions': reactions
            }, f, indent=2)
        
        logger.info(f"JSON data saved to: {json_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare Phase 2 scenarios"
    )
    parser.add_argument(
        'directories',
        nargs='+',
        type=str,
        help="Scenario directories to compare"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='analysis/scenario_comparison',
        help="Output directory (default: analysis/scenario_comparison)"
    )
    
    args = parser.parse_args()
    
    # Convert to Path objects
    scenario_dirs = [Path(d) for d in args.directories]
    output_dir = Path(args.output)
    
    # Validate directories
    for d in scenario_dirs:
        if not d.exists():
            logger.error(f"Directory not found: {d}")
            return 1
    
    logger.info("=" * 70)
    logger.info("SCENARIO COMPARATOR")
    logger.info("=" * 70)
    logger.info(f"Scenarios: {len(scenario_dirs)}")
    for d in scenario_dirs:
        logger.info(f"  - {d}")
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 70)
    
    # Run comparison
    comparator = ScenarioComparator(scenario_dirs, output_dir)
    comparator.load_results()
    comparator.generate_report()
    
    logger.info("\n" + "=" * 70)
    logger.info("COMPARISON COMPLETE!")
    logger.info("=" * 70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

