"""
Phase 2 Batch Results Analyzer
================================

Analyzes results from multiple Phase 2 simulation runs.
Extracts molecules, matches with PubChem, generates reports.

Usage:
    python scripts/analyze_phase2_batch.py \
        --input results/phase2/miller_urey \
        --output analysis/miller_urey \
        --use-matcher

For all scenarios:
    python scripts/analyze_phase2_batch.py \
        --input results/phase2 \
        --output analysis/phase2_complete \
        --recursive \
        --use-matcher
"""

import sys
import argparse
import logging
import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict
from collections import defaultdict

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.sim.molecule_extractor import MoleculeExtractor, extract_molecules_from_results

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Phase2BatchAnalyzer:
    """Analyzes batch Phase 2 simulation results"""
    
    def __init__(self, input_dir: str, output_dir: str, use_matcher: bool = True):
        """
        Initialize batch analyzer
        
        Args:
            input_dir: Directory containing simulation results
            output_dir: Output directory for analysis
            use_matcher: Use MatcherV2 for PubChem matching
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.use_matcher = use_matcher
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("=" * 70)
        logger.info("PHASE 2 BATCH ANALYZER")
        logger.info("=" * 70)
        logger.info(f"Input: {self.input_dir}")
        logger.info(f"Output: {self.output_dir}")
        logger.info(f"Use MatcherV2: {self.use_matcher}")
        logger.info("=" * 70)
    
    def find_result_directories(self, recursive: bool = False) -> List[Path]:
        """
        Find all result directories containing results.json
        
        Args:
            recursive: Search recursively
        
        Returns:
            List of result directory paths
        """
        logger.info("Finding result directories...")
        
        if recursive:
            # Find all directories containing results.json
            result_dirs = []
            for results_file in self.input_dir.rglob("results.json"):
                result_dirs.append(results_file.parent)
        else:
            # Only direct subdirectories
            result_dirs = []
            for subdir in self.input_dir.iterdir():
                if subdir.is_dir() and (subdir / "results.json").exists():
                    result_dirs.append(subdir)
        
        logger.info(f"  Found {len(result_dirs)} result directories")
        for d in result_dirs:
            logger.info(f"    - {d.relative_to(self.input_dir)}")
        
        return sorted(result_dirs)
    
    def analyze_single_run(self, results_dir: Path) -> Dict:
        """
        Analyze single simulation run
        
        Args:
            results_dir: Directory containing simulation results
        
        Returns:
            Analysis results dictionary
        """
        logger.info(f"\nAnalyzing: {results_dir.name}")
        logger.info("-" * 70)
        
        try:
            # Extract molecules
            extraction = extract_molecules_from_results(
                str(results_dir),
                output_dir=str(results_dir / "analysis"),
                export_for_matcher=self.use_matcher
            )
            
            # Load results.json for metadata
            with open(results_dir / "results.json", 'r') as f:
                results = json.load(f)
            
            analysis = {
                'directory': str(results_dir),
                'run_name': results_dir.name,
                'scenario': results.get('scenario', 'unknown'),
                'config': results.get('configuration', {}),
                'molecules': extraction['molecules'],
                'statistics': extraction['statistics'],
                'analysis_dir': extraction['output_dir'],
                'report_file': extraction['report_file'],
                'timestamp': datetime.now().isoformat()
            }
            
            logger.info(f"  Molecules found: {len(extraction['molecules'])}")
            logger.info(f"  Analysis saved to: {extraction['output_dir']}")
            
            return analysis
        
        except Exception as e:
            logger.error(f"  Failed to analyze {results_dir.name}: {e}")
            return {
                'directory': str(results_dir),
                'run_name': results_dir.name,
                'error': str(e)
            }
    
    def analyze_batch(self, recursive: bool = False) -> Dict:
        """
        Analyze all simulation runs in batch
        
        Args:
            recursive: Search recursively for result directories
        
        Returns:
            Batch analysis results
        """
        # Find all result directories
        result_dirs = self.find_result_directories(recursive=recursive)
        
        if not result_dirs:
            logger.warning("No result directories found!")
            return {'runs': [], 'summary': {}}
        
        # Analyze each run
        run_analyses = []
        for i, results_dir in enumerate(result_dirs, 1):
            logger.info(f"\n[{i}/{len(result_dirs)}] Processing: {results_dir.name}")
            analysis = self.analyze_single_run(results_dir)
            run_analyses.append(analysis)
        
        # Aggregate results
        logger.info("\n" + "=" * 70)
        logger.info("AGGREGATING RESULTS")
        logger.info("=" * 70)
        
        batch_results = {
            'runs': run_analyses,
            'summary': self._compute_batch_summary(run_analyses),
            'timestamp': datetime.now().isoformat()
        }
        
        # Save batch results
        self._save_batch_results(batch_results)
        
        # Generate batch report
        self._generate_batch_report(batch_results)
        
        return batch_results
    
    def _compute_batch_summary(self, run_analyses: List[Dict]) -> Dict:
        """Compute summary statistics across all runs"""
        logger.info("Computing batch summary...")
        
        # Aggregate by scenario
        by_scenario = defaultdict(list)
        for run in run_analyses:
            if 'error' not in run:
                scenario = run.get('scenario', 'unknown')
                by_scenario[scenario].append(run)
        
        # Compute statistics per scenario
        scenario_stats = {}
        for scenario, runs in by_scenario.items():
            total_molecules = sum(len(r.get('molecules', [])) for r in runs)
            total_instances = sum(
                r.get('statistics', {}).get('total_instances', 0) 
                for r in runs
            )
            
            # Collect all formulas
            all_formulas = set()
            formula_counts = defaultdict(int)
            for run in runs:
                for mol in run.get('molecules', []):
                    formula = mol.get('formula', 'UNKNOWN')
                    all_formulas.add(formula)
                    formula_counts[formula] += mol.get('count', 1)
            
            scenario_stats[scenario] = {
                'n_runs': len(runs),
                'total_unique_molecules': total_molecules,
                'total_instances': total_instances,
                'unique_formulas': len(all_formulas),
                'top_molecules': dict(sorted(
                    formula_counts.items(),
                    key=lambda x: x[1],
                    reverse=True
                )[:20])
            }
        
        summary = {
            'total_runs': len(run_analyses),
            'successful_runs': len([r for r in run_analyses if 'error' not in r]),
            'failed_runs': len([r for r in run_analyses if 'error' in r]),
            'scenarios': list(by_scenario.keys()),
            'by_scenario': scenario_stats,
            'overall': {
                'total_unique_molecules': sum(len(r.get('molecules', [])) for r in run_analyses if 'error' not in r),
                'total_unique_formulas': len(set(
                    mol.get('formula', 'UNKNOWN')
                    for r in run_analyses if 'error' not in r
                    for mol in r.get('molecules', [])
                ))
            }
        }
        
        logger.info(f"  Total runs: {summary['total_runs']}")
        logger.info(f"  Successful: {summary['successful_runs']}")
        logger.info(f"  Scenarios: {', '.join(summary['scenarios'])}")
        logger.info(f"  Total unique molecules: {summary['overall']['total_unique_molecules']}")
        logger.info(f"  Total unique formulas: {summary['overall']['total_unique_formulas']}")
        
        return summary
    
    def _save_batch_results(self, batch_results: Dict):
        """Save batch results to JSON"""
        output_file = self.output_dir / "batch_analysis.json"
        
        logger.info(f"Saving batch results to: {output_file}")
        
        with open(output_file, 'w') as f:
            json.dump(batch_results, f, indent=2)
    
    def _generate_batch_report(self, batch_results: Dict):
        """Generate human-readable batch report"""
        output_file = self.output_dir / "batch_report.txt"
        
        logger.info(f"Generating batch report: {output_file}")
        
        summary = batch_results['summary']
        runs = batch_results['runs']
        
        lines = []
        lines.append("=" * 70)
        lines.append("PHASE 2 BATCH ANALYSIS REPORT")
        lines.append("=" * 70)
        lines.append(f"Generated: {batch_results['timestamp']}")
        lines.append(f"Input directory: {self.input_dir}")
        lines.append(f"Output directory: {self.output_dir}")
        lines.append("")
        
        # Overall summary
        lines.append("OVERALL SUMMARY")
        lines.append("-" * 70)
        lines.append(f"Total runs analyzed: {summary['total_runs']}")
        lines.append(f"  Successful: {summary['successful_runs']}")
        lines.append(f"  Failed: {summary['failed_runs']}")
        lines.append(f"Scenarios: {', '.join(summary['scenarios'])}")
        lines.append(f"Total unique molecules: {summary['overall']['total_unique_molecules']}")
        lines.append(f"Total unique formulas: {summary['overall']['total_unique_formulas']}")
        lines.append("")
        
        # Per-scenario statistics
        lines.append("BY SCENARIO")
        lines.append("-" * 70)
        for scenario, stats in summary.get('by_scenario', {}).items():
            lines.append(f"\n{scenario.upper()}:")
            lines.append(f"  Runs: {stats['n_runs']}")
            lines.append(f"  Unique molecules: {stats['total_unique_molecules']}")
            lines.append(f"  Unique formulas: {stats['unique_formulas']}")
            lines.append(f"  Total instances: {stats['total_instances']}")
            lines.append(f"\n  Top 10 molecules:")
            for i, (formula, count) in enumerate(list(stats['top_molecules'].items())[:10], 1):
                lines.append(f"    {i:2d}. {formula:15s} : {count:6d}")
        lines.append("")
        
        # Individual runs
        lines.append("INDIVIDUAL RUNS")
        lines.append("-" * 70)
        for i, run in enumerate(runs, 1):
            lines.append(f"\n{i}. {run['run_name']}")
            if 'error' in run:
                lines.append(f"   ERROR: {run['error']}")
            else:
                lines.append(f"   Scenario: {run.get('scenario', 'unknown')}")
                lines.append(f"   Molecules: {len(run.get('molecules', []))}")
                stats = run.get('statistics', {})
                lines.append(f"   Instances: {stats.get('total_instances', 0)}")
        
        lines.append("")
        lines.append("=" * 70)
        
        # Write report
        report_text = "\n".join(lines)
        
        with open(output_file, 'w') as f:
            f.write(report_text)
        
        logger.info(f"Report saved to: {output_file}")
        
        # Also log
        logger.info("\n" + report_text)
    
    def run_matcher_on_batch(self):
        """
        Run MatcherV2 on all extracted molecules
        
        NOTE: Requires MatcherV2 to be available
        """
        if not self.use_matcher:
            logger.info("MatcherV2 disabled, skipping")
            return
        
        logger.info("\n" + "=" * 70)
        logger.info("RUNNING MATCHERV2 ON BATCH")
        logger.info("=" * 70)
        
        try:
            # Import MatcherV2
            from matcher.matcher_v2 import MatcherV2
            
            matcher = MatcherV2()
            
            # Find all molecules_for_matcher.json files
            matcher_files = list(self.input_dir.rglob("molecules_for_matcher.json"))
            logger.info(f"Found {len(matcher_files)} molecule files for matching")
            
            all_matches = []
            for i, mol_file in enumerate(matcher_files, 1):
                logger.info(f"\n[{i}/{len(matcher_files)}] Matching: {mol_file.parent.name}")
                
                # Load molecules
                with open(mol_file, 'r') as f:
                    molecules = json.load(f)
                
                # Match each molecule
                for mol in molecules:
                    try:
                        matches = matcher.match_cluster_to_pubchem(mol, top_n=5)
                        all_matches.append({
                            'source': str(mol_file.parent),
                            'molecule': mol,
                            'matches': matches
                        })
                    except Exception as e:
                        logger.error(f"Failed to match {mol.get('formula', 'UNKNOWN')}: {e}")
            
            # Save all matches
            matches_file = self.output_dir / "pubchem_matches.json"
            with open(matches_file, 'w') as f:
                json.dump(all_matches, f, indent=2)
            
            logger.info(f"\nAll matches saved to: {matches_file}")
            logger.info(f"Total molecules matched: {len(all_matches)}")
        
        except ImportError:
            logger.warning("MatcherV2 not available - skipping PubChem matching")
        except Exception as e:
            logger.error(f"Failed to run MatcherV2: {e}")


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(description="Analyze Phase 2 batch results")
    parser.add_argument('--input', required=True,
                       help='Input directory containing simulation results')
    parser.add_argument('--output', required=True,
                       help='Output directory for analysis')
    parser.add_argument('--recursive', action='store_true',
                       help='Search recursively for result directories')
    parser.add_argument('--use-matcher', action='store_true',
                       help='Use MatcherV2 for PubChem matching')
    parser.add_argument('--no-matcher', action='store_true',
                       help='Skip MatcherV2 (faster)')
    
    args = parser.parse_args()
    
    # Determine matcher usage
    use_matcher = args.use_matcher and not args.no_matcher
    
    # Create analyzer
    analyzer = Phase2BatchAnalyzer(
        input_dir=args.input,
        output_dir=args.output,
        use_matcher=use_matcher
    )
    
    # Run batch analysis
    try:
        batch_results = analyzer.analyze_batch(recursive=args.recursive)
        
        # Run MatcherV2 if requested
        if use_matcher:
            analyzer.run_matcher_on_batch()
        
        logger.info("\n" + "=" * 70)
        logger.info("[SUCCESS] Batch analysis complete!")
        logger.info(f"  Results: {args.output}/batch_analysis.json")
        logger.info(f"  Report: {args.output}/batch_report.txt")
        logger.info("=" * 70)
        
        sys.exit(0)
    
    except Exception as e:
        logger.error(f"\n[FAILED] Batch analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

