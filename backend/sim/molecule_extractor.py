"""
Molecule Extractor for Phase 2 Results
=======================================

Extracts detected molecules from simulation results.
Prepares them for MatcherV2 analysis.

Integrates with:
- GraphProcessor (cluster detection)
- MatcherV2 (PubChem matching)
- Results analyzer
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import numpy as np

logger = logging.getLogger(__name__)


class MoleculeExtractor:
    """Extracts molecules from simulation results"""
    
    def __init__(self, results_dir: str):
        """
        Initialize molecule extractor
        
        Args:
            results_dir: Directory containing simulation results
        """
        self.results_dir = Path(results_dir)
        self.molecules = []
        self.statistics = {}
        
        logger.info(f"MoleculeExtractor initialized for: {self.results_dir}")
    
    def extract_from_snapshots(self, snapshot_interval: int = 1) -> List[Dict]:
        """
        Extract molecules from saved snapshots
        
        Args:
            snapshot_interval: Process every Nth snapshot (1 = all)
        
        Returns:
            List of detected molecules
        """
        logger.info("Extracting molecules from snapshots...")
        
        snapshot_dir = self.results_dir / "snapshots"
        if not snapshot_dir.exists():
            logger.warning(f"No snapshots directory found: {snapshot_dir}")
            return []
        
        # Get all snapshot files
        snapshot_files = sorted(snapshot_dir.glob("step_*.json"))
        logger.info(f"  Found {len(snapshot_files)} snapshot files")
        
        # Process snapshots
        all_molecules = []
        for i, snapshot_file in enumerate(snapshot_files[::snapshot_interval]):
            if i % 10 == 0:
                logger.info(f"  Processing snapshot {i+1}/{len(snapshot_files)//snapshot_interval}")
            
            molecules = self._extract_from_snapshot(snapshot_file)
            all_molecules.extend(molecules)
        
        logger.info(f"Extracted {len(all_molecules)} molecule instances")
        
        # Deduplicate and aggregate
        self.molecules = self._aggregate_molecules(all_molecules)
        logger.info(f"  Unique molecules: {len(self.molecules)}")
        
        return self.molecules
    
    def _extract_from_snapshot(self, snapshot_file: Path) -> List[Dict]:
        """Extract molecules from a single snapshot"""
        try:
            with open(snapshot_file, 'r') as f:
                snapshot = json.load(f)
            
            # For now, return empty - would need cluster detection integration
            # In full implementation:
            # 1. Load particle positions/types
            # 2. Run cluster detection
            # 3. Convert clusters to molecules
            # 4. Return molecule dicts
            
            return []
        
        except Exception as e:
            logger.error(f"Failed to process {snapshot_file.name}: {e}")
            return []
    
    def extract_from_final_state(self, results_file: Optional[Path] = None) -> List[Dict]:
        """
        Extract molecules from final simulation state
        
        Args:
            results_file: Path to results.json (default: results_dir/results.json)
        
        Returns:
            List of detected molecules
        """
        if results_file is None:
            results_file = self.results_dir / "results.json"
        
        logger.info(f"Extracting molecules from final state: {results_file}")
        
        try:
            with open(results_file, 'r') as f:
                results = json.load(f)
            
            # Get molecules from results
            molecules = results.get('molecules_detected', [])
            
            if not molecules:
                logger.warning("No molecules found in results file")
                logger.info("  This is expected - cluster detection needs integration")
            
            self.molecules = molecules
            logger.info(f"Found {len(molecules)} molecules in final state")
            
            return molecules
        
        except FileNotFoundError:
            logger.error(f"Results file not found: {results_file}")
            return []
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in results file: {e}")
            return []
    
    def _aggregate_molecules(self, molecules: List[Dict]) -> List[Dict]:
        """
        Aggregate molecules by formula, count occurrences
        
        Args:
            molecules: List of molecule instances
        
        Returns:
            List of unique molecules with counts
        """
        # Group by formula
        formula_groups = defaultdict(list)
        for mol in molecules:
            formula = mol.get('formula', 'UNKNOWN')
            formula_groups[formula].append(mol)
        
        # Create aggregated molecules
        aggregated = []
        for formula, instances in formula_groups.items():
            # Find most complete instance
            best_instance = max(instances, 
                              key=lambda m: len(m.get('atoms', [])))
            
            # Add count
            best_instance['count'] = len(instances)
            best_instance['first_seen'] = min(m.get('step', 0) for m in instances)
            best_instance['last_seen'] = max(m.get('step', 0) for m in instances)
            
            aggregated.append(best_instance)
        
        # Sort by count (most common first)
        aggregated.sort(key=lambda m: m.get('count', 0), reverse=True)
        
        return aggregated
    
    def compute_statistics(self) -> Dict:
        """
        Compute statistics about extracted molecules
        
        Returns:
            Dictionary of statistics
        """
        if not self.molecules:
            logger.warning("No molecules to compute statistics for")
            return {}
        
        logger.info("Computing molecule statistics...")
        
        stats = {
            'total_unique': len(self.molecules),
            'total_instances': sum(m.get('count', 1) for m in self.molecules),
            'formulas': {},
            'size_distribution': defaultdict(int),
            'complexity_distribution': defaultdict(int),
        }
        
        # Analyze each molecule
        for mol in self.molecules:
            formula = mol.get('formula', 'UNKNOWN')
            count = mol.get('count', 1)
            
            # Formula frequency
            stats['formulas'][formula] = count
            
            # Size distribution (number of atoms)
            n_atoms = len(mol.get('atoms', []))
            stats['size_distribution'][n_atoms] += count
            
            # Complexity (number of bonds)
            n_bonds = len(mol.get('bonds', []))
            stats['complexity_distribution'][n_bonds] += count
        
        # Convert defaultdicts to regular dicts
        stats['size_distribution'] = dict(stats['size_distribution'])
        stats['complexity_distribution'] = dict(stats['complexity_distribution'])
        
        # Sort formulas by count
        stats['formulas'] = dict(sorted(
            stats['formulas'].items(),
            key=lambda x: x[1],
            reverse=True
        ))
        
        self.statistics = stats
        
        logger.info(f"  Unique molecules: {stats['total_unique']}")
        logger.info(f"  Total instances: {stats['total_instances']}")
        logger.info(f"  Most common: {list(stats['formulas'].keys())[:5]}")
        
        return stats
    
    def filter_by_novelty(self, min_complexity: int = 2) -> List[Dict]:
        """
        Filter molecules by novelty criteria
        
        Args:
            min_complexity: Minimum number of bonds
        
        Returns:
            List of novel molecules
        """
        logger.info(f"Filtering molecules (min_complexity={min_complexity})")
        
        novel = []
        for mol in self.molecules:
            n_bonds = len(mol.get('bonds', []))
            
            if n_bonds >= min_complexity:
                novel.append(mol)
        
        logger.info(f"  Novel molecules: {len(novel)}/{len(self.molecules)}")
        
        return novel
    
    def export_for_matcher(self, output_file: str, 
                          include_simple: bool = False) -> int:
        """
        Export molecules in format for MatcherV2
        
        Args:
            output_file: Path to output JSON file
            include_simple: Include simple molecules (H2, H2O, etc.)
        
        Returns:
            Number of molecules exported
        """
        logger.info(f"Exporting molecules for MatcherV2: {output_file}")
        
        # Filter molecules
        if include_simple:
            export_molecules = self.molecules
        else:
            # Exclude very simple molecules
            simple_formulas = {'H2', 'H2O', 'NH3', 'CH4', 'CO2', 'O2', 'N2'}
            export_molecules = [
                m for m in self.molecules 
                if m.get('formula', '') not in simple_formulas
            ]
        
        logger.info(f"  Exporting {len(export_molecules)} molecules")
        
        # Convert to MatcherV2 format
        matcher_input = []
        for mol in export_molecules:
            matcher_mol = {
                'formula': mol.get('formula', 'UNKNOWN'),
                'atoms': mol.get('atoms', []),
                'bonds': mol.get('bonds', []),
                'count': mol.get('count', 1),
                'energy': mol.get('energy', 0.0),
                'metadata': {
                    'first_seen': mol.get('first_seen', 0),
                    'last_seen': mol.get('last_seen', 0),
                    'source': str(self.results_dir)
                }
            }
            matcher_input.append(matcher_mol)
        
        # Save to file
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(matcher_input, f, indent=2)
        
        logger.info(f"Exported to: {output_path}")
        
        return len(matcher_input)
    
    def generate_report(self, output_file: str) -> None:
        """
        Generate human-readable report
        
        Args:
            output_file: Path to output text file
        """
        logger.info(f"Generating molecule report: {output_file}")
        
        # Compute statistics if not done
        if not self.statistics:
            self.compute_statistics()
        
        report_lines = []
        report_lines.append("=" * 70)
        report_lines.append("PHASE 2 MOLECULE EXTRACTION REPORT")
        report_lines.append("=" * 70)
        report_lines.append("")
        
        # Summary
        report_lines.append("SUMMARY")
        report_lines.append("-" * 70)
        report_lines.append(f"Results directory: {self.results_dir}")
        report_lines.append(f"Unique molecules: {self.statistics.get('total_unique', 0)}")
        report_lines.append(f"Total instances: {self.statistics.get('total_instances', 0)}")
        report_lines.append("")
        
        # Top molecules
        report_lines.append("TOP 20 MOLECULES (by occurrence)")
        report_lines.append("-" * 70)
        formulas = self.statistics.get('formulas', {})
        for i, (formula, count) in enumerate(list(formulas.items())[:20], 1):
            report_lines.append(f"{i:2d}. {formula:15s} : {count:6d} instances")
        report_lines.append("")
        
        # Size distribution
        report_lines.append("SIZE DISTRIBUTION (atoms per molecule)")
        report_lines.append("-" * 70)
        size_dist = self.statistics.get('size_distribution', {})
        for size in sorted(size_dist.keys()):
            count = size_dist[size]
            bar = "█" * min(int(count / max(size_dist.values()) * 50), 50)
            report_lines.append(f"{size:2d} atoms: {count:6d} {bar}")
        report_lines.append("")
        
        # Complexity distribution
        report_lines.append("COMPLEXITY DISTRIBUTION (bonds per molecule)")
        report_lines.append("-" * 70)
        comp_dist = self.statistics.get('complexity_distribution', {})
        for complexity in sorted(comp_dist.keys()):
            count = comp_dist[complexity]
            bar = "█" * min(int(count / max(comp_dist.values()) * 50), 50)
            report_lines.append(f"{complexity:2d} bonds: {count:6d} {bar}")
        report_lines.append("")
        
        report_lines.append("=" * 70)
        
        # Write report
        report_text = "\n".join(report_lines)
        
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(report_text)
        
        logger.info(f"Report saved to: {output_path}")
        
        # Also log summary
        logger.info("\n" + report_text)


def extract_molecules_from_results(results_dir: str,
                                   output_dir: Optional[str] = None,
                                   export_for_matcher: bool = True) -> Dict:
    """
    Convenience function to extract molecules from results
    
    Args:
        results_dir: Directory containing simulation results
        output_dir: Output directory (default: results_dir/analysis)
        export_for_matcher: Export molecules for MatcherV2
    
    Returns:
        Dictionary with extraction results
    """
    if output_dir is None:
        output_dir = str(Path(results_dir) / "analysis")
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create extractor
    extractor = MoleculeExtractor(results_dir)
    
    # Extract from final state
    molecules = extractor.extract_from_final_state()
    
    # Compute statistics
    stats = extractor.compute_statistics()
    
    # Generate report
    report_file = Path(output_dir) / "molecule_report.txt"
    extractor.generate_report(str(report_file))
    
    # Export for MatcherV2
    if export_for_matcher and molecules:
        matcher_file = Path(output_dir) / "molecules_for_matcher.json"
        n_exported = extractor.export_for_matcher(str(matcher_file))
        logger.info(f"Exported {n_exported} molecules for MatcherV2")
    
    return {
        'molecules': molecules,
        'statistics': stats,
        'report_file': str(report_file),
        'output_dir': output_dir
    }


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python molecule_extractor.py <results_dir>")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    
    logging.basicConfig(level=logging.INFO)
    
    print(f"Extracting molecules from: {results_dir}")
    result = extract_molecules_from_results(results_dir)
    
    print(f"\nExtraction complete!")
    print(f"  Molecules: {len(result['molecules'])}")
    print(f"  Report: {result['report_file']}")

