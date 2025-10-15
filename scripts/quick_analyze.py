"""
Quick Analysis Helper
======================

Quick wrapper for analyzing Phase 2 results with smart defaults.

Usage:
    # Analyze specific run
    python scripts/quick_analyze.py results/overnight_test_2025-10-13_18-17-09
    
    # Analyze all runs in directory
    python scripts/quick_analyze.py results/ --recursive
    
    # Full analysis with PubChem matching
    python scripts/quick_analyze.py results/ --full
"""

import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.sim.molecule_extractor import extract_molecules_from_results

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def quick_analyze(result_dir: Path, output_dir: Path = None, use_matcher: bool = False):
    """
    Quick analysis of a single result directory
    
    Args:
        result_dir: Directory containing simulation results
        output_dir: Output directory (default: result_dir/analysis)
        use_matcher: Use PubChem matching (slower but more detailed)
    """
    logger.info("=" * 70)
    logger.info("QUICK ANALYSIS")
    logger.info("=" * 70)
    logger.info(f"Input: {result_dir}")
    
    # Default output directory
    if output_dir is None:
        output_dir = result_dir / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Output: {output_dir}")
    logger.info(f"PubChem matching: {use_matcher}")
    logger.info("=" * 70)
    
    # Extract molecules
    logger.info("\n[1/3] Extracting molecules from results...")
    try:
        molecules = extract_molecules_from_results(result_dir)
        logger.info(f"Found {len(molecules)} unique molecules")
        
        # Save molecule list
        mol_file = output_dir / "molecules.txt"
        with open(mol_file, 'w', encoding='utf-8') as f:
            f.write(f"Molecules found: {len(molecules)}\n")
            f.write(f"Analysis time: {datetime.now()}\n")
            f.write("=" * 70 + "\n\n")
            
            for i, mol in enumerate(molecules, 1):
                f.write(f"{i}. {mol.get('formula', 'Unknown')}\n")
                if 'name' in mol:
                    f.write(f"   Name: {mol['name']}\n")
                if 'atoms' in mol:
                    f.write(f"   Atoms: {len(mol['atoms'])}\n")
                f.write("\n")
        
        logger.info(f"Saved molecule list to: {mol_file}")
        
    except Exception as e:
        logger.error(f"Failed to extract molecules: {e}")
        return False
    
    # Optional: PubChem matching
    if use_matcher and molecules:
        logger.info("\n[2/3] Matching with PubChem (this may take a while)...")
        try:
            from matcher.matcher import MatcherV2
            from matcher.compose import create_molecule_from_atoms
            
            matcher = MatcherV2()
            matches = []
            
            for mol in molecules[:20]:  # Top 20 only to save time
                if 'atoms' in mol and 'bonds' in mol:
                    try:
                        mol_obj = create_molecule_from_atoms(
                            mol['atoms'],
                            mol['bonds']
                        )
                        match = matcher.match(mol_obj)
                        if match:
                            matches.append({
                                'formula': mol.get('formula', 'Unknown'),
                                'pubchem_name': match.get('name', 'Unknown'),
                                'confidence': match.get('confidence', 0.0)
                            })
                    except Exception as e:
                        logger.debug(f"Failed to match molecule: {e}")
            
            # Save matches
            match_file = output_dir / "pubchem_matches.txt"
            with open(match_file, 'w', encoding='utf-8') as f:
                f.write(f"PubChem matches: {len(matches)}\n")
                f.write("=" * 70 + "\n\n")
                
                for i, match in enumerate(matches, 1):
                    f.write(f"{i}. {match['formula']} -> {match['pubchem_name']}\n")
                    f.write(f"   Confidence: {match['confidence']:.2%}\n\n")
            
            logger.info(f"Saved {len(matches)} PubChem matches to: {match_file}")
            
        except Exception as e:
            logger.error(f"PubChem matching failed: {e}")
    
    # Summary
    logger.info("\n[3/3] Generating summary...")
    summary_file = output_dir / "summary.txt"
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("QUICK ANALYSIS SUMMARY\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Result directory: {result_dir}\n")
        f.write(f"Analysis time: {datetime.now()}\n\n")
        f.write(f"Unique molecules found: {len(molecules)}\n")
        if use_matcher:
            f.write(f"PubChem matches: {len(matches)}\n")
        f.write("\n" + "=" * 70 + "\n")
    
    logger.info(f"Summary saved to: {summary_file}")
    logger.info("\n" + "=" * 70)
    logger.info("ANALYSIS COMPLETE!")
    logger.info("=" * 70)
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Quick analysis of Phase 2 simulation results"
    )
    parser.add_argument(
        'result_dir',
        type=str,
        help="Directory containing simulation results"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        help="Output directory (default: result_dir/analysis)"
    )
    parser.add_argument(
        '--full',
        action='store_true',
        help="Full analysis with PubChem matching (slower)"
    )
    parser.add_argument(
        '--recursive',
        action='store_true',
        help="Analyze all subdirectories (batch mode)"
    )
    
    args = parser.parse_args()
    
    result_dir = Path(args.result_dir)
    output_dir = Path(args.output) if args.output else None
    
    if not result_dir.exists():
        logger.error(f"Result directory not found: {result_dir}")
        return 1
    
    # Batch mode
    if args.recursive:
        logger.info("Batch mode: analyzing all subdirectories...")
        result_dirs = [d for d in result_dir.rglob('*') if d.is_dir() and (d / 'results.json').exists()]
        logger.info(f"Found {len(result_dirs)} result directories")
        
        for i, rd in enumerate(result_dirs, 1):
            logger.info(f"\n[{i}/{len(result_dirs)}] Analyzing: {rd.name}")
            quick_analyze(rd, use_matcher=args.full)
    else:
        # Single directory
        quick_analyze(result_dir, output_dir, use_matcher=args.full)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

