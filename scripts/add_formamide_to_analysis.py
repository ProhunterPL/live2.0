"""
Add Formamide Extended results to Phase 2B analysis

This script can be run after formamide_extended runs complete to:
1. Re-run complete analysis (will automatically include formamide)
2. Or update only formamide results if other scenarios already analyzed

Usage:
    python scripts/add_formamide_to_analysis.py
    python scripts/add_formamide_to_analysis.py --update-only  # Only update formamide
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.analyze_phase2b_complete import Phase2BAnalyzer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def check_formamide_ready(input_dir: Path) -> bool:
    """Check if formamide_extended results are ready"""
    formamide_dir = input_dir / "formamide_extended"
    
    if not formamide_dir.exists():
        logger.warning(f"Formamide directory not found: {formamide_dir}")
        return False
    
    # Count completed runs (have results.json)
    run_dirs = list(formamide_dir.glob("run_*"))
    completed = sum(1 for run_dir in run_dirs if (run_dir / "results.json").exists())
    
    logger.info(f"Formamide Extended: {completed} runs found in {formamide_dir}")
    
    if completed == 0:
        logger.warning("No completed runs found. Make sure simulations have finished.")
        return False
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Add Formamide Extended to Phase 2B analysis"
    )
    parser.add_argument(
        '--input',
        default='results/phase2b_additional',
        help='Input directory (default: results/phase2b_additional)'
    )
    parser.add_argument(
        '--output',
        default='paper/results_data',
        help='Output directory (default: paper/results_data)'
    )
    parser.add_argument(
        '--update-only',
        action='store_true',
        help='Only update formamide (assumes other scenarios already analyzed)'
    )
    
    args = parser.parse_args()
    
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    
    logger.info("="*80)
    logger.info("ADDING FORMAMIDE EXTENDED TO PHASE 2B ANALYSIS")
    logger.info("="*80)
    logger.info("")
    
    # Check if formamide is ready
    if not check_formamide_ready(input_dir):
        logger.error("Formamide Extended results not ready. Exiting.")
        return 1
    
    # Run analysis
    logger.info("Running complete Phase 2B analysis (includes all scenarios)...")
    logger.info("")
    
    analyzer = Phase2BAnalyzer(input_dir, output_dir)
    analyzer.run_complete_analysis()
    analyzer.print_summary()
    
    logger.info("")
    logger.info("="*80)
    logger.info("FORMAMIDE EXTENDED ADDED TO ANALYSIS!")
    logger.info(f"Results updated in: {output_dir}")
    logger.info("="*80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

