#!/usr/bin/env python3
"""
Master Figure Generator
=======================

Generates all figures for the publication in one go.

Usage:
    python scripts/generate_all_figures.py --input results/phase2_full --output figures/
"""

import sys
import argparse
import logging
import subprocess
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run_figure_script(script_name: str, input_dir: Path, output_dir: Path):
    """Run a figure generation script"""
    script_path = Path(__file__).parent / script_name
    
    if not script_path.exists():
        logger.error(f"Script not found: {script_path}")
        return False
    
    logger.info(f"Running {script_name}...")
    
    try:
        result = subprocess.run([
            sys.executable, str(script_path),
            '--input', str(input_dir),
            '--output', str(output_dir)
        ], capture_output=True, text=True, check=True)
        
        logger.info(f"✅ {script_name} completed successfully")
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ {script_name} failed:")
        logger.error(f"   stdout: {e.stdout}")
        logger.error(f"   stderr: {e.stderr}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Generate all publication figures')
    parser.add_argument('--input', type=Path, required=True, help='Results directory')
    parser.add_argument('--output', type=Path, required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("MASTER FIGURE GENERATOR")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    # List of figure scripts to run
    figure_scripts = [
        'plot_molecular_diversity.py',      # Figure 3
        'plot_reaction_networks.py',        # Figure 4
        'plot_autocatalytic_cycles.py',     # Figure 5
        'plot_top_molecules.py',            # Figure 6
        'plot_emergence_timeline.py'        # Figure 7
    ]
    
    # Run each script
    success_count = 0
    total_count = len(figure_scripts)
    
    for script in figure_scripts:
        if run_figure_script(script, args.input, args.output):
            success_count += 1
    
    # Summary
    logger.info("=" * 70)
    logger.info("FIGURE GENERATION SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Successfully generated: {success_count}/{total_count} figures")
    
    if success_count == total_count:
        logger.info("🎉 ALL FIGURES GENERATED SUCCESSFULLY!")
        return 0
    else:
        logger.warning(f"⚠️  {total_count - success_count} figures failed to generate")
        return 1

if __name__ == "__main__":
    sys.exit(main())
