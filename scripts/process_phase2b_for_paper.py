"""
Master Script: Process Phase 2B Results for Paper 1

Runs complete pipeline:
1. Analysis (autocatalysis + complexity)
2. Figure generation
3. Table generation

ONE COMMAND to go from AWS results → Paper-ready outputs!

Usage:
    python scripts/process_phase2b_for_paper.py \
        --input results/phase2b_additional

Author: Live 2.0 Team
Date: November 2025
"""

import argparse
import sys
import subprocess
import logging
from pathlib import Path
import time

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_command(cmd: list, description: str) -> bool:
    """Run command and handle errors"""
    logger.info(f"\n{'='*80}")
    logger.info(f"STEP: {description}")
    logger.info(f"{'='*80}")
    logger.info(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        if result.stdout:
            logger.info(result.stdout)
            
        logger.info(f"[OK] {description} - COMPLETE")
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"[ERROR] {description} - FAILED")
        logger.error(f"Error: {e.stderr}")
        return False


def check_input_directory(input_dir: Path) -> bool:
    """Verify input directory structure"""
    logger.info("\n" + "="*80)
    logger.info("CHECKING INPUT DIRECTORY")
    logger.info("="*80)
    
    if not input_dir.exists():
        logger.error(f"[ERROR] Input directory not found: {input_dir}")
        return False
        
    logger.info(f"[OK] Input directory exists: {input_dir}")
    
    # Check for scenario directories
    scenarios = ['miller_urey_extended', 'hydrothermal_extended', 'formamide_extended']
    found_scenarios = []
    
    for scenario in scenarios:
        scenario_dir = input_dir / scenario
        if scenario_dir.exists():
            # Check for run directories
            runs = list(scenario_dir.glob('run_*'))
            logger.info(f"[OK] {scenario}: {len(runs)} runs found")
            found_scenarios.append(scenario)
        else:
            logger.warning(f"[WARNING] {scenario}: directory not found")
            
    if len(found_scenarios) == 0:
        logger.error("[ERROR] No scenario directories found!")
        return False
        
    logger.info(f"\n[OK] Found {len(found_scenarios)}/3 scenarios")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Master script: Process Phase 2B results for paper"
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input directory (results/phase2b_additional)'
    )
    parser.add_argument(
        '--skip-analysis',
        action='store_true',
        help='Skip analysis step (if already done)'
    )
    parser.add_argument(
        '--skip-figures',
        action='store_true',
        help='Skip figure generation'
    )
    parser.add_argument(
        '--skip-tables',
        action='store_true',
        help='Skip table generation'
    )
    
    args = parser.parse_args()
    
    # Setup paths
    input_dir = Path(args.input)
    data_dir = Path('paper/results_data')
    figures_dir = Path('paper/figures')
    tables_dir = Path('paper/tables')
    
    # Header
    print("\n" + "="*80)
    print("PHASE 2B -> PAPER 1: COMPLETE PROCESSING PIPELINE")
    print("="*80)
    print(f"\nInput:   {input_dir}")
    print(f"Data:    {data_dir}")
    print(f"Figures: {figures_dir}")
    print(f"Tables:  {tables_dir}")
    print()
    
    start_time = time.time()
    
    # Check input
    if not check_input_directory(input_dir):
        logger.error("\n[ERROR] Input validation failed. Exiting.")
        sys.exit(1)
        
    success = True
    
    # Step 1: Analysis
    if not args.skip_analysis:
        cmd = [
            sys.executable,
            'scripts/analyze_phase2b_complete.py',
            '--input', str(input_dir),
            '--output', str(data_dir)
        ]
        success = run_command(cmd, "ANALYSIS (autocatalysis + complexity)")
        
        if not success:
            logger.error("\n[ERROR] Analysis failed. Cannot proceed.")
            sys.exit(1)
    else:
        logger.info("\n[SKIP] Skipping analysis (--skip-analysis)")
        
    # Step 2: Figures
    if not args.skip_figures:
        cmd = [
            sys.executable,
            'scripts/generate_all_figures.py',
            '--data', str(data_dir),
            '--output', str(figures_dir)
        ]
        success = run_command(cmd, "FIGURE GENERATION")
        
        if not success:
            logger.warning("\n[WARNING] Figure generation failed (continuing...)")
    else:
        logger.info("\n[SKIP] Skipping figures (--skip-figures)")
        
    # Step 3: Tables
    if not args.skip_tables:
        cmd = [
            sys.executable,
            'scripts/generate_all_tables.py',
            '--data', str(data_dir),
            '--output', str(tables_dir)
        ]
        success = run_command(cmd, "TABLE GENERATION")
        
        if not success:
            logger.warning("\n[WARNING] Table generation failed (continuing...)")
    else:
        logger.info("\n[SKIP] Skipping tables (--skip-tables)")
        
    # Summary
    elapsed = time.time() - start_time
    
    print("\n" + "="*80)
    print("PIPELINE COMPLETE!")
    print("="*80)
    print(f"\nElapsed time: {elapsed/60:.1f} minutes")
    print(f"\nGenerated outputs:")
    print(f"  Analysis data: {data_dir}/")
    print(f"  Figures:       {figures_dir}/")
    print(f"  Tables:        {tables_dir}/")
    
    print(f"\nNext steps:")
    print(f"  1. Review outputs in {data_dir}/")
    print(f"  2. Check figures in {figures_dir}/")
    print(f"  3. Check tables in {tables_dir}/")
    print(f"  4. Fill Results section in manuscript_draft.tex")
    print(f"  5. Fill Discussion section in manuscript_draft.tex")
    
    print("\n" + "="*80)
    print("See TIER1_IMPLEMENTATION_GUIDE.md for detailed instructions")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()

