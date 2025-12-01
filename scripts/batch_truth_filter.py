#!/usr/bin/env python3
"""
Truth-Filter: CLI script for batch filtering
"""
import argparse
import sys
import logging
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.validation import TruthFilter, ValidationLevel

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Truth-Filter: Batch validate simulation results"
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to results directory (e.g., results/phase2b_additional)'
    )
    parser.add_argument(
        '--scenario',
        type=str,
        help='Optional scenario filter (e.g., miller_urey_extended)'
    )
    parser.add_argument(
        '--level',
        type=str,
        choices=['STRICT', 'MEDIUM', 'LENIENT'],
        default='MEDIUM',
        help='Validation level (default: MEDIUM)'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Output directory for reports (default: validation_reports/)'
    )
    
    args = parser.parse_args()
    
    # Validate input
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input directory does not exist: {input_path}")
        sys.exit(1)
    
    # Set output directory
    if args.output:
        output_dir = args.output
    else:
        output_dir = "validation_reports"
    
    # Convert level string to enum
    validation_level = ValidationLevel[args.level]
    
    # Initialize filter
    logger.info(f"Initializing Truth-Filter with level: {validation_level.value}")
    filter = TruthFilter(
        validation_level=validation_level,
        output_dir=output_dir
    )
    
    # Filter batch
    logger.info(f"Filtering batch: {input_path}")
    if args.scenario:
        logger.info(f"  Scenario filter: {args.scenario}")
    
    try:
        results = filter.filter_batch(str(input_path), scenario=args.scenario)
        
        if not results:
            logger.warning("No runs found to filter")
            sys.exit(1)
        
        # Generate summary
        logger.info("Generating summary...")
        summary = filter.generate_summary(results)
        
        # Save summary
        summary_path = filter.report_generator.save_summary_markdown(
            summary,
            "truth_filter_summary.md"
        )
        
        # Print summary
        print("\n" + "=" * 70)
        print("TRUTH-FILTER BATCH VALIDATION SUMMARY")
        print("=" * 70)
        print(f"Total Runs: {summary['total_runs']}")
        print(f"Passed: {summary['passed_runs']} ({summary['pass_rate']:.1%})")
        print(f"Warnings: {summary['warning_runs']}")
        print(f"Failed: {summary['failed_runs']}")
        print()
        
        print("Filter Statistics:")
        for filter_name, stats in summary.get('filter_statistics', {}).items():
            print(f"  {filter_name}:")
            print(f"    Pass: {stats['pass']}")
            print(f"    Warning: {stats['warning']}")
            print(f"    Fail: {stats['fail']}")
        
        print()
        mol_stats = summary.get('molecule_statistics', {})
        if mol_stats:
            print("Molecule Statistics:")
            print(f"  Total Original: {mol_stats['total_original']}")
            print(f"  Total Filtered: {mol_stats['total_filtered']}")
            print(f"  Retention Rate: {mol_stats['retention_rate']:.1%}")
        
        print()
        print(f"Summary saved to: {summary_path}")
        print("=" * 70)
        
        # Exit code based on results
        if summary['failed_runs'] > 0:
            sys.exit(1)
        elif summary['warning_runs'] > 0:
            sys.exit(2)
        else:
            sys.exit(0)
    
    except Exception as e:
        logger.error(f"Error filtering batch: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

