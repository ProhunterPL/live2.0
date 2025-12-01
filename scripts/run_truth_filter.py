#!/usr/bin/env python3
"""
Truth-Filter: CLI script for filtering single run
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
        description="Truth-Filter: Validate simulation results"
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to run directory (e.g., results/phase2b_additional/miller_urey_extended/run_1)'
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
    parser.add_argument(
        '--format',
        type=str,
        choices=['json', 'markdown', 'both'],
        default='both',
        help='Report format (default: both)'
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
    
    # Filter run
    logger.info(f"Filtering run: {input_path}")
    try:
        result = filter.filter_run(str(input_path))
        
        # Generate report
        logger.info("Generating report...")
        output_files = filter.generate_report(result, format=args.format)
        
        # Print summary
        print("\n" + "=" * 70)
        print("TRUTH-FILTER VALIDATION RESULT")
        print("=" * 70)
        print(f"Run: {result.run_id}")
        print(f"Validation Level: {result.validation_level.value}")
        print(f"Overall Status: {result.overall_status.value}")
        print()
        
        print("Filter Results:")
        for filter_name, filter_result in result.filters.items():
            status_icon = {
                'PASS': '✅',
                'WARNING': '⚠️',
                'FAIL': '❌'
            }
            icon = status_icon.get(filter_result.status.value, '?')
            print(f"  {icon} {filter_name}: {filter_result.status.value}")
        
        print()
        print("Summary:")
        print(f"  Total Molecules: {result.summary.get('total_molecules', 0)}")
        print(f"  Filtered Molecules: {result.summary.get('filtered_molecules', 0)}")
        print(f"  Retention Rate: {result.summary.get('retention_rate', 0):.1%}")
        
        if result.get_warnings():
            print()
            print("Warnings:")
            for warning in result.get_warnings():
                print(f"  ⚠️ {warning}")
        
        if result.get_errors():
            print()
            print("Errors:")
            for error in result.get_errors():
                print(f"  ❌ {error}")
        
        print()
        print("Report Files:")
        for fmt, path in output_files.items():
            print(f"  {fmt}: {path}")
        print("=" * 70)
        
        # Exit code based on status
        if result.overall_status.value == 'FAIL':
            sys.exit(1)
        elif result.overall_status.value == 'WARNING':
            sys.exit(2)
        else:
            sys.exit(0)
    
    except Exception as e:
        logger.error(f"Error filtering run: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

