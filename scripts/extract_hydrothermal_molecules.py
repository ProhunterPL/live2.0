#!/usr/bin/env python3
"""
Extract molecules from Phase 2B hydrothermal snapshots
"""

import sys
from pathlib import Path
import json
import logging

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from backend.sim.molecule_extractor import extract_molecules_from_results

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    base_dir = Path("results/phase2b_additional/hydrothermal_extended")
    
    if not base_dir.exists():
        logger.error(f"Directory not found: {base_dir}")
        return
    
    # Find all run directories
    run_dirs = sorted(base_dir.glob("run_*"), key=lambda x: int(x.name.split("_")[1]))
    
    logger.info(f"Found {len(run_dirs)} runs to process")
    
    for run_dir in run_dirs:
        run_id = run_dir.name
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing {run_id}")
        logger.info(f"{'='*60}")
        
        # Check if results.json exists
        results_file = run_dir / "results.json"
        if not results_file.exists():
            logger.warning(f"  Skipping {run_id}: results.json not found")
            continue
        
        # Check if molecules.json already exists and has data
        molecules_file = run_dir / "molecules.json"
        if molecules_file.exists():
            try:
                with open(molecules_file) as f:
                    molecules = json.load(f)
                if len(molecules) > 0:
                    logger.info(f"  {run_id}: molecules.json already exists with {len(molecules)} molecules")
                    continue
            except:
                pass
        
        # Extract molecules
        try:
            result = extract_molecules_from_results(
                str(run_dir),
                output_dir=str(run_dir / "analysis"),
                export_for_matcher=False
            )
            
            molecules = result['molecules']
            stats = result['statistics']
            
            # Save to molecules.json
            with open(molecules_file, 'w') as f:
                json.dump(molecules, f, indent=2)
            
            logger.info(f"  ✓ Extracted {len(molecules)} molecules")
            logger.info(f"  ✓ Statistics: {stats.get('total_molecules', 0)} total, {stats.get('unique_smiles', 0)} unique")
            
        except Exception as e:
            logger.error(f"  ✗ Error processing {run_id}: {e}")
            import traceback
            traceback.print_exc()
    
    logger.info(f"\n{'='*60}")
    logger.info("Extraction complete!")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()

