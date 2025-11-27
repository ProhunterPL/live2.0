#!/usr/bin/env python3
"""
Add Abundance History to Existing Reaction Networks
===================================================

Updates existing reaction_network.json files with abundance_history
extracted from snapshots. This enables autocatalysis detection.

Usage:
    # Single run
    python scripts/add_abundance_history.py \
        --run results/phase2b_additional/hydrothermal_extended/run_1
    
    # All runs in scenario
    python scripts/add_abundance_history.py \
        --scenario hydrothermal_extended \
        --base-dir results/phase2b_additional
"""

import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Dict, List
from collections import defaultdict

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from scripts.build_reaction_network_from_snapshots import extract_molecules_from_snapshot

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def extract_abundance_history(run_dir: Path) -> Dict[str, List[float]]:
    """
    Extract abundance history from snapshots.
    
    Returns:
        {molecule_id: [abundance_at_step_0, abundance_at_step_1, ...]}
    """
    snapshot_dir = run_dir / "snapshots"
    if not snapshot_dir.exists():
        logger.warning(f"Snapshots directory not found: {snapshot_dir}")
        return {}
    
    snapshot_files = sorted(snapshot_dir.glob("step_*.json"))
    if not snapshot_files:
        logger.warning(f"No snapshot files found in {snapshot_dir}")
        return {}
    
    logger.info(f"  Processing {len(snapshot_files)} snapshots...")
    
    # Extract molecules from each snapshot
    temporal_molecules = []
    all_molecules = set()
    
    for snapshot_file in snapshot_files:
        molecules = extract_molecules_from_snapshot(snapshot_file)
        temporal_molecules.append(molecules)
        all_molecules.update(molecules.keys())
    
    # Create molecule ID mapping (must match reaction_network.json format)
    # Load existing network to get ID mapping
    network_file = run_dir / "reaction_network.json"
    if network_file.exists():
        with open(network_file) as f:
            network = json.load(f)
        
        # Build mapping: formula -> id
        formula_to_id = {}
        for mol in network.get('molecules', []):
            formula = mol.get('formula', '')
            mol_id = mol.get('id', '')
            if formula and mol_id:
                formula_to_id[formula] = mol_id
    else:
        # Fallback: create IDs like in build_reaction_network
        formula_to_id = {mol: f"MOL_{i}" for i, mol in enumerate(sorted(all_molecules))}
    
    # Build abundance history: {molecule_id: [count_at_snapshot_0, count_at_snapshot_1, ...]}
    abundance_history = {}
    
    for molecules in temporal_molecules:
        for formula, count in molecules.items():
            mol_id = formula_to_id.get(formula, formula)
            if mol_id not in abundance_history:
                abundance_history[mol_id] = []
            abundance_history[mol_id].append(float(count))
    
    # Ensure all molecules in network have history (pad with zeros if missing)
    if network_file.exists():
        with open(network_file) as f:
            network = json.load(f)
        
        n_snapshots = len(snapshot_files)
        for mol in network.get('molecules', []):
            mol_id = mol.get('id', '')
            if mol_id and mol_id not in abundance_history:
                abundance_history[mol_id] = [0.0] * n_snapshots
    
    logger.info(f"  Extracted abundance history for {len(abundance_history)} molecules")
    return abundance_history


def update_reaction_network(run_dir: Path) -> bool:
    """Update reaction_network.json with abundance_history"""
    network_file = run_dir / "reaction_network.json"
    
    if not network_file.exists():
        logger.warning(f"  Reaction network not found: {network_file}")
        return False
    
    # Load existing network
    with open(network_file) as f:
        network = json.load(f)
    
    # Extract abundance history
    abundance_history = extract_abundance_history(run_dir)
    
    if not abundance_history:
        logger.warning(f"  No abundance history extracted")
        return False
    
    # Add abundance_history to network
    network['abundance_history'] = abundance_history
    
    # Also add molecule_names mapping for autocatalysis detector
    molecule_names = {}
    for mol in network.get('molecules', []):
        mol_id = mol.get('id', '')
        formula = mol.get('formula', '')
        if mol_id and formula:
            molecule_names[mol_id] = formula
    network['molecule_names'] = molecule_names
    
    # Update metadata
    if 'metadata' not in network:
        network['metadata'] = {}
    network['metadata']['has_abundance_history'] = True
    network['metadata']['n_snapshots'] = len(abundance_history.get(list(abundance_history.keys())[0], [])) if abundance_history else 0
    
    # Save updated network
    with open(network_file, 'w') as f:
        json.dump(network, f, indent=2)
    
    logger.info(f"  [OK] Updated {network_file.name} with abundance_history")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Add abundance history to existing reaction networks"
    )
    parser.add_argument(
        '--run',
        help='Single run directory (e.g., results/phase2b_additional/hydrothermal_extended/run_1)'
    )
    parser.add_argument(
        '--scenario',
        help='Scenario name (e.g., hydrothermal_extended) - processes all runs'
    )
    parser.add_argument(
        '--base-dir',
        type=str,
        default='results/phase2b_additional',
        help='Base directory (default: results/phase2b_additional)'
    )
    
    args = parser.parse_args()
    
    if not args.run and not args.scenario:
        parser.error("Must specify either --run or --scenario")
    
    if args.run:
        # Single run
        run_dir = Path(args.run)
        if not run_dir.exists():
            logger.error(f"Run directory not found: {run_dir}")
            sys.exit(1)
        
        logger.info(f"Updating {run_dir.name}...")
        success = update_reaction_network(run_dir)
        sys.exit(0 if success else 1)
    
    else:
        # All runs in scenario
        base_dir = Path(args.base_dir)
        scenario_dir = base_dir / args.scenario
        
        if not scenario_dir.exists():
            logger.error(f"Scenario directory not found: {scenario_dir}")
            sys.exit(1)
        
        run_dirs = sorted([
            d for d in scenario_dir.iterdir() 
            if d.is_dir() and d.name.startswith('run_')
        ], key=lambda x: int(x.name.split('_')[1]))
        
        logger.info("="*70)
        logger.info("ADDING ABUNDANCE HISTORY TO REACTION NETWORKS")
        logger.info("="*70)
        logger.info(f"Scenario: {args.scenario}")
        logger.info(f"Runs: {len(run_dirs)}")
        logger.info("="*70)
        
        results = []
        for run_dir in run_dirs:
            logger.info(f"\nProcessing {run_dir.name}...")
            success = update_reaction_network(run_dir)
            results.append((run_dir.name, success))
        
        # Summary
        successful = [r for r in results if r[1]]
        failed = [r for r in results if not r[1]]
        
        logger.info("\n" + "="*70)
        logger.info("SUMMARY")
        logger.info("="*70)
        logger.info(f"Successful: {len(successful)}/{len(results)}")
        if failed:
            logger.warning(f"Failed: {len(failed)}/{len(results)}")
            for run_name, _ in failed:
                logger.warning(f"  - {run_name}")
        logger.info("="*70)
        
        sys.exit(0 if len(successful) == len(results) else 1)


if __name__ == "__main__":
    main()

