#!/usr/bin/env python3
"""
Extract Real Molecular Formulas from Snapshots
==============================================

Uses EnhancedMoleculeExtractor to extract molecules with real chemical formulas
from snapshots (using atom masses to infer types).

Usage:
    python scripts/extract_real_molecular_formulas.py \
        --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
        --output molecules_with_formulas.json
"""

import sys
import argparse
import json
import logging
from pathlib import Path
from collections import Counter

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from backend.sim.molecule_extractor_enhanced import (
    EnhancedMoleculeExtractor,
    infer_atom_type_from_mass,
    calculate_molecular_formula
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def extract_molecules_from_snapshot_with_formulas(snapshot_file: Path) -> list:
    """Extract molecules with real formulas from a single snapshot"""
    try:
        with open(snapshot_file, 'r') as f:
            data = json.load(f)
        
        bonds = data.get('bonds', [])
        attributes = data.get('attributes', [])
        positions = data.get('positions', [])
        
        if not bonds or not attributes:
            return []
        
        # Build graph of connected particles
        from collections import defaultdict
        graph = defaultdict(set)
        for bond in bonds:
            if len(bond) >= 2:
                i, j = bond[0], bond[1]
                graph[i].add(j)
                graph[j].add(i)
        
        # Find connected components (molecules)
        visited = set()
        components = []
        
        for node in graph:
            if node not in visited:
                component = []
                stack = [node]
                while stack:
                    n = stack.pop()
                    if n not in visited:
                        visited.add(n)
                        component.append(n)
                        stack.extend(graph[n])
                if len(component) >= 2:  # At least 2 atoms for a molecule
                    components.append(component)
        
        # Extract atom types and calculate formulas
        molecules = []
        for component in components:
            # Get atom types from attributes (mass is first element)
            atoms = []
            for idx in component:
                if idx < len(attributes):
                    attr = attributes[idx]
                    if isinstance(attr, list) and len(attr) > 0:
                        mass = attr[0]  # First element is mass
                        atom_type = infer_atom_type_from_mass(mass)
                        atoms.append(atom_type)
                    else:
                        atoms.append('X')  # Unknown
            
            # Calculate formula
            if atoms and not all(a == 'X' for a in atoms):
                formula = calculate_molecular_formula(atoms)
                if formula and formula != 'X':
                    molecules.append({
                        'formula': formula,
                        'atoms': atoms,
                        'size': len(component),
                        'num_atoms': len(component),
                        'bonds': len([b for b in bonds if len(b) >= 2 and b[0] in component and b[1] in component])
                    })
        
        return molecules
        
    except Exception as e:
        logger.error(f"Error processing {snapshot_file.name}: {e}")
        return []


def extract_all_molecules_with_formulas(results_dir: Path) -> list:
    """Extract all molecules with real formulas from all snapshots"""
    snapshot_dir = results_dir / "snapshots"
    if not snapshot_dir.exists():
        logger.error(f"Snapshots directory not found: {snapshot_dir}")
        return []
    
    snapshot_files = sorted(snapshot_dir.glob("step_*.json"))
    if not snapshot_files:
        logger.error(f"No snapshot files found")
        return []
    
    logger.info(f"Processing {len(snapshot_files)} snapshots...")
    
    all_molecules = []
    for snapshot_file in snapshot_files:
        molecules = extract_molecules_from_snapshot_with_formulas(snapshot_file)
        all_molecules.extend(molecules)
    
    # Aggregate by formula
    formula_counts = Counter(m['formula'] for m in all_molecules)
    
    # Get unique molecules with counts
    unique_molecules = {}
    for mol in all_molecules:
        formula = mol['formula']
        if formula not in unique_molecules:
            unique_molecules[formula] = {
                'formula': formula,
                'atoms': mol['atoms'],
                'num_atoms': mol['num_atoms'],
                'count': formula_counts[formula],
                'abundance': formula_counts[formula]
            }
    
    molecules_list = sorted(unique_molecules.values(), key=lambda x: x['count'], reverse=True)
    
    logger.info(f"Extracted {len(molecules_list)} unique molecules with real formulas")
    logger.info(f"Top 5: {[m['formula'] for m in molecules_list[:5]]}")
    
    return molecules_list


def main():
    parser = argparse.ArgumentParser(
        description="Extract real molecular formulas from snapshots"
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help='Path to simulation results directory'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='molecules_with_formulas.json',
        help='Output JSON file'
    )
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    output_file = Path(args.output)
    
    if not results_dir.exists():
        logger.error(f"Results directory not found: {results_dir}")
        return 1
    
    logger.info("="*70)
    logger.info("EXTRACTING REAL MOLECULAR FORMULAS FROM SNAPSHOTS")
    logger.info("="*70)
    logger.info(f"Results directory: {results_dir}")
    logger.info(f"Output file: {output_file}")
    logger.info("="*70)
    
    # Extract molecules
    molecules = extract_all_molecules_with_formulas(results_dir)
    
    if not molecules:
        logger.error("No molecules extracted!")
        return 1
    
    # Save to file
    with open(output_file, 'w') as f:
        json.dump(molecules, f, indent=2)
    
    logger.info(f"\n✅ Saved {len(molecules)} molecules to: {output_file}")
    logger.info(f"\nTop 10 molecules:")
    for i, mol in enumerate(molecules[:10], 1):
        logger.info(f"  {i}. {mol['formula']} (count: {mol['count']}, atoms: {mol['atoms']})")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

