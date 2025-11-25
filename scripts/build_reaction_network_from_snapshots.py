#!/usr/bin/env python3
"""
Build Reaction Network from Snapshots
=====================================

Generates reaction_network.json from temporal snapshot analysis.
This enables autocatalysis detection when reaction networks weren't saved during simulation.

Usage:
    python scripts/build_reaction_network_from_snapshots.py \
        --run results/phase2b_additional/hydrothermal_extended/run_1
"""

import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict, Counter
import networkx as nx

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def extract_molecules_from_snapshot(snapshot_file: Path) -> Dict[str, int]:
    """Extract molecules and their counts from a snapshot"""
    try:
        with open(snapshot_file) as f:
            data = json.load(f)
        
        molecules = {}
        
        # Get bonds to identify connected components (molecules)
        bonds = data.get('bonds', [])
        attributes = data.get('attributes', [])
        
        if not bonds:
            return molecules
        
        # Build graph of connected particles
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
        
        # Create molecule signatures (simplified: use size + bond pattern)
        for component in components:
            # Create a simple signature based on component size
            # In real implementation, would use molecular formula
            size = len(component)
            bond_count = sum(len(graph[n]) for n in component) // 2
            
            # Create a simple formula-like identifier
            formula = f"MOL_{size}_{bond_count}"
            molecules[formula] = molecules.get(formula, 0) + 1
        
        return molecules
        
    except Exception as e:
        logger.error(f"Error processing {snapshot_file.name}: {e}")
        return {}


def infer_reactions(temporal_molecules: List[Dict[str, int]]) -> List[Dict]:
    """
    Infer reactions from temporal changes in molecule abundances.
    
    Simple heuristic: if molecule A decreases and B increases, 
    infer A -> B (or A + C -> B + D)
    """
    reactions = []
    
    if len(temporal_molecules) < 2:
        return reactions
    
    # Compare consecutive snapshots
    for i in range(len(temporal_molecules) - 1):
        prev = temporal_molecules[i]
        curr = temporal_molecules[i + 1]
        
        # Find molecules that decreased (reactants)
        decreased = {}
        for mol, count in prev.items():
            if mol in curr and curr[mol] < count:
                decreased[mol] = count - curr[mol]
            elif mol not in curr:
                decreased[mol] = count
        
        # Find molecules that increased (products)
        increased = {}
        for mol, count in curr.items():
            if mol in prev and count > prev[mol]:
                increased[mol] = count - prev[mol]
            elif mol not in prev:
                increased[mol] = count
        
        # Create reactions (simplified: 1 reactant -> 1 product)
        # In reality, would need more sophisticated matching
        for reactant, r_count in decreased.items():
            for product, p_count in increased.items():
                if reactant != product:
                    reactions.append({
                        'reactants': [reactant],
                        'products': [product],
                        'step': i * 50000,  # Approximate step
                        'count': min(r_count, p_count)
                    })
    
    return reactions


def build_reaction_network(run_dir: Path) -> Dict:
    """Build reaction network from snapshots"""
    logger.info(f"Building reaction network for {run_dir.name}")
    
    snapshot_dir = run_dir / "snapshots"
    if not snapshot_dir.exists():
        logger.error(f"Snapshots directory not found: {snapshot_dir}")
        return None
    
    snapshot_files = sorted(snapshot_dir.glob("step_*.json"))
    if not snapshot_files:
        logger.error(f"No snapshot files found in {snapshot_dir}")
        return None
    
    logger.info(f"  Found {len(snapshot_files)} snapshots")
    
    # Extract molecules from each snapshot
    temporal_molecules = []
    all_molecules = set()
    
    for snapshot_file in snapshot_files:
        molecules = extract_molecules_from_snapshot(snapshot_file)
        temporal_molecules.append(molecules)
        all_molecules.update(molecules.keys())
    
    logger.info(f"  Found {len(all_molecules)} unique molecule types")
    
    # Infer reactions from temporal changes
    reactions = infer_reactions(temporal_molecules)
    logger.info(f"  Inferred {len(reactions)} reactions")
    
    # Build network structure
    # Create molecule ID mapping
    molecule_ids = {mol: f"MOL_{i}" for i, mol in enumerate(sorted(all_molecules))}
    
    # Convert reactions to edges format (for autocatalysis detector)
    edges = []
    for rxn in reactions:
        reactants = rxn.get('reactants', [])
        products = rxn.get('products', [])
        
        # Create edges: reactant -> product
        for reactant in reactants:
            reactant_id = molecule_ids.get(reactant, reactant)
            for product in products:
                product_id = molecule_ids.get(product, product)
                edges.append({
                    'source': reactant_id,
                    'target': product_id,
                    'reaction': rxn
                })
    
    network = {
        'molecules': [
            {'formula': mol, 'id': molecule_ids[mol]} 
            for mol in sorted(all_molecules)
        ],
        'reactions': reactions,  # Keep for reference
        'edges': edges,  # For autocatalysis detector
        'metadata': {
            'source': 'snapshots',
            'n_snapshots': len(snapshot_files),
            'n_molecules': len(all_molecules),
            'n_reactions': len(reactions),
            'n_edges': len(edges)
        }
    }
    
    return network


def main():
    parser = argparse.ArgumentParser(
        description="Build reaction network from snapshots for autocatalysis detection"
    )
    parser.add_argument(
        '--run',
        required=True,
        help='Run directory (e.g., results/phase2b_additional/hydrothermal_extended/run_1)'
    )
    parser.add_argument(
        '--output',
        help='Output file (default: run_dir/reaction_network.json)'
    )
    
    args = parser.parse_args()
    
    run_dir = Path(args.run)
    if not run_dir.exists():
        logger.error(f"Run directory not found: {run_dir}")
        sys.exit(1)
    
    # Build network
    network = build_reaction_network(run_dir)
    
    if not network:
        logger.error("Failed to build reaction network")
        sys.exit(1)
    
    # Save network
    output_file = Path(args.output) if args.output else run_dir / "reaction_network.json"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(network, f, indent=2)
    
    logger.info(f"  [OK] Reaction network saved to {output_file}")
    logger.info(f"  Molecules: {network['metadata']['n_molecules']}")
    logger.info(f"  Reactions: {network['metadata']['n_reactions']}")
    
    print(f"\n[OK] Reaction network generated successfully!")
    print(f"  File: {output_file}")
    print(f"  Next step: Run autocatalysis detection with this network file")


if __name__ == "__main__":
    main()

