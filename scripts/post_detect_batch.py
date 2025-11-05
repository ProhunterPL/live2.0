#!/usr/bin/env python3
"""
Batch Post-Detection Analysis
=============================

Runs novelty detection offline on saved snapshots.
Can be run in parallel on multiple machines/cores.

Usage:
    # Single snapshot
    python scripts/post_detect_batch.py --input results/snapshot_50000.npz --output results/clusters_50000.json
    
    # Batch all snapshots (parallel)
    python scripts/post_detect_batch.py --dir results/phase2b_local/miller_urey/run_01 --parallel 32
    
    # With GNU parallel
    find results/ -name "state_*.npz" | \
      xargs -n1 -P32 python scripts/post_detect_batch.py --input {}
"""

import sys
import argparse
import json
import logging
from pathlib import Path
from datetime import datetime
import multiprocessing

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from backend.sim.core.binding import BindingSystem
from backend.sim.core.graphs import GraphProcessor, MolecularGraph
from backend.sim.core.catalog import SubstanceCatalog

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def process_snapshot(snapshot_path: Path, output_path: Path, config):
    """Process a single snapshot for novelty detection"""
    
    logger.info(f"Processing: {snapshot_path}")
    
    try:
        # Load snapshot (supports both JSON and NPZ formats)
        if snapshot_path.suffix == '.json':
            with open(snapshot_path, 'r') as f:
                data = json.load(f)
            
            # Extract data from JSON
            step = data.get('step', 0)
            positions = np.array(data.get('positions', []))
            attributes = np.array(data.get('attributes', []))
            bonds = data.get('bonds', [])  # List of (i, j, strength) tuples
            clusters = data.get('clusters', [])  # List of particle index lists
            
            # If clusters not in snapshot, use bonds to infer them
            if not clusters and bonds:
                # Build cluster graph from bonds
                from collections import defaultdict
                cluster_graph = defaultdict(set)
                for bond in bonds:
                    i, j = bond[0], bond[1]
                    cluster_graph[i].add(j)
                    cluster_graph[j].add(i)
                
                # Find connected components
                visited = set()
                clusters = []
                for node in cluster_graph:
                    if node not in visited:
                        cluster = []
                        stack = [node]
                        while stack:
                            n = stack.pop()
                            if n not in visited:
                                visited.add(n)
                                cluster.append(n)
                                stack.extend(cluster_graph[n])
                        if len(cluster) >= 2:
                            clusters.append(cluster)
        else:
            # NPZ format (legacy)
            data = np.load(snapshot_path)
            positions = data['positions']
            attributes = data['attributes']
            bonds = data['bonds']
            step = int(data.get('step', 0))
            
            # Detect clusters (would need BindingSystem, but for now use bonds)
            from collections import defaultdict
            cluster_graph = defaultdict(set)
            for bond in bonds:
                i, j = bond[0], bond[1]
                cluster_graph[i].add(j)
                cluster_graph[j].add(i)
            
            visited = set()
            clusters = []
            for node in cluster_graph:
                if node not in visited:
                    cluster = []
                    stack = [node]
                    while stack:
                        n = stack.pop()
                        if n not in visited:
                            visited.add(n)
                            cluster.append(n)
                            stack.extend(cluster_graph[n])
                    if len(cluster) >= 2:
                        clusters.append(cluster)
        
        # Process each cluster for novelty detection
        novel_substances = []
        catalog = SubstanceCatalog()
        
        for cluster in clusters:
            if len(cluster) >= 3:  # Minimum size for molecules
                # Get cluster bonds
                cluster_bonds = []
                bond_dict = {(b[0], b[1]): b[2] if len(b) > 2 else 1.0 for b in bonds}
                for i in cluster:
                    for j in cluster:
                        if i < j and (i, j) in bond_dict:
                            cluster_bonds.append((i, j))
                
                # Get attributes
                particle_attributes = {}
                for idx in cluster:
                    if idx < len(attributes):
                        particle_attributes[idx] = attributes[idx].tolist() if hasattr(attributes[idx], 'tolist') else attributes[idx]
                
                # Create graph
                try:
                    graph = MolecularGraph(cluster, cluster_bonds, particle_attributes)
                    
                    # Add to catalog
                    is_novel, substance_id = catalog.add_substance(graph, 0.0, {})
                    
                    if is_novel:
                        novel_substances.append({
                            'substance_id': substance_id,
                            'cluster': cluster,
                            'size': len(cluster),
                            'bonds': cluster_bonds,
                            'step': step
                        })
                except Exception as e:
                    logger.debug(f"Failed to process cluster {cluster}: {e}")
                    continue
        
        # Save results
        results = {
            'snapshot': str(snapshot_path),
            'step': step,
            'n_clusters': len(clusters),
            'n_novel': len(novel_substances),
            'novel_substances': novel_substances,
            'timestamp': datetime.now().isoformat()
        }
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Processed: {len(novel_substances)} novel substances from {len(clusters)} clusters")
        
        return results
        
    except Exception as e:
        logger.error(f"Failed to process {snapshot_path}: {e}")
        import traceback
        traceback.print_exc()
        return None


def process_directory(directory: Path, parallel: int = None):
    """Process all snapshots in a directory"""
    
    logger.info(f"Processing directory: {directory}")
    
    # Check if directory has snapshots subdirectory
    snapshots_dir = directory / "snapshots" if (directory / "snapshots").exists() else directory
    
    # Find all snapshots (both JSON and NPZ formats)
    snapshots_json = list(snapshots_dir.glob("step_*.json"))
    snapshots_npz = list(snapshots_dir.glob("**/state_*.npz"))
    snapshots = snapshots_json + snapshots_npz
    
    if not snapshots:
        logger.warning(f"No snapshots found in {directory}")
        logger.info(f"Searched in: {snapshots_dir}")
        logger.info(f"Patterns: step_*.json and state_*.npz")
        return
    
    logger.info(f"Found {len(snapshots)} snapshots")
    
    # Create output directory
    output_dir = directory / "post_detect"
    output_dir.mkdir(exist_ok=True)
    
    # Process snapshots
    results = []
    
    if parallel and parallel > 1:
        # Parallel processing
        logger.info(f"Running {parallel} parallel workers")
        
        with multiprocessing.Pool(parallel) as pool:
            tasks = [
                (snapshot, output_dir / f"{snapshot.stem}.json", None)
                for snapshot in snapshots
            ]
            results = pool.starmap(process_snapshot, tasks)
    else:
        # Sequential processing
        for snapshot in snapshots:
            output_path = output_dir / f"{snapshot.stem}.json"
            result = process_snapshot(snapshot, output_path, None)
            results.append(result)
    
    # Summary
    logger.info("=" * 70)
    logger.info(f"BATCH PROCESSING COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Total snapshots: {len(snapshots)}")
    logger.info(f"Successful: {sum(1 for r in results if r)}")
    logger.info(f"Failed: {sum(1 for r in results if not r)}")
    
    # Save summary
    summary_path = output_dir / "batch_summary.json"
    with open(summary_path, 'w') as f:
        json.dump({
            'total_snapshots': len(snapshots),
            'successful': sum(1 for r in results if r),
            'failed': sum(1 for r in results if not r),
            'results': results,
            'timestamp': datetime.now().isoformat()
        }, f, indent=2)
    
    logger.info(f"Summary saved to: {summary_path}")


def main():
    parser = argparse.ArgumentParser(description="Batch Post-Detection Analysis")
    parser.add_argument('--input', type=Path, help='Single snapshot file or directory')
    parser.add_argument('--output', type=Path, help='Output file (optional, auto-generated if not provided)')
    parser.add_argument('--dir', type=Path, help='Directory with snapshots (alternative to --input)')
    parser.add_argument('--parallel', type=int, default=1,
                       help='Number of parallel workers (default: 1)')
    
    args = parser.parse_args()
    
    # Handle --dir (explicit directory)
    if args.dir:
        process_directory(args.dir, args.parallel)
    # Handle --input
    elif args.input:
        input_path = Path(args.input)
        
        # Check if input is a directory
        if input_path.is_dir():
            process_directory(input_path, args.parallel)
        # Check if input is a file
        elif input_path.is_file():
            # Generate output path if not provided
            if args.output:
                output_path = Path(args.output)
            else:
                # Auto-generate output path next to input file
                output_dir = input_path.parent / "post_detect"
                output_dir.mkdir(exist_ok=True)
                output_path = output_dir / f"{input_path.stem}_detected.json"
            
            process_snapshot(input_path, output_path, None)
        else:
            parser.error(f"Input path does not exist: {input_path}")
    else:
        parser.error("Must specify --input (file or directory) or --dir")


if __name__ == "__main__":
    main()

