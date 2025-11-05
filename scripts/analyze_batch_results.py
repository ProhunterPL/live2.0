#!/usr/bin/env python3
"""Quick analysis of batch detection results"""

import json
from pathlib import Path

results_dir = Path("results/phase2b_local/miller_urey/cpu_run_02/post_detect")

# Load all results
all_substances = []
for f in results_dir.glob("step_*.json"):
    with open(f) as file:
        data = json.load(file)
        all_substances.extend(data.get('novel_substances', []))

# Get unique substances
unique = {s['substance_id']: s for s in all_substances}

print("=" * 70)
print("BATCH ANALYSIS SUMMARY")
print("=" * 70)
print(f"Total novel substances (with duplicates): {len(all_substances)}")
print(f"Unique novel substances: {len(unique)}")
print(f"Largest cluster: {max([s['size'] for s in unique.values()])} particles")
print(f"Average cluster size: {sum([s['size'] for s in unique.values()])/len(unique):.1f} particles")
print(f"Smallest cluster: {min([s['size'] for s in unique.values()])} particles")
print("=" * 70)

# Distribution by size
size_dist = {}
for s in unique.values():
    size = s['size']
    size_dist[size] = size_dist.get(size, 0) + 1

print("\nCluster Size Distribution:")
for size in sorted(size_dist.keys()):
    print(f"  {size} particles: {size_dist[size]} substances")

