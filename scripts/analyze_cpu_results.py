#!/usr/bin/env python3
"""
Analyze CPU Run Results
========================

Quick analysis of completed CPU simulation.
"""

import json
from pathlib import Path

results_dir = Path("results/phase2b_local/miller_urey/cpu_run_01")

print("=" * 70)
print("CPU SIMULATION RESULTS ANALYSIS")
print("=" * 70)

# Load results
results_file = results_dir / "results.json"
if results_file.exists():
    with open(results_file) as f:
        results = json.load(f)
    
    print(f"\nScenario: {results['scenario']}")
    print(f"Description: {results['description']}")
    
    print(f"\nConfiguration:")
    print(f"  Steps: {results['configuration']['max_steps']:,}")
    print(f"  Temperature: {results['configuration']['temperature']}K")
    print(f"  Seed: {results['configuration']['seed']}")
    
    print(f"\nInitialization:")
    print(f"  Molecules placed: {results['initialization']['molecules_placed']}")
    print(f"  Atoms placed: {results['initialization']['atoms_placed']}")
    
    print(f"\nFinal State:")
    print(f"  Particles: {results['final_state']['n_particles']}")
    print(f"  Simulation time: {results['final_state']['time']:.2f}")
    print(f"  Steps completed: {results['final_state']['step']:,}")
    
    print(f"\nMolecules Detected:")
    print(f"  Total: {len(results.get('molecules_detected', []))}")
    print(f"  Novel: {len(results.get('novel_molecules', []))}")
    print(f"  (Note: 0 is expected - novelty detection was disabled)")
    
    print(f"\nSnapshots:")
    snapshots_dir = results_dir / "snapshots"
    if snapshots_dir.exists():
        snapshots = list(snapshots_dir.glob("*.json"))
        print(f"  Total snapshots: {len(snapshots)}")
        print(f"  Ready for offline batch analysis!")
    
    print(f"\nCheckpoints:")
    checkpoints_dir = results_dir / "checkpoints"
    if checkpoints_dir.exists():
        checkpoints = list(checkpoints_dir.glob("*.json"))
        print(f"  Total checkpoints: {len(checkpoints)}")
    
    print("\n" + "=" * 70)
    print("NEXT STEPS:")
    print("=" * 70)
    print("1. Run offline batch analysis on snapshots:")
    print("   python scripts/post_detect_batch.py --dir results/phase2b_local/miller_urey/cpu_run_01/snapshots --parallel 16")
    print("\n2. Or analyze manually using existing analysis scripts")
    print("\n3. Run more simulations:")
    print("   python run_phase2_cpu_test.py --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml --output results/phase2b_local/miller_urey/cpu_run_02 --steps 500000 --seed 101")
    print("=" * 70)

else:
    print(f"Results file not found: {results_file}")

