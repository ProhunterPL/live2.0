"""
Proof-of-Concept Benchmark: Adaptive vs Fixed Spatial Hash
===========================================================

Tests adaptive spatial hashing performance for patent documentation.

Usage:
    python3 scripts/poc_adaptive_hash_benchmark.py

Output:
    - Performance comparison data
    - Cell size evolution plots
    - Statistics for patent application
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import taichi as ti
import numpy as np
import time
import json
from pathlib import Path

# Initialize Taichi (CPU for PoC)
ti.init(arch=ti.cpu, cpu_max_num_threads=4)

from backend.sim.core.adaptive_spatial_hash import (
    AdaptiveSpatialHash,
    compare_fixed_vs_adaptive,
    init_adaptive_fields
)
from backend.sim.core.spatial_hash import (
    init_spatial_hash_fields,
    build_spatial_hash,
    compute_forces_spatial_simple,
    get_stats as get_fixed_stats
)

# Initialize fields
from backend.sim.core.particles import init_particle_fields
from backend.sim.core.binding import init_binding_fields

init_particle_fields()
init_binding_fields()
init_adaptive_fields()
init_spatial_hash_fields()


def create_test_system(n_particles: int, box_size: float = 256.0, 
                       bond_fraction: float = 0.1):
    """
    Create a simple test system
    
    Args:
        n_particles: Number of particles
        box_size: Simulation box size
        bond_fraction: Fraction of particles that are bonded
    
    Returns:
        (positions, attributes, active, bond_count)
    """
    # Create Taichi fields
    positions = ti.Vector.field(2, dtype=ti.f32, shape=(n_particles,))
    attributes = ti.Vector.field(4, dtype=ti.f32, shape=(n_particles,))
    active = ti.field(dtype=ti.i32, shape=(n_particles,))
    forces = ti.Vector.field(2, dtype=ti.f32, shape=(n_particles,))
    
    # Random positions
    pos_np = np.random.uniform(0, box_size, size=(n_particles, 2))
    attr_np = np.random.uniform(-1, 1, size=(n_particles, 4))
    attr_np[:, 0] = 1.0  # mass = 1.0
    
    positions.from_numpy(pos_np)
    attributes.from_numpy(attr_np)
    active.fill(1)
    
    # Estimate bond count
    bond_count = int(n_particles * bond_fraction)
    
    return positions, attributes, active, forces, bond_count


def benchmark_fixed_hash(positions, attributes, active, particle_count: int,
                        forces, box_size: float, n_iterations: int = 100):
    """Benchmark fixed spatial hash (baseline)"""
    
    print(f"\n{'='*60}")
    print(f"Benchmarking FIXED spatial hash (baseline)")
    print(f"{'='*60}")
    
    times = []
    
    for i in range(n_iterations):
        start = time.perf_counter()
        
        # Build hash (cell_size = 10.0)
        build_spatial_hash(positions, active, particle_count, box_size, box_size)
        
        # Compute forces
        compute_forces_spatial_simple(positions, attributes, active, particle_count, forces)
        
        end = time.perf_counter()
        times.append(end - start)
    
    stats = get_fixed_stats()
    
    results = {
        'method': 'fixed',
        'mean_time': np.mean(times),
        'std_time': np.std(times),
        'min_time': np.min(times),
        'max_time': np.max(times),
        'total_time': np.sum(times),
        'iterations': n_iterations,
        'stats': stats
    }
    
    print(f"Mean time: {results['mean_time']*1000:.3f} ms/iteration")
    print(f"Total time: {results['total_time']:.3f} s")
    print(f"Grid stats: {stats}")
    
    return results


def benchmark_adaptive_hash(positions, attributes, active, particle_count: int,
                           forces, bond_count: int, box_size: float, 
                           n_iterations: int = 100):
    """Benchmark adaptive spatial hash"""
    
    print(f"\n{'='*60}")
    print(f"Benchmarking ADAPTIVE spatial hash")
    print(f"{'='*60}")
    
    adaptive_hash = AdaptiveSpatialHash(box_size, box_size)
    
    times_rebuild = []
    times_forces = []
    times_total = []
    
    for i in range(n_iterations):
        # Rebuild
        start = time.perf_counter()
        cell_size = adaptive_hash.rebuild_grid(positions, active, particle_count, bond_count)
        rebuild_time = time.perf_counter() - start
        
        # Forces
        start = time.perf_counter()
        adaptive_hash.compute_forces(positions, attributes, active, particle_count, forces)
        force_time = time.perf_counter() - start
        
        times_rebuild.append(rebuild_time)
        times_forces.append(force_time)
        times_total.append(rebuild_time + force_time)
        
        # Vary bond count slightly to test adaptation
        bond_count = int(particle_count * (0.05 + 0.1 * (i / n_iterations)))
    
    stats = adaptive_hash.get_stats()
    
    results = {
        'method': 'adaptive',
        'mean_time': np.mean(times_total),
        'std_time': np.std(times_total),
        'min_time': np.min(times_total),
        'max_time': np.max(times_total),
        'total_time': np.sum(times_total),
        'mean_rebuild_time': np.mean(times_rebuild),
        'mean_force_time': np.mean(times_forces),
        'iterations': n_iterations,
        'stats': stats,
        'cell_size_evolution': adaptive_hash.get_cell_size_evolution().tolist()
    }
    
    print(f"Mean time: {results['mean_time']*1000:.3f} ms/iteration")
    print(f"  Rebuild: {results['mean_rebuild_time']*1000:.3f} ms")
    print(f"  Forces:  {results['mean_force_time']*1000:.3f} ms")
    print(f"Total time: {results['total_time']:.3f} s")
    print(f"Adaptive stats: {stats}")
    
    return results


def run_benchmark_suite():
    """Run complete benchmark suite for patent documentation"""
    
    print("\n" + "="*70)
    print("ADAPTIVE SPATIAL HASH - PROOF-OF-CONCEPT BENCHMARK")
    print("Patent Application Support")
    print("="*70)
    
    # Test configurations
    configs = [
        {'n_particles': 500, 'bond_fraction': 0.05, 'name': 'sparse_small'},
        {'n_particles': 500, 'bond_fraction': 0.15, 'name': 'dense_small'},
        {'n_particles': 1000, 'bond_fraction': 0.05, 'name': 'sparse_medium'},
        {'n_particles': 1000, 'bond_fraction': 0.15, 'name': 'dense_medium'},
    ]
    
    box_size = 256.0
    n_iterations = 50
    
    all_results = []
    
    for config in configs:
        print(f"\n{'#'*70}")
        print(f"Configuration: {config['name']}")
        print(f"  Particles: {config['n_particles']}")
        print(f"  Bond fraction: {config['bond_fraction']}")
        print(f"{'#'*70}")
        
        # Create test system
        positions, attributes, active, forces, bond_count = create_test_system(
            config['n_particles'], box_size, config['bond_fraction']
        )
        
        # Benchmark fixed
        fixed_results = benchmark_fixed_hash(
            positions, attributes, active, config['n_particles'],
            forces, box_size, n_iterations
        )
        
        # Benchmark adaptive
        adaptive_results = benchmark_adaptive_hash(
            positions, attributes, active, config['n_particles'],
            forces, bond_count, box_size, n_iterations
        )
        
        # Compute speedup
        speedup = fixed_results['mean_time'] / adaptive_results['mean_time']
        
        # Theoretical comparison
        theoretical = compare_fixed_vs_adaptive(
            config['n_particles'], bond_count, box_size
        )
        
        result = {
            'config': config,
            'fixed': fixed_results,
            'adaptive': adaptive_results,
            'speedup': speedup,
            'theoretical': theoretical
        }
        
        all_results.append(result)
        
        print(f"\n{'='*60}")
        print(f"SPEEDUP: {speedup:.3f}x")
        print(f"Theoretical speedup: {theoretical['expected_speedup']:.3f}x")
        print(f"{'='*60}")
    
    # Save results
    output_dir = Path('results/poc_adaptive_hash')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / 'benchmark_results.json'
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n✅ Results saved to: {output_file}")
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY FOR PATENT APPLICATION")
    print(f"{'='*70}")
    
    speedups = [r['speedup'] for r in all_results]
    print(f"\nSpeedup range: {np.min(speedups):.3f}x - {np.max(speedups):.3f}x")
    print(f"Average speedup: {np.mean(speedups):.3f}x")
    print(f"Speedup std dev: {np.std(speedups):.3f}x")
    
    print("\nCell size adaptation:")
    for i, result in enumerate(all_results):
        config_name = result['config']['name']
        stats = result['adaptive']['stats']
        print(f"  {config_name}:")
        print(f"    Min cell size: {stats['min_cell_size']:.2f}")
        print(f"    Max cell size: {stats['max_cell_size']:.2f}")
        print(f"    Mean cell size: {stats['mean_cell_size']:.2f}")
        print(f"    Std cell size: {stats['std_cell_size']:.2f}")
    
    return all_results


if __name__ == '__main__':
    results = run_benchmark_suite()
    
    print("\n" + "="*70)
    print("✅ BENCHMARK COMPLETE")
    print("="*70)
    print("\nData ready for patent application!")
    print("See: results/poc_adaptive_hash/benchmark_results.json")

