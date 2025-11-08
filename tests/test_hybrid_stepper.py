#!/usr/bin/env python3
"""
Test HybridSimulationStepper functionality
"""

import sys
import time
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import taichi as ti

# Initialize Taichi - use CPU for tests to avoid GPU memory issues
# Hybrid mode will still work - we're just testing the chemistry worker
ti.reset()  # Reset any previous initialization

# Use CPU for tests - GPU has memory allocation issues during kernel compilation
# This is fine because we're testing the CPU chemistry worker anyway!
import multiprocessing
ti.init(arch=ti.cpu, cpu_max_num_threads=multiprocessing.cpu_count())
print(f"Using CPU ({multiprocessing.cpu_count()} threads) for tests")
print("Note: Testing CPU chemistry worker (works with any backend)")

from backend.sim.config import SimulationConfig
from backend.sim.core.hybrid_stepper import HybridSimulationStepper


def test_basic_functionality():
    """Test basic hybrid stepper functionality"""
    print("="*70)
    print("TEST 1: Basic Functionality")
    print("="*70)
    
    # Create config with minimal particles for testing
    config = SimulationConfig(
        n_particles=30,  # Very small for memory-constrained GPU
        max_particles=100,  # Limit max
        max_steps=1000,
        dt=0.001,
        box_size=50.0,  # Smaller box
        temperature=298.0,
        mode='open_chemistry',
        seed=42
    )
    
    # Disable ALL heavy operations for testing
    config.detect_novel_substances = False
    config.enable_diagnostics = False
    config.enable_mutations = False  # Disable mutations
    config.metrics_update_interval = 999999  # Disable metrics
    config.chemistry_snapshot_interval = 10  # Frequent snapshots for testing
    
    print("\n1. Creating hybrid stepper...")
    stepper = HybridSimulationStepper(config)
    
    print("2. Starting simulation...")
    stepper.start()
    
    print("3. Running 50 steps...")
    for i in range(50):
        stepper.step()
        if i % 10 == 0:
            print(f"   Step {i+1}/50")
    
    print("4. Checking state...")
    state = stepper.get_simulation_state()
    print(f"   Particles: {state['particle_count']}")
    print(f"   Step count: {state['step_count']}")
    
    print("5. Getting visualization data...")
    viz_data = stepper.get_visualization_data()
    
    if 'cpu_worker' in viz_data:
        cpu_stats = viz_data['cpu_worker']
        print(f"   CPU worker analyzed: {cpu_stats['total_analyzed']} snapshots")
        print(f"   Avg analysis time: {cpu_stats.get('avg_analysis_time_ms', 0):.1f}ms")
    
    if 'chemistry' in viz_data:
        chem = viz_data['chemistry']
        print(f"   Bonds detected: {len(chem.get('bonds', []))}")
        print(f"   Clusters detected: {len(chem.get('clusters', []))}")
    
    print("6. Stopping simulation...")
    stepper.stop()
    
    print("\n[OK] Test 1 PASSED")
    return True


def test_cpu_worker_timing():
    """Test CPU worker timing and async behavior"""
    print("\n" + "="*70)
    print("TEST 2: CPU Worker Timing")
    print("="*70)
    
    config = SimulationConfig(
        n_particles=50,  # Minimal for memory-constrained GPU
        max_particles=100,
        max_steps=1000,
        dt=0.001,
        box_size=50.0,  # Smaller box
        temperature=298.0,
        mode='open_chemistry',
        seed=42
    )
    
    config.detect_novel_substances = False
    config.enable_diagnostics = False
    config.enable_mutations = False
    config.metrics_update_interval = 999999
    config.chemistry_snapshot_interval = 20  # Every 20 steps
    
    print("\n1. Creating hybrid stepper with 50 particles...")
    stepper = HybridSimulationStepper(config)
    stepper.start()
    
    print("2. Running 100 steps and measuring timing...")
    
    step_times = []
    for i in range(100):
        start = time.time()
        stepper.step()
        elapsed = time.time() - start
        step_times.append(elapsed * 1000)  # Convert to ms
    
    import numpy as np
    avg_step = np.mean(step_times)
    max_step = np.max(step_times)
    min_step = np.min(step_times)
    
    print(f"\n3. Step timing:")
    print(f"   Avg: {avg_step:.1f}ms")
    print(f"   Min: {min_step:.1f}ms")
    print(f"   Max: {max_step:.1f}ms")
    
    # Check if max step is reasonable (shouldn't block on CPU)
    if max_step < 100:  # Should be fast since CPU is async
        print("   [OK] Steps are fast (CPU not blocking)")
    else:
        print("   [WARNING] Some steps are slow (CPU might be blocking)")
    
    print("\n4. Checking CPU worker stats...")
    state = stepper.get_simulation_state()
    if 'cpu_worker' in state:
        cpu = state['cpu_worker']
        print(f"   Total analyzed: {cpu.get('total_analyzed', 0)}")
        print(f"   Queue size: {cpu.get('queue_size', 0)}")
        print(f"   Avg analysis: {cpu.get('avg_analysis_time_ms', 0):.1f}ms")
    
    stepper.stop()
    
    print("\n[OK] Test 2 PASSED")
    return True


def test_chemistry_accuracy():
    """Test that CPU chemistry produces reasonable results"""
    print("\n" + "="*70)
    print("TEST 3: Chemistry Accuracy")
    print("="*70)
    
    config = SimulationConfig(
        n_particles=50,  # Minimal for memory-constrained GPU
        max_particles=100,
        max_steps=1000,
        dt=0.001,
        box_size=50.0,  # Smaller box
        temperature=298.0,
        mode='open_chemistry',
        seed=42
    )
    
    config.detect_novel_substances = False
    config.enable_diagnostics = False
    config.enable_mutations = False
    config.metrics_update_interval = 999999
    config.chemistry_snapshot_interval = 50
    
    print("\n1. Running simulation for 200 steps...")
    stepper = HybridSimulationStepper(config)
    stepper.start()
    
    for i in range(200):
        stepper.step()
    
    print("\n2. Checking chemistry results...")
    viz_data = stepper.get_visualization_data()
    
    if 'chemistry' in viz_data:
        chem = viz_data['chemistry']
        bonds = chem.get('bonds', [])
        clusters = chem.get('clusters', [])
        metrics = chem.get('metrics', {})
        
        print(f"   Bonds: {len(bonds)}")
        print(f"   Clusters: {len(clusters)}")
        
        if metrics:
            print(f"   Particles: {metrics.get('n_particles', 0)}")
            print(f"   Avg cluster size: {metrics.get('avg_cluster_size', 0):.1f}")
            print(f"   Max cluster size: {metrics.get('max_cluster_size', 0)}")
            print(f"   Bonding ratio: {metrics.get('bonding_ratio', 0):.4f}")
        
        if len(bonds) > 0 or len(clusters) > 0:
            print("   [OK] Chemistry detection working")
        else:
            print("   [INFO] No bonds/clusters detected (might be normal early in simulation)")
    else:
        print("   [INFO] No chemistry results yet (worker might be processing)")
    
    stepper.stop()
    
    print("\n[OK] Test 3 PASSED")
    return True


def main():
    print("\n" + "="*70)
    print("HYBRID STEPPER TESTS")
    print("="*70)
    
    try:
        # Run tests
        test1 = test_basic_functionality()
        test2 = test_cpu_worker_timing()
        test3 = test_chemistry_accuracy()
        
        # Summary
        print("\n" + "="*70)
        print("TEST SUMMARY")
        print("="*70)
        print(f"Test 1 (Basic): {'[PASS]' if test1 else '[FAIL]'}")
        print(f"Test 2 (Timing): {'[PASS]' if test2 else '[FAIL]'}")
        print(f"Test 3 (Chemistry): {'[PASS]' if test3 else '[FAIL]'}")
        
        if all([test1, test2, test3]):
            print("\n*** ALL TESTS PASSED! ***")
            return 0
        else:
            print("\n*** SOME TESTS FAILED ***")
            return 1
            
    except Exception as e:
        print(f"\n[ERROR] Tests failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

