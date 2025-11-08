#!/usr/bin/env python3
"""
Benchmark GPU vs CPU Performance for Live 2.0 Simulation
Tests both simulation steps and visualization performance
"""

import sys
import time
import argparse
import numpy as np
from pathlib import Path
import multiprocessing

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import taichi as ti


def test_taichi_backend(backend_name: str, max_threads: int = None):
    """Test a specific Taichi backend"""
    print(f"\n{'='*70}")
    print(f"TESTING: {backend_name}")
    print(f"{'='*70}")
    
    # Reset Taichi
    ti.reset()
    
    # Initialize backend
    try:
        if backend_name == "cuda":
            ti.init(arch=ti.cuda, device_memory_GB=4.0)
        elif backend_name == "vulkan":
            ti.init(arch=ti.vulkan, device_memory_GB=4.0)
        elif backend_name == "cpu":
            if max_threads is None:
                max_threads = multiprocessing.cpu_count()
            ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)
            print(f"  CPU threads: {max_threads}")
        else:
            raise ValueError(f"Unknown backend: {backend_name}")
            
        print(f"  ✅ Backend initialized successfully")
    except Exception as e:
        print(f"  ❌ Backend initialization failed: {e}")
        return None
    
    # Import after Taichi is initialized
    from backend.sim.config import SimulationConfig
    from backend.sim.core.stepper import SimulationStepper
    
    # Create small simulation for testing
    config = SimulationConfig(
        n_particles=500,  # Small for quick test
        max_steps=1000,
        dt=0.001,
        box_size=100.0,
        temperature=298.0,
        mode='open_chemistry',
        seed=42
    )
    
    # Disable heavy operations for pure performance test
    config.detect_novel_substances = False
    config.enable_diagnostics = False
    config.metrics_update_interval = 999999
    config.novelty_check_interval = 999999
    config.mutation_interval = 999999
    
    print(f"  Creating simulation with {config.n_particles} particles...")
    
    try:
        stepper = SimulationStepper(config)
        stepper.start()
        print(f"  ✅ Simulation created successfully")
    except Exception as e:
        print(f"  ❌ Simulation creation failed: {e}")
        return None
    
    # Warm-up (5 steps)
    print(f"  Warming up...")
    for _ in range(5):
        stepper.step()
    
    # Benchmark simulation steps
    print(f"  Benchmarking simulation steps...")
    num_steps = 100
    start_time = time.time()
    
    for i in range(num_steps):
        stepper.step()
    
    elapsed = time.time() - start_time
    steps_per_sec = num_steps / elapsed
    
    print(f"  ✅ Simulation benchmark complete")
    print(f"     Steps: {num_steps}")
    print(f"     Time: {elapsed:.2f}s")
    print(f"     Speed: {steps_per_sec:.1f} steps/sec")
    print(f"     Per step: {elapsed/num_steps*1000:.1f}ms")
    
    # Benchmark visualization
    print(f"  Benchmarking visualization...")
    num_viz = 50
    viz_times = []
    
    for i in range(num_viz):
        start = time.time()
        viz_data = stepper.get_visualization_data()
        viz_times.append(time.time() - start)
        
        # Do a step between visualizations
        stepper.step()
    
    avg_viz_time = np.mean(viz_times)
    max_viz_time = np.max(viz_times)
    min_viz_time = np.min(viz_times)
    
    print(f"  ✅ Visualization benchmark complete")
    print(f"     Calls: {num_viz}")
    print(f"     Avg time: {avg_viz_time*1000:.1f}ms")
    print(f"     Min time: {min_viz_time*1000:.1f}ms")
    print(f"     Max time: {max_viz_time*1000:.1f}ms")
    
    # Memory usage
    try:
        import psutil
        process = psutil.Process()
        memory_mb = process.memory_info().rss / 1024 / 1024
        print(f"  Memory usage: {memory_mb:.1f} MB")
    except:
        pass
    
    # Cleanup
    stepper.stop()
    
    return {
        'backend': backend_name,
        'steps_per_sec': steps_per_sec,
        'step_time_ms': elapsed/num_steps*1000,
        'viz_avg_ms': avg_viz_time*1000,
        'viz_max_ms': max_viz_time*1000,
        'viz_min_ms': min_viz_time*1000,
        'success': True
    }


def main():
    parser = argparse.ArgumentParser(description='Benchmark GPU vs CPU performance')
    parser.add_argument('--backends', nargs='+', 
                       default=['cuda', 'cpu'],
                       choices=['cuda', 'vulkan', 'cpu'],
                       help='Backends to test (default: cuda cpu)')
    parser.add_argument('--cpu-threads', type=int, nargs='+',
                       default=None,
                       help='CPU thread counts to test (default: max, max/2, max/4)')
    args = parser.parse_args()
    
    print("="*70)
    print("LIVE 2.0 - GPU vs CPU PERFORMANCE BENCHMARK")
    print("="*70)
    
    results = []
    
    # Test requested backends
    for backend in args.backends:
        if backend == 'cpu':
            # Test different CPU thread counts
            max_cores = multiprocessing.cpu_count()
            
            if args.cpu_threads:
                thread_counts = args.cpu_threads
            else:
                # Default: test max, half, quarter
                thread_counts = [max_cores, max(1, max_cores//2), max(1, max_cores//4)]
            
            for threads in thread_counts:
                result = test_taichi_backend(f"cpu-{threads}threads", max_threads=threads)
                if result:
                    result['backend'] = f"CPU ({threads} threads)"
                    results.append(result)
        else:
            result = test_taichi_backend(backend)
            if result:
                results.append(result)
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    if not results:
        print("❌ No successful tests")
        return
    
    # Sort by steps per second (descending)
    results.sort(key=lambda x: x['steps_per_sec'], reverse=True)
    
    print("\n📊 Simulation Performance (steps/sec):")
    print("-"*70)
    for i, r in enumerate(results):
        speedup = r['steps_per_sec'] / results[-1]['steps_per_sec']
        print(f"{i+1}. {r['backend']:20s} - {r['steps_per_sec']:6.1f} steps/sec "
              f"({r['step_time_ms']:5.1f}ms/step) [{speedup:.1f}x faster]")
    
    print("\n📊 Visualization Performance (lower is better):")
    print("-"*70)
    results_viz = sorted(results, key=lambda x: x['viz_avg_ms'])
    for i, r in enumerate(results_viz):
        speedup = results_viz[-1]['viz_avg_ms'] / r['viz_avg_ms']
        print(f"{i+1}. {r['backend']:20s} - {r['viz_avg_ms']:6.1f}ms avg "
              f"(min: {r['viz_min_ms']:5.1f}ms, max: {r['viz_max_ms']:6.1f}ms) [{speedup:.1f}x faster]")
    
    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    
    best_sim = results[0]
    best_viz = results_viz[0]
    
    print(f"\n🚀 Best for simulation steps: {best_sim['backend']}")
    print(f"   Speed: {best_sim['steps_per_sec']:.1f} steps/sec")
    
    print(f"\n🎨 Best for visualization: {best_viz['backend']}")
    print(f"   Speed: {best_viz['viz_avg_ms']:.1f}ms per frame")
    
    # Overall recommendation
    if best_sim['backend'] == best_viz['backend']:
        print(f"\n✅ WINNER: {best_sim['backend']}")
        print(f"   Best for both simulation and visualization!")
    else:
        # Calculate weighted score (70% simulation, 30% visualization)
        for r in results:
            sim_score = r['steps_per_sec'] / results[0]['steps_per_sec']
            viz_score = results_viz[0]['viz_avg_ms'] / r['viz_avg_ms']
            r['combined_score'] = 0.7 * sim_score + 0.3 * viz_score
        
        best_overall = max(results, key=lambda x: x['combined_score'])
        print(f"\n✅ WINNER (weighted): {best_overall['backend']}")
        print(f"   Combined score: {best_overall['combined_score']:.2f}")
        print(f"   (70% simulation performance, 30% visualization performance)")
    
    print("\n" + "="*70)
    
    # Save results to file
    import json
    output_file = Path(__file__).parent / "benchmark_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n📝 Results saved to: {output_file}")


if __name__ == "__main__":
    main()

