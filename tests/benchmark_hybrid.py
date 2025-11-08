#!/usr/bin/env python3
"""
Benchmark: Pure GPU vs Pure CPU vs Hybrid GPU+CPU

Tests which approach is fastest for LIVE 2.0 simulation
"""

import sys
import time
import argparse
from pathlib import Path
import numpy as np

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import taichi as ti


def benchmark_stepper(stepper_class, config, name: str, num_steps: int = 200):
    """Benchmark a specific stepper implementation"""
    print(f"\n{'='*70}")
    print(f"TESTING: {name}")
    print(f"{'='*70}")
    
    try:
        # Create stepper
        print(f"  Creating stepper...")
        stepper = stepper_class(config)
        stepper.start()
        
        # Warm-up
        print(f"  Warming up (10 steps)...")
        for _ in range(10):
            stepper.step()
        
        # Benchmark simulation
        print(f"  Benchmarking simulation ({num_steps} steps)...")
        start_time = time.time()
        
        step_times = []
        for i in range(num_steps):
            step_start = time.time()
            stepper.step()
            step_times.append(time.time() - step_start)
        
        sim_elapsed = time.time() - start_time
        steps_per_sec = num_steps / sim_elapsed
        avg_step_ms = np.mean(step_times) * 1000
        max_step_ms = np.max(step_times) * 1000
        
        print(f"  ✅ Simulation benchmark complete")
        print(f"     Steps: {num_steps}")
        print(f"     Time: {sim_elapsed:.2f}s")
        print(f"     Speed: {steps_per_sec:.1f} steps/sec")
        print(f"     Avg step: {avg_step_ms:.1f}ms")
        print(f"     Max step: {max_step_ms:.1f}ms")
        
        # Benchmark visualization
        print(f"  Benchmarking visualization (50 calls)...")
        viz_times = []
        
        for i in range(50):
            viz_start = time.time()
            viz_data = stepper.get_visualization_data()
            viz_times.append(time.time() - viz_start)
            
            # Do a step
            stepper.step()
        
        avg_viz_ms = np.mean(viz_times) * 1000
        max_viz_ms = np.max(viz_times) * 1000
        min_viz_ms = np.min(viz_times) * 1000
        
        print(f"  ✅ Visualization benchmark complete")
        print(f"     Avg: {avg_viz_ms:.1f}ms")
        print(f"     Min: {min_viz_ms:.1f}ms")
        print(f"     Max: {max_viz_ms:.1f}ms")
        
        # Get state
        state = stepper.get_simulation_state()
        
        # Stop
        stepper.stop()
        
        return {
            'name': name,
            'success': True,
            'steps_per_sec': steps_per_sec,
            'avg_step_ms': avg_step_ms,
            'max_step_ms': max_step_ms,
            'viz_avg_ms': avg_viz_ms,
            'viz_max_ms': max_viz_ms,
            'viz_min_ms': min_viz_ms,
            'particle_count': state.get('particle_count', 0),
            'bond_count': state.get('bond_count', 0),
            'cluster_count': state.get('cluster_count', 0)
        }
        
    except Exception as e:
        print(f"  ❌ Benchmark failed: {e}")
        import traceback
        traceback.print_exc()
        return {
            'name': name,
            'success': False,
            'error': str(e)
        }


def main():
    parser = argparse.ArgumentParser(description='Benchmark Pure vs Hybrid steppers')
    parser.add_argument('--steps', type=int, default=200,
                       help='Number of steps to benchmark (default: 200)')
    parser.add_argument('--particles', type=int, default=500,
                       help='Number of particles (default: 500)')
    parser.add_argument('--modes', nargs='+',
                       default=['pure_gpu', 'pure_cpu', 'hybrid'],
                       choices=['pure_gpu', 'pure_cpu', 'hybrid'],
                       help='Modes to test (default: all)')
    args = parser.parse_args()
    
    print("="*70)
    print("LIVE 2.0 - PURE vs HYBRID BENCHMARK")
    print("="*70)
    print(f"Steps: {args.steps}")
    print(f"Particles: {args.particles}")
    print(f"Modes: {', '.join(args.modes)}")
    
    # Import after printing (to avoid Taichi init messages)
    from backend.sim.config import SimulationConfig
    from backend.sim.core.stepper import SimulationStepper
    from backend.sim.core.hybrid_stepper import HybridSimulationStepper
    
    results = []
    
    # Test Pure GPU
    if 'pure_gpu' in args.modes:
        ti.reset()
        try:
            ti.init(arch=ti.cuda, device_memory_GB=4.0)
            
            config = SimulationConfig(
                n_particles=args.particles,
                max_steps=10000,
                dt=0.001,
                box_size=100.0,
                temperature=298.0,
                mode='open_chemistry',
                seed=42
            )
            
            # Disable heavy operations for fair comparison
            config.detect_novel_substances = False
            config.enable_diagnostics = False
            config.metrics_update_interval = 999999
            config.novelty_check_interval = 999999
            config.mutation_interval = 999999
            
            result = benchmark_stepper(SimulationStepper, config, "Pure GPU (CUDA)", args.steps)
            results.append(result)
            
        except Exception as e:
            print(f"\n❌ Pure GPU not available: {e}")
    
    # Test Pure CPU
    if 'pure_cpu' in args.modes:
        ti.reset()
        import multiprocessing
        max_threads = multiprocessing.cpu_count()
        ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)
        
        config = SimulationConfig(
            n_particles=args.particles,
            max_steps=10000,
            dt=0.001,
            box_size=100.0,
            temperature=298.0,
            mode='open_chemistry',
            seed=42
        )
        
        config.detect_novel_substances = False
        config.enable_diagnostics = False
        config.metrics_update_interval = 999999
        config.novelty_check_interval = 999999
        config.mutation_interval = 999999
        
        result = benchmark_stepper(SimulationStepper, config, f"Pure CPU ({max_threads} threads)", args.steps)
        results.append(result)
    
    # Test Hybrid
    if 'hybrid' in args.modes:
        ti.reset()
        try:
            ti.init(arch=ti.cuda, device_memory_GB=4.0)
            
            config = SimulationConfig(
                n_particles=args.particles,
                max_steps=10000,
                dt=0.001,
                box_size=100.0,
                temperature=298.0,
                mode='open_chemistry',
                seed=42
            )
            
            # Hybrid mode: GPU physics, CPU chemistry
            config.detect_novel_substances = False
            config.enable_diagnostics = False
            config.metrics_update_interval = 999999
            config.novelty_check_interval = 999999
            config.mutation_interval = 999999
            config.chemistry_snapshot_interval = 100  # Send to CPU every 100 steps
            
            result = benchmark_stepper(HybridSimulationStepper, config, "Hybrid (GPU physics + CPU chemistry)", args.steps)
            results.append(result)
            
        except Exception as e:
            print(f"\n❌ Hybrid mode not available: {e}")
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    if not results:
        print("❌ No successful tests")
        return
    
    # Filter successful results
    successful = [r for r in results if r.get('success', False)]
    
    if not successful:
        print("❌ No successful tests")
        return
    
    # Sort by steps per second
    successful.sort(key=lambda x: x['steps_per_sec'], reverse=True)
    
    print("\n📊 Simulation Performance (steps/sec):")
    print("-"*70)
    for i, r in enumerate(successful):
        speedup = r['steps_per_sec'] / successful[-1]['steps_per_sec']
        print(f"{i+1}. {r['name']:40s} - {r['steps_per_sec']:6.1f} steps/sec "
              f"({r['avg_step_ms']:5.1f}ms/step) [{speedup:.1f}x]")
    
    print("\n📊 Visualization Performance (lower is better):")
    print("-"*70)
    viz_sorted = sorted(successful, key=lambda x: x['viz_avg_ms'])
    for i, r in enumerate(viz_sorted):
        speedup = viz_sorted[-1]['viz_avg_ms'] / r['viz_avg_ms']
        print(f"{i+1}. {r['name']:40s} - {r['viz_avg_ms']:6.1f}ms avg "
              f"[{speedup:.1f}x faster]")
    
    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    
    best_sim = successful[0]
    best_viz = viz_sorted[0]
    
    print(f"\n🚀 Best for simulation: {best_sim['name']}")
    print(f"   Speed: {best_sim['steps_per_sec']:.1f} steps/sec")
    
    print(f"\n🎨 Best for visualization: {best_viz['name']}")
    print(f"   Speed: {best_viz['viz_avg_ms']:.1f}ms per frame")
    
    # Overall winner
    if best_sim['name'] == best_viz['name']:
        print(f"\n✅ WINNER: {best_sim['name']}")
    else:
        # Calculate weighted score
        for r in successful:
            sim_score = r['steps_per_sec'] / successful[0]['steps_per_sec']
            viz_score = viz_sorted[0]['viz_avg_ms'] / r['viz_avg_ms']
            r['combined_score'] = 0.7 * sim_score + 0.3 * viz_score
        
        best_overall = max(successful, key=lambda x: x['combined_score'])
        print(f"\n✅ OVERALL WINNER: {best_overall['name']}")
        print(f"   Combined score: {best_overall['combined_score']:.2f}")
    
    # Save results
    import json
    output_file = Path(__file__).parent / "hybrid_benchmark_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n📝 Results saved to: {output_file}")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    main()

