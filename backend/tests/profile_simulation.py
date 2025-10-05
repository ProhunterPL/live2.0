#!/usr/bin/env python3
"""
Profile simulation performance to identify bottlenecks
"""

import cProfile
import pstats
import io
import sys
sys.path.insert(0, 'backend')

import taichi as ti

# Initialize Taichi FIRST - force CUDA for performance
try:
    ti.init(arch=ti.cuda)
    print("Using CUDA backend")
except Exception as e:
    print(f"CUDA failed: {e}")
    ti.init(arch=ti.cpu)
    print("Using CPU backend (fallback)")

# Import AFTER Taichi initialization
from sim.config import SimulationConfig
from sim.core.stepper import SimulationStepper

def profile_simulation_steps(num_particles=1000, num_steps=100):
    """Profile a simulation with specified particles and steps"""
    
    print(f"Profiling simulation: {num_particles} particles, {num_steps} steps")
    
    # Create configuration
    config = SimulationConfig(
        max_particles=num_particles,
        grid_width=128,
        grid_height=128,
        dt=0.01,
        mode="open_chemistry"
    )
    
    # Create simulation
    sim = SimulationStepper(config)
    sim.initialize_simulation()
    sim.start()
    
    # Profile the steps
    profiler = cProfile.Profile()
    profiler.enable()
    
    for i in range(num_steps):
        sim.step()
    
    profiler.disable()
    
    # Print profiling results
    s = io.StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
    ps.print_stats(30)  # Top 30 functions
    
    print("\n" + "="*80)
    print("PROFILING RESULTS (Top 30 by cumulative time)")
    print("="*80)
    print(s.getvalue())
    
    # Calculate performance metrics
    total_time = sum(ps.stats[k][3] for k in ps.stats)
    steps_per_second = num_steps / total_time if total_time > 0 else 0
    
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY")
    print("="*80)
    print(f"Total time:       {total_time:.3f} seconds")
    print(f"Steps completed:  {num_steps}")
    print(f"Steps/second:     {steps_per_second:.2f}")
    print(f"Target FPS:       60")
    print(f"Performance:      {(steps_per_second / 60.0) * 100:.1f}% of target")
    print("="*80)
    
    sim.stop()

if __name__ == "__main__":
    num_particles = 1000
    num_steps = 100
    
    if len(sys.argv) > 1:
        try:
            num_particles = int(sys.argv[1])
        except ValueError:
            print(f"Invalid particle count: {sys.argv[1]}")
            sys.exit(1)
    
    if len(sys.argv) > 2:
        try:
            num_steps = int(sys.argv[2])
        except ValueError:
            print(f"Invalid step count: {sys.argv[2]}")
            sys.exit(1)
    
    profile_simulation_steps(num_particles, num_steps)
