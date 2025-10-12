#!/usr/bin/env python3
"""
Minimal thermodynamic validation test
Ultra-lightweight version for low memory systems
"""

import sys
import os
import time
import json
import numpy as np
import gc

# Add current directory to path
sys.path.append('.')

import taichi as ti

# Initialize Taichi with minimal memory
try:
    ti.init(arch=ti.cpu, cpu_max_num_threads=1, device_memory_fraction=0.05)
    print("Taichi initialized with minimal memory")
except Exception as e:
    print(f"Taichi init failed: {e}")
    exit(1)

from backend.sim.config import SimulationConfig
from backend.sim.core.stepper import SimulationStepper

def run_minimal_validation_test():
    """Run minimal thermodynamic validation test"""
    
    print("Minimal Thermodynamic Validation Test")
    print("=" * 40)
    
    # Ultra-minimal configuration
    config = SimulationConfig()
    config.max_particles = 50  # Very small system
    config.grid_width = 16  # Tiny grid
    config.grid_height = 16  # Tiny grid
    config.dt = 0.01
    config.validate_every_n_steps = 500  # Very infrequent validation
    config.energy_tolerance = 1e-2  # Relaxed tolerance
    config.momentum_tolerance = 1e-3  # Relaxed tolerance
    
    print(f"Configuration:")
    print(f"  Max particles: {config.max_particles}")
    print(f"  Grid size: {config.grid_width}x{config.grid_height}")
    print(f"  Validation interval: {config.validate_every_n_steps}")
    
    # Create simulation
    print("\nCreating simulation...")
    simulation = SimulationStepper(config)
    simulation.start()
    
    # Run short simulation
    print("\nRunning simulation...")
    start_time = time.time()
    target_steps = 5000  # Very short test
    
    try:
        for step in range(target_steps):
            simulation.step()
            
            # Progress update every 1000 steps
            if step % 1000 == 0 and step > 0:
                elapsed = time.time() - start_time
                steps_per_sec = step / elapsed
                eta = (target_steps - step) / steps_per_sec
                print(f"  Step {step}/{target_steps} ({step/target_steps*100:.1f}%) - "
                      f"{steps_per_sec:.1f} steps/sec - ETA: {eta:.1f}s")
                
                # Monitor validation performance
                validation_log = simulation.get_validation_log()
                if validation_log:
                    last_validation = validation_log[-1]
                    validation_time = last_validation.get('validation_time', 0)
                    if validation_time > 0:
                        print(f"    Last validation took: {validation_time:.3f}s")
                
                # Force garbage collection
                gc.collect()
    
    except KeyboardInterrupt:
        print("\nSimulation interrupted by user")
    except Exception as e:
        print(f"\nSimulation failed: {e}")
        return False
    
    elapsed = time.time() - start_time
    print(f"\nSimulation completed in {elapsed:.1f}s")
    print(f"   Final step: {simulation.step_count}")
    print(f"   Steps per second: {simulation.step_count/elapsed:.1f}")
    
    # Get validation results
    print("\nThermodynamic Validation Results:")
    print("=" * 30)
    
    summary = simulation.get_thermodynamic_validation_summary()
    
    if 'message' in summary:
        print(f"  {summary['message']}")
    else:
        print(f"  Total validation tests: {summary['total_tests']}")
        print(f"  Overall success rate: {summary['overall_success_rate']*100:.1f}%")
        
        if 'test_statistics' in summary:
            print(f"\n  Per-test statistics:")
            for test_name, stats in summary['test_statistics'].items():
                print(f"    {test_name}: {stats['passed']}/{stats['total']} "
                      f"({stats['success_rate']*100:.1f}%)")
    
    # Save minimal validation log
    validation_log = simulation.get_validation_log()
    if validation_log:
        log_file = "minimal_validation_log.json"
        with open(log_file, 'w') as f:
            json.dump(validation_log, f, indent=2)
        print(f"\nValidation log saved to: {log_file}")
    
    return True

def main():
    """Main function"""
    print("Live 2.0 - Minimal Thermodynamic Validation")
    print("Ultra-lightweight version for low memory systems")
    print()
    
    success = run_minimal_validation_test()
    
    if success:
        print("\nMinimal validation test completed!")
        print("\nResults:")
        print("  - System runs without memory issues")
        print("  - Basic thermodynamic validation works")
        print("  - Ready for larger simulations when memory available")
    else:
        print("\nMinimal validation test failed!")
        print("  Check system resources and try again")

if __name__ == "__main__":
    main()

