#!/usr/bin/env python3
"""
Test script for thermodynamic validation
Runs a long simulation and validates thermodynamic laws
"""

import sys
import os
import time
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gc  # Garbage collection
import psutil  # Memory monitoring

# Add current directory to path
sys.path.append('.')

import taichi as ti

# Initialize Taichi with memory limits
try:
    ti.init(arch=ti.cpu, cpu_max_num_threads=2, device_memory_fraction=0.3)
except Exception as e:
    print(f"Taichi init failed: {e}")
    try:
        ti.init(arch=ti.cpu, cpu_max_num_threads=1, device_memory_fraction=0.1)
    except Exception as e2:
        print(f"Fallback Taichi init failed: {e2}")
        exit(1)

from backend.sim.config import SimulationConfig
from backend.sim.core.stepper import SimulationStepper

def run_thermodynamic_validation_test():
    """Run a long simulation to test thermodynamic validation"""
    
    print("Starting Thermodynamic Validation Test")
    print("=" * 50)
    
    # Create configuration - ULTRA OPTIMIZED FOR MEMORY
    config = SimulationConfig()
    config.max_particles = 200  # Much smaller system
    config.grid_width = 32  # Much smaller grid
    config.grid_height = 32  # Much smaller grid
    config.dt = 0.01
    config.validate_every_n_steps = 200  # Less frequent validation
    config.energy_tolerance = 1e-3  # 0.1% energy tolerance
    config.momentum_tolerance = 1e-4  # 0.01% momentum tolerance
    
    print(f"Configuration:")
    print(f"  Max particles: {config.max_particles}")
    print(f"  Grid size: {config.grid_width}x{config.grid_height}")
    print(f"  Timestep: {config.dt}")
    print(f"  Validation interval: {config.validate_every_n_steps}")
    print(f"  Energy tolerance: {config.energy_tolerance}")
    print(f"  Momentum tolerance: {config.momentum_tolerance}")
    
    # Create simulation
    print("\nCreating simulation...")
    simulation = SimulationStepper(config)
    simulation.start()
    
    # Run simulation
    print("\nRunning simulation...")
    start_time = time.time()
    target_steps = 10000  # Reduced to 10k steps for memory optimization
    
    try:
        for step in range(target_steps):
            simulation.step()
            
            # Progress update every 5000 steps
            if step % 5000 == 0 and step > 0:
                elapsed = time.time() - start_time
                steps_per_sec = step / elapsed
                eta = (target_steps - step) / steps_per_sec
                print(f"  Step {step}/{target_steps} ({step/target_steps*100:.1f}%) - "
                      f"{steps_per_sec:.1f} steps/sec - ETA: {eta:.1f}s")
                
                # Print current validation summary
                summary = simulation.get_thermodynamic_validation_summary()
                if 'overall_success_rate' in summary:
                    print(f"    Validation success rate: {summary['overall_success_rate']*100:.1f}%")
                
                # Monitor validation performance
                validation_log = simulation.get_validation_log()
                if validation_log:
                    last_validation = validation_log[-1]
                    validation_time = last_validation.get('validation_time', 0)
                    if validation_time > 0:
                        print(f"    Last validation took: {validation_time:.3f}s")
                
                # Print current energy drift
                try:
                    energy_drift = simulation._get_energy_drift()
                    print(f"    Energy drift: {energy_drift:.4f}%")
                except:
                    pass
                
                # Force garbage collection to free memory
                gc.collect()
                
                # Monitor memory usage
                memory_percent = psutil.virtual_memory().percent
                print(f"    Memory usage: {memory_percent:.1f}%")
                
                if memory_percent > 90:
                    print("    WARNING: High memory usage! Consider stopping simulation.")
    
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
    print("=" * 40)
    
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
    
    # Save validation log
    validation_log = simulation.get_validation_log()
    if validation_log:
        log_file = "thermodynamic_validation_log.json"
        with open(log_file, 'w') as f:
            json.dump(validation_log, f, indent=2)
        print(f"\nValidation log saved to: {log_file}")
    
    # Generate Figure 1 (Energy Conservation)
    print("\nGenerating Figure 1 (Energy Conservation)...")
    generate_energy_conservation_plot(validation_log)
    
    return True

def generate_energy_conservation_plot(validation_log):
    """Generate Figure 1: Energy conservation plot"""
    
    if not validation_log:
        print("  No validation data available for plotting")
        return
    
    # Extract energy conservation data
    steps = []
    energy_errors = []
    energy_before = []
    energy_after = []
    energy_expected = []
    
    for entry in validation_log:
        if 'energy' in entry['results']:
            steps.append(entry['step'])
            energy_data = entry['results']['energy']['details']
            energy_errors.append(entry['results']['energy']['error'])
            energy_before.append(energy_data['E_before'])
            energy_after.append(energy_data['E_after'])
            energy_expected.append(energy_data['E_expected'])
    
    if not steps:
        print("  No energy conservation data found")
        return
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Panel A: Energy conservation
    axes[0].plot(steps, energy_after, label='Actual', alpha=0.7, linewidth=1)
    axes[0].plot(steps, energy_expected, '--', label='Expected', alpha=0.7, linewidth=1)
    
    # Add tolerance band
    if energy_expected:
        tolerance = 0.001  # 0.1%
        upper_bound = np.array(energy_expected) * (1 + tolerance)
        lower_bound = np.array(energy_expected) * (1 - tolerance)
        axes[0].fill_between(steps, lower_bound, upper_bound, 
                           alpha=0.2, label='±0.1% tolerance')
    
    axes[0].set_xlabel('Simulation Step')
    axes[0].set_ylabel('Total Energy')
    axes[0].legend()
    axes[0].set_title('A) Energy Conservation')
    axes[0].grid(True, alpha=0.3)
    
    # Panel B: Relative error
    axes[1].semilogy(steps, energy_errors, linewidth=1)
    axes[1].axhline(1e-3, color='r', linestyle='--', 
                    label='Tolerance threshold (0.1%)')
    axes[1].set_xlabel('Simulation Step')
    axes[1].set_ylabel('Relative Error')
    axes[1].legend()
    axes[1].set_title('B) Relative Error')
    axes[1].grid(True, alpha=0.3)
    
    # Panel C: Cumulative drift
    cumulative_drift = np.cumsum(energy_errors)
    axes[2].plot(steps, cumulative_drift, linewidth=1)
    axes[2].set_xlabel('Simulation Step')
    axes[2].set_ylabel('Cumulative Energy Drift')
    axes[2].set_title('C) Cumulative Drift')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_path = "fig1_energy_conservation.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Figure saved to: {output_path}")
    
    # Also save as PDF for publication
    pdf_path = "fig1_energy_conservation.pdf"
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"  PDF saved to: {pdf_path}")
    
    plt.close()

def main():
    """Main function"""
    print("Live 2.0 - Thermodynamic Validation Test")
    print("Week 1: Fundamental Physics Validation")
    print()
    
    success = run_thermodynamic_validation_test()
    
    if success:
        print("\nThermodynamic validation test completed successfully!")
        print("\nNext steps:")
        print("  1. Review Figure 1 (energy conservation)")
        print("  2. Check validation log for any failures")
        print("  3. Adjust tolerances if needed")
        print("  4. Proceed to Week 2 (Literature Parameters)")
    else:
        print("\nThermodynamic validation test failed!")
        print("  Check logs for errors and fix before proceeding")

if __name__ == "__main__":
    main()
