#!/usr/bin/env python3
"""
Quick Phase 2 Runner - Bez kosztownej detekcji klastrów
"""

import sys
import time
import argparse
from pathlib import Path
from backend.sim.phase2_config import load_phase2_config_from_yaml
import taichi as ti
import logging

logging.basicConfig(level=logging.INFO)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--steps', type=int, default=500000)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    
    # Initialize Taichi - SPRÓBUJ GPU
    try:
        ti.init(arch=ti.cuda)
        print("✅ Using GPU acceleration")
    except:
        import multiprocessing
        max_threads = multiprocessing.cpu_count()
        ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)
        print(f"⚠️ Using CPU with {max_threads} threads")
    
    # Load config
    phase2_config = load_phase2_config_from_yaml(args.config)
    
    # Import stepper
    from backend.sim.core.stepper import SimulationStepper
    from backend.sim.run_simulation import create_simulation
    
    # Create simulation
    sim_config = create_simulation()
    stepper = SimulationStepper(sim_config)
    
    # WYŁĄCZ detekcję klastrów (wewnątrz stepper)
    stepper.detection_enabled = False  # Jeśli możliwe
    # Lub zwiększ interwał:
    # stepper.detection_interval = 1000000
    
    # Run
    stepper.start()
    start_time = time.time()
    
    for step in range(args.steps):
        stepper.step()
        
        if step % 10000 == 0:
            elapsed = time.time() - start_time
            print(f"Step {step:,}/{args.steps:,} | Speed: {step/elapsed:.0f} steps/s | ETA: {(args.steps-step)/(step/elapsed)/60:.1f} min")
    
    print(f"\n✅ Completed {args.steps} steps in {time.time()-start_time:.0f} seconds")

if __name__ == "__main__":
    main()

