#!/usr/bin/env python3
"""
Performance test for Live 2.0
Tests 10,000 particles @ target 60 FPS
"""

import requests
import json
import time
import sys
from typing import Dict, Any

API_BASE = "http://localhost:8001"

def create_simulation(particle_count: int = 10000) -> str:
    """Create a simulation with specified particle count"""
    config = {
        "config": {
            "max_particles": particle_count,
            "grid_width": 256,
            "grid_height": 256,
            "dt": 0.01,
            "mode": "open_chemistry"
        },
        "mode": "open_chemistry"
    }
    
    print(f"Creating simulation with {particle_count} max particles...")
    response = requests.post(f"{API_BASE}/simulation/create", json=config)
    response.raise_for_status()
    data = response.json()
    
    simulation_id = data.get("simulation_id")
    print(f"Simulation created: {simulation_id}")
    return simulation_id

def start_simulation(simulation_id: str):
    """Start the simulation"""
    print(f"Starting simulation {simulation_id}...")
    response = requests.post(f"{API_BASE}/simulation/{simulation_id}/start")
    response.raise_for_status()
    print("Simulation started")

def get_status(simulation_id: str) -> Dict[str, Any]:
    """Get simulation status"""
    response = requests.get(f"{API_BASE}/simulation/{simulation_id}/status")
    response.raise_for_status()
    return response.json()

def measure_performance(simulation_id: str, duration_seconds: int = 10):
    """Measure simulation performance over a period"""
    print(f"\nMeasuring performance for {duration_seconds} seconds...")
    
    initial_status = get_status(simulation_id)
    initial_step = initial_status['step_count']
    initial_time = time.time()
    
    time.sleep(duration_seconds)
    
    final_status = get_status(simulation_id)
    final_step = final_status['step_count']
    final_time = time.time()
    
    steps_per_second = (final_step - initial_step) / (final_time - initial_time)
    
    print(f"\n{'='*60}")
    print(f"PERFORMANCE TEST RESULTS")
    print(f"{'='*60}")
    print(f"Particles:       {final_status['particle_count']:,} / {final_status.get('max_particles', 10000):,}")
    print(f"Duration:        {final_time - initial_time:.2f} seconds")
    print(f"Steps completed: {final_step - initial_step:,}")
    print(f"Steps/second:    {steps_per_second:.2f}")
    print(f"Target FPS:      60")
    print(f"Actual FPS:      {steps_per_second:.2f}")
    print(f"Performance:     {(steps_per_second / 60.0) * 100:.1f}% of target")
    print(f"")
    print(f"Sim time:        {final_status['current_time']:.2f}")
    print(f"Novelty rate:    {final_status['novelty_rate']:.4f}")
    print(f"Health score:    {final_status['health_score']:.4f}")
    
    # Check energy conservation if available
    response = requests.get(f"{API_BASE}/simulation/{simulation_id}/status")
    if response.status_code == 200:
        full_status = response.json()
        if 'energy_drift' in full_status:
            print(f"Energy drift:    {full_status['energy_drift']:.2f}%")
        if 'adaptive_dt' in full_status:
            print(f"Adaptive dt:     {full_status['adaptive_dt']:.6f}")
    
    print(f"{'='*60}\n")
    
    return steps_per_second

def stop_simulation(simulation_id: str):
    """Stop the simulation"""
    print(f"Stopping simulation {simulation_id}...")
    response = requests.post(f"{API_BASE}/simulation/{simulation_id}/stop")
    response.raise_for_status()
    print("Simulation stopped")

def main():
    particle_count = 10000
    test_duration = 10
    
    if len(sys.argv) > 1:
        try:
            particle_count = int(sys.argv[1])
        except ValueError:
            print(f"Invalid particle count: {sys.argv[1]}")
            sys.exit(1)
    
    if len(sys.argv) > 2:
        try:
            test_duration = int(sys.argv[2])
        except ValueError:
            print(f"Invalid test duration: {sys.argv[2]}")
            sys.exit(1)
    
    print("="*60)
    print("Live 2.0 Performance Test")
    print("="*60)
    print(f"Target: {particle_count:,} particles @ 60 FPS")
    print(f"Test duration: {test_duration} seconds")
    print("")
    
    try:
        # Create simulation
        sim_id = create_simulation(particle_count)
        
        # Start simulation
        start_simulation(sim_id)
        
        # Wait for initialization
        print("Waiting for simulation to initialize...")
        time.sleep(2)
        
        # Measure performance
        fps = measure_performance(sim_id, test_duration)
        
        # Stop simulation
        stop_simulation(sim_id)
        
        # Return exit code based on performance
        if fps >= 60:
            print("PASS: Performance target met!")
            sys.exit(0)
        elif fps >= 30:
            print("PARTIAL: Performance acceptable but below target")
            sys.exit(0)
        else:
            print("FAIL: Performance below acceptable threshold")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n\nTest interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
