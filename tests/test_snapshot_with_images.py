#!/usr/bin/env python3
"""
Test snapshot generation with image creation
"""

import sys
import os
sys.path.insert(0, 'backend')

import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend

import requests
import json
import time

API_BASE = "http://localhost:8001"

def test_snapshot_with_images():
    """Test snapshot creation with image generation"""
    print("Testing snapshot with image generation...")
    
    try:
        # Check if backend is running
        response = requests.get(f"{API_BASE}/")
        if response.status_code != 200:
            print("Backend not running! Start it with: cd backend && python -m api.server")
            return False
            
        # Create simulation
        config = {
            "config": {
                "grid_height": 64,
                "grid_width": 64,
                "max_particles": 100,
                "dt": 0.035,
                "mode": "open_chemistry"
            },
            "mode": "open_chemistry"
        }
        
        print("Creating simulation...")
        response = requests.post(f"{API_BASE}/simulation/create", json=config)
        response.raise_for_status()
        data = response.json()
        simulation_id = data.get("simulation_id")
        print(f"Created simulation: {simulation_id}")
        
        # Start simulation
        print("Starting simulation...")
        response = requests.post(f"{API_BASE}/simulation/{simulation_id}/start")
        response.raise_for_status()
        
        # Wait for some activity
        print("Running for 10 seconds...")
        time.sleep(10)
        
        # Check status
        response = requests.get(f"{API_BASE}/simulation/{simulation_id}/status")
        response.raise_for_status()
        status = response.json()
        print(f"Simulation status: {status}")
        
        # Save snapshot with images
        print("Saving snapshot with images...")
        response = requests.post(f"{API_BASE}/simulation/{simulation_id}/snapshot/save?filename=test_snapshot.json&save_images=true")
        response.raise_for_status()
        snapshot_data = response.json()
        print(f"Snapshot saved: {snapshot_data}")
        
        # Check files created
        snapshots_dir = "snapshots"
        if os.path.exists(snapshots_dir):
            files = os.listdir(snapshots_dir)
            print(f"\nFiles in {snapshots_dir}:")
            for file in files:
                filepath = os.path.join(snapshots_dir, file)
                size = os.path.getsize(filepath)
                print(f"  {file} ({size} bytes)")
        
        # Stop simulation
        print("Stopping simulation...")
        requests.post(f"{API_BASE}/simulation/{simulation_id}/stop")
        
        return True
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_snapshot_with_images()
    if success:
        print("\nSnapshots test PASSED!")
    else:
        print("\nSnapshots test FAILED!")
        sys.exit(1)
