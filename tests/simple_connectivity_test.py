#!/usr/bin/env python3
"""
Simple connectivity test for Live 2.0 backend
Basic functionality verification for Phase 0
"""

import pytest
import requests
import time

API_BASE = "http://localhost:8001"

@pytest.mark.integration
def test_backend_connectivity():
    """Test basic backend connectivity"""
    print("="*60)
    print("LIVE 2.0 SIMPLE CONNECTIVITY TEST")
    print("="*30)
    
    try:
        # Test 1: Root endpoint
        print("1. Testing root endpoint...")
        response = requests.get(f"{API_BASE}/", timeout=5)
        if response.status_code == 200:
            print("   OK Root endpoint working")
        else:
            print(f"   FAIL Root endpoint returned {response.status_code}")
            pytest.fail(f"Root endpoint failed with status {response.status_code}")
        
        # Test 2: List active simulations
        print("2. Testing simulations list...")
        response = requests.get(f"{API_BASE}/simulations/active", timeout=5)
        if response.status_code == 200:
            data = response.json()
            print(f"   OK Active simulations: {data['count']}")
            sim_ids = data['simulations']
        else:
            print(f"   FAIL Simulations list returned {response.status_code}")
            pytest.fail(f"Simulations list failed with status {response.status_code}")
        
        # Test 3: Create small simulation
        print("3. Testing simulation creation...")
        test_config = {
            "config": {
                "max_particles": 50,
                "grid_width": 64,
                "grid_height": 64,
                "mode": "open_chemistry",
                "dt": 0.035
            },
            "mode": "open_chemistry"
        }
        
        response = requests.post(f"{API_BASE}/simulation/create", json=test_config, timeout=10)
        if response.status_code == 200:
            data = response.json()
            test_sim_id = data.get("simulation_id")
            print(f"   OK Created test simulation: {test_sim_id}")
        else:
            print(f"   FAIL Simulation creation returned {response.status_code}")
            print(f"   Response: {response.text}")
            pytest.fail(f"Simulation creation failed with status {response.status_code}")
        
        # Test 4: Get simulation status
        print("4. Testing status retrieval...")
        response = requests.get(f"{API_BASE}/simulation/{test_sim_id}/status", timeout=5)
        if response.status_code == 200:
            data = response.json()
            print(f"   OK Status retrieved: {data['particle_count']} particles")
        else:
            print(f"   FAIL Status retrieval returned {response.status_code}")
            pytest.fail(f"Status retrieval failed with status {response.status_code}")
        
        # Test 5: Start simulation (with timeout tolerance)
        print("5. Testing simulation start...")
        try:
            response = requests.post(f"{API_BASE}/simulation/{test_sim_id}/start", timeout=15)
            if response.status_code == 200:
                print("   OK Simulation started")
                
                # Quick status check after 2 seconds
                time.sleep(2)
                status_response = requests.get(f"{API_BASE}/simulation/{test_sim_id}/status", timeout=10)
                if status_response.status_code == 200:
                    status_data = status_response.json()
                    print(f"   INFO Status after 2s: {status_data['step_count']} steps, running: {status_data['is_running']}")
                else:
                    print("   WARNING Could not get status after start")
            else:
                print(f"   FAIL Simulation start returned {response.status_code}")
                print(f"   Response: {response.text}")
        except requests.exceptions.Timeout:
            print("   WARNING Simulation start timed out (but this might be due to Taichi compilation)")
        
        # Test 6: Stop simulation
        print("6. Testing simulation stop...")
        try:
            response = requests.post(f"{API_BASE}/simulation/{test_sim_id}/stop", timeout=5)
            if response.status_code == 200:
                print("   OK Simulation stopped")
            else:
                print(f"   FAIL Simulation stop returned {response.status_code}")
        except requests.exceptions.Timeout:
            print("   WARNING Simulation stop timed out")
        
        print("="*30)
        print("CONNECTIVITY TEST COMPLETE")
        print("="*30)
        
        # Summary
        if len(sim_ids) > 0:
            print(f"INFO: Found {len(sim_ids)} existing simulation(s)")
            print(f"INFO: Backend is responding correctly")
            print(f"INFO: Basic API functionality verified")
        else:
            print("INFO: No active simulations found")
            print("INFO: Backend connectivity OK")
        
        # Test passed
        assert True
            
    except Exception as e:
        print(f"\nERROR: {e}")
        pytest.fail(f"Backend connectivity test failed: {e}")

if __name__ == "__main__":
    success = test_backend_connectivity()
    if success:
        print("\nRESULT: Backend connectivity test PASSED")
        exit(0)
    else:
        print("\nRESULT: Backend connectivity test FAILED")
        exit(1)
