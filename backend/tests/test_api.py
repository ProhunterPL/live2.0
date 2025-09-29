"""
Test suite for Live 2.0 API endpoints
"""

import pytest
import asyncio
import json
from fastapi.testclient import TestClient
from api.server import app

client = TestClient(app)

class TestSimulationAPI:
    """Test simulation API endpoints"""
    
    def test_create_simulation_open_chemistry(self):
        """Test creating open chemistry simulation"""
        response = client.post("/simulation/create", json={
            "config": {
                "grid_height": 128,
                "grid_width": 128,
                "mode": "open_chemistry",
                "max_particles": 1000,
                "dt": 0.01,
                "energy_decay": 0.95,
                "energy_threshold": 0.1,
                "particle_radius": 0.5,
                "binding_threshold": 0.8,
                "unbinding_threshold": 0.2,
                "novelty_window": 100,
                "min_cluster_size": 2,
                "vis_frequency": 10,
                "log_frequency": 100,
                "seed": 42
            },
            "mode": "open_chemistry"
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert "simulation_id" in data
        assert data["simulation_id"].startswith("sim_")
        
        return data["simulation_id"]
    
    def test_create_simulation_preset_prebiotic(self):
        """Test creating preset prebiotic simulation"""
        response = client.post("/simulation/create", json={
            "config": {
                "species": {
                    "HCN": 0.1,
                    "NH2CHO": 0.0,
                    "H2O": 0.5
                },
                "reaction_rates": {
                    "HCN_to_NH2CHO": 0.01
                },
                "diffusion_coeffs": {
                    "HCN": 0.1,
                    "NH2CHO": 0.05,
                    "H2O": 0.2
                }
            },
            "mode": "preset_prebiotic"
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert "simulation_id" in data
        
        return data["simulation_id"]
    
    def test_get_simulation_status(self):
        """Test getting simulation status"""
        # First create a simulation
        sim_id = self.test_create_simulation_open_chemistry()
        
        response = client.get(f"/simulation/{sim_id}/status")
        
        assert response.status_code == 200
        data = response.json()
        assert "simulation_id" in data
        assert "is_running" in data
        assert "is_paused" in data
        assert "current_time" in data
        assert "step_count" in data
        assert "particle_count" in data
        assert "novelty_rate" in data
        assert "health_score" in data
    
    def test_start_simulation(self):
        """Test starting simulation"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        response = client.post(f"/simulation/{sim_id}/start")
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert data["message"] == "Simulation started"
    
    def test_pause_simulation(self):
        """Test pausing simulation"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        # Start simulation first
        client.post(f"/simulation/{sim_id}/start")
        
        response = client.post(f"/simulation/{sim_id}/pause")
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert data["message"] == "Simulation paused"
    
    def test_resume_simulation(self):
        """Test resuming simulation"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        # Start and pause simulation first
        client.post(f"/simulation/{sim_id}/start")
        client.post(f"/simulation/{sim_id}/pause")
        
        response = client.post(f"/simulation/{sim_id}/resume")
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert data["message"] == "Simulation resumed"
    
    def test_stop_simulation(self):
        """Test stopping simulation"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        # Start simulation first
        client.post(f"/simulation/{sim_id}/start")
        
        response = client.post(f"/simulation/{sim_id}/stop")
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert data["message"] == "Simulation stopped"
    
    def test_reset_simulation(self):
        """Test resetting simulation"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        response = client.post(f"/simulation/{sim_id}/reset")
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert data["message"] == "Simulation reset"
    
    def test_get_novel_substances(self):
        """Test getting novel substances"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        response = client.get(f"/simulation/{sim_id}/novel-substances?count=5")
        
        assert response.status_code == 200
        data = response.json()
        assert "substances" in data
        assert isinstance(data["substances"], list)
    
    def test_get_metrics(self):
        """Test getting simulation metrics"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        response = client.get(f"/simulation/{sim_id}/metrics")
        
        assert response.status_code == 200
        data = response.json()
        assert "metrics" in data
        metrics = data["metrics"]
        
        # Check required metrics
        required_metrics = [
            "particle_count", "total_energy", "total_mass",
            "bond_count", "cluster_count", "novelty_rate",
            "health_score"
        ]
        
        for metric in required_metrics:
            assert metric in metrics
    
    def test_save_snapshot(self):
        """Test saving snapshot"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        # Test with custom filename
        test_filename = "test_snapshot.json"
        response = client.post(f"/simulation/{sim_id}/snapshot/save", json={
            "filename": test_filename
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert data["filename"] == test_filename
        
        # Verify file was actually written to disk
        import os
        assert os.path.exists(test_filename), f"Snapshot file {test_filename} was not created"
        
        # Verify file contains valid JSON
        import json
        with open(test_filename, 'r') as f:
            snapshot_data = json.load(f)
            assert isinstance(snapshot_data, dict), "Snapshot file should contain valid JSON"
            assert "simulation" in snapshot_data, "Snapshot should contain simulation data"
        
        # Clean up the test file
        os.remove(test_filename)
        assert not os.path.exists(test_filename), f"Failed to clean up test file {test_filename}"
    
    def test_load_snapshot(self):
        """Test loading snapshot"""
        sim_id = self.test_create_simulation_open_chemistry()
        
        # First save a snapshot
        client.post(f"/simulation/{sim_id}/snapshot/save", json={
            "filename": "test_load_snapshot.json"
        })
        
        response = client.post(f"/simulation/{sim_id}/snapshot/load", json={
            "filename": "test_load_snapshot.json"
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] == True
        assert data["message"] == "Snapshot loaded"
    
    def test_invalid_simulation_id(self):
        """Test operations with invalid simulation ID"""
        invalid_id = "invalid_simulation_id"
        
        # Test status endpoint
        response = client.get(f"/simulation/{invalid_id}/status")
        assert response.status_code == 404
        
        # Test start endpoint
        response = client.post(f"/simulation/{invalid_id}/start")
        assert response.status_code == 404
        
        # Test metrics endpoint
        response = client.get(f"/simulation/{invalid_id}/metrics")
        assert response.status_code == 404
    
    def test_invalid_config(self):
        """Test creating simulation with invalid config"""
        response = client.post("/simulation/create", json={
            "config": {
                "grid_height": -1,  # Invalid negative value
                "grid_width": 128,
                "mode": "open_chemistry"
            },
            "mode": "open_chemistry"
        })
        
        assert response.status_code == 400
    
    def test_invalid_mode(self):
        """Test creating simulation with invalid mode"""
        response = client.post("/simulation/create", json={
            "config": {
                "grid_height": 128,
                "grid_width": 128
            },
            "mode": "invalid_mode"
        })
        
        assert response.status_code == 400

class TestWebSocketAPI:
    """Test WebSocket functionality"""
    
    def test_websocket_connection(self):
        """Test WebSocket connection"""
        sim_id = TestSimulationAPI().test_create_simulation_open_chemistry()
        
        with client.websocket_connect(f"/simulation/{sim_id}/stream") as websocket:
            # Should be able to connect
            assert websocket is not None
    
    def test_websocket_invalid_simulation(self):
        """Test WebSocket connection with invalid simulation ID"""
        invalid_id = "invalid_simulation_id"
        
        with pytest.raises(Exception):  # Should raise an exception
            with client.websocket_connect(f"/simulation/{invalid_id}/stream") as websocket:
                pass
    
    def test_websocket_data_format(self):
        """Test WebSocket data format"""
        sim_id = TestSimulationAPI().test_create_simulation_open_chemistry()
        
        # Start simulation to generate data
        client.post(f"/simulation/{sim_id}/start")
        
        with client.websocket_connect(f"/simulation/{sim_id}/stream") as websocket:
            # Wait for data
            data = websocket.receive_bytes()
            
            # Data should be binary (msgpack encoded)
            assert isinstance(data, bytes)
            assert len(data) > 0

class TestErrorHandling:
    """Test error handling"""
    
    def test_malformed_json(self):
        """Test malformed JSON request"""
        response = client.post(
            "/simulation/create",
            data="invalid json",
            headers={"Content-Type": "application/json"}
        )
        
        assert response.status_code == 422  # Unprocessable Entity
    
    def test_missing_required_fields(self):
        """Test missing required fields"""
        response = client.post("/simulation/create", json={
            "mode": "open_chemistry"
            # Missing config field
        })
        
        assert response.status_code == 422
    
    def test_invalid_endpoint(self):
        """Test invalid endpoint"""
        response = client.get("/invalid/endpoint")
        assert response.status_code == 404

class TestConcurrentSimulations:
    """Test multiple concurrent simulations"""
    
    def test_multiple_simulations(self):
        """Test creating multiple simulations"""
        sim_ids = []
        
        # Create multiple simulations
        for i in range(3):
            response = client.post("/simulation/create", json={
                "config": {
                    "grid_height": 64,
                    "grid_width": 64,
                    "mode": "open_chemistry",
                    "max_particles": 100,
                    "dt": 0.01,
                    "energy_decay": 0.95,
                    "energy_threshold": 0.1,
                    "particle_radius": 0.5,
                    "binding_threshold": 0.8,
                    "unbinding_threshold": 0.2,
                    "novelty_window": 100,
                    "min_cluster_size": 2,
                    "vis_frequency": 10,
                    "log_frequency": 100,
                    "seed": 42 + i
                },
                "mode": "open_chemistry"
            })
            
            assert response.status_code == 200
            sim_ids.append(response.json()["simulation_id"])
        
        # All simulations should have different IDs
        assert len(set(sim_ids)) == 3
        
        # All simulations should be accessible
        for sim_id in sim_ids:
            response = client.get(f"/simulation/{sim_id}/status")
            assert response.status_code == 200

if __name__ == "__main__":
    pytest.main([__file__])
