"""
Live 2.0 Snapshot Tests
Tests for snapshot save/load functionality
"""

import unittest
import tempfile
import os
import json
from sim.io.snapshot import SnapshotManager, SnapshotSerializer, SnapshotAPI
from sim.core.catalog import SubstanceCatalog
from sim.core.graphs import MolecularGraph

class TestSnapshotManager(unittest.TestCase):
    """Test snapshot management functionality"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.snapshot_manager = SnapshotManager(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_snapshot_creation(self):
        """Test snapshot creation"""
        # Create mock simulation data
        simulation_data = {
            "config": {"mode": "open_chemistry", "grid_width": 128, "grid_height": 128},
            "current_time": 100.5,
            "step_count": 5000,
            "particles": {
                "positions": [[10.0, 20.0], [30.0, 40.0]],
                "attributes": [[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]],
                "active_mask": [1, 1]
            },
            "energy_field": [[0.1, 0.2], [0.3, 0.4]],
            "bonds": [[0, 1, 0.5]],
            "clusters": [[0, 1]],
            "metrics": {"particle_count": 2, "bond_count": 1}
        }
        
        filename = self.snapshot_manager.create_snapshot(simulation_data, "test_snapshot")
        
        # Check that file was created
        self.assertIsInstance(filename, str)
        self.assertTrue(filename.endswith('.json'))
        
        # Check metadata was updated
        snapshots = self.snapshot_manager.list_snapshots()
        self.assertEqual(len(snapshots), 1)
        self.assertEqual(snapshots[0]['name'], 'test_snapshot')
        self.assertEqual(snapshots[0]['simulation_mode'], 'open_chemistry')
        self.assertEqual(snapshots[0]['particle_count'], 2)
    
    def test_snapshot_loading(self):
        """Test snapshot loading"""
        # First create a snapshot
        simulation_data = {
            "config": {"mode": "open_chemistry", "grid_width": 128, "grid_height": 128},
            "current_time": 50.0,
            "step_count": 2500,
            "particles": {
                "positions": [[5.0, 10.0]],
                "attributes": [[2.0, 0.5, 0.0, 0.0]],
                "active_mask": [1]
            },
            "energy_field": [[0.5]],
            "bonds": [],
            "clusters": [],
            "metrics": {"particle_count": 1}
        }
        
        self.snapshot_manager.create_snapshot(simulation_data, "load_test")
        
        # Now load it back
        loaded_data = self.snapshot_manager.load_snapshot("load_test")
        
        self.assertEqual(loaded_data["current_time"], 50.0)
        self.assertEqual(loaded_data["step_count"], 2500)
        self.assertEqual(len(loaded_data["particles"]["positions"]), 1)
        self.assertEqual(loaded_data["particles"]["positions"][0], [5.0, 10.0])
    
    def test_snapshot_validation(self):
        """Test snapshot validation"""
        # Valid snapshot
        valid_data = {
            "config": {"mode": "open_chemistry"},
            "current_time": 0.0,
            "step_count": 0,
            "mode": "open_chemistry",
            "particles": {
                "positions": [[1.0, 2.0]],
                "attributes": [[1.0, 0.0, 0.0, 0.0]],
                "active_mask": [1]
            },
            "energy_field": [[0.0]],
            "bonds": [],
            "clusters": []
        }
        
        isValid = SnapshotSerializer.validate_snapshot(valid_data)
        self.assertTrue(isValid)
        
        # Invalid snapshot (missing field)
        invalid_data = valid_data.copy()
        del invalid_data["step_count"]
        
        isInvalid = SnapshotSerializer.validate_snapshot(invalid_data)
        self.assertFalse(isInvalid)
        
        # Invalid snapshot (inconsistent arrays)
        invalid_data2 = valid_data.copy()
        invalid_data2["particles"]["positions"] = [[1.0, 2.0]]
        invalid_data2["particles"]["attributes"] = [[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]]
        invalid_data2["particles"]["active_mask"] = [1]
        
        isInvalid2 = SnapshotSerializer.validate_snapshot(invalid_data2)
        self.assertFalse(isInvalid2)
    
    def test_snapshot_deletion(self):
        """Test snapshot deletion"""
        # Create multiple snapshots
        for i in range(3):
            data = {"config": {}, "current_time": i, "step_count": 100 * i,
                   "particles": {"positions": [], "attributes": [], "active_mask": []},
                   "energy_field": [[0.0]], "bonds": [], "clusters": []}
            self.snapshot_manager.create_snapshot(data, f"snapshot_{i}")
        
        # Check snapshots exist
        snapshots = self.snapshot_manager.list_snapshots()
        self.assertEqual(len(snapshots), 3)
        
        # Delete one
        success = self.snapshot_manager.delete_snapshot("snapshot_1")
        self.assertTrue(success)
        
        # Check snapshot was removed
        snapshots = self.snapshot_manager.list_snapshots()
        self.assertEqual(len(snapshots), 2)
        
        # Check remaining snapshots
        names = [s['name'] for s in snapshots]
        self.assertIn('snapshot_0', names)
        self.assertNotIn('snapshot_1', names)
        self.assertIn('snapshot_2', names)

class TestSnapshotAPI(unittest.TestCase):
    """Test snapshot API functionality"""
    
    def setUp(self):
        """Set up test environment"""
        self temp_dir = tempfile.mkdtemp()
        self.snapshot_manager = SnapshotManager(self.temp_dir)
        self.snapshot_api = SnapshotAPI(self.snapshot_manager)
    
    def tearDown(self):
        """Clean up test environment"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_snapshot_api_operations(self):
        """Test snapshot API operations"""
        # Create mock simulation object
        class MockSimulation:
            def __init__(self):
                self.config = {"mode": "test"}
                self.current_time = 25.0
                self.step_count = 1000
        
        simulation = MockSimulation()
        
        # Test snapshot creation through API
        filename = self.snapshot_api.create_snapshot(simulation, "api_test")
        
        self.assertIsInstance(filename, str)
        
        # Test listing snapshots
        snapshots = self.snapshot_api.list_snapshots()
        self.assertEqual(len(snapshots), 1)
        self.assertEqual(snapshots[0]['name'], 'api_test')
        
        # Test snapshot info
        info = self.snapshot_api.get_snapshot_info('api_test')
        self.assertIsNotNone(info)
        self.assertEqual(info['name'], 'api_test')
        
        # Test snapshot deletion
        success = self.snapshot_api.delete_snapshot('api_test')
        self.assertTrue(success)

if __name__ == '__main__':
    unittest.main()
