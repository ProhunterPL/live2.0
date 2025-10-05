"""
Live 2.0 Unit Tests
Basic unit tests for core functionality
"""

import unittest
import numpy as np
import taichi as ti
from sim.core.graphs import MolecularGraph, GraphCatalog
from sim.core.catalog import SubstanceCatalog, SubstanceRecord
from sim.config import SimulationConfig, PresetPrebioticConfig

class TestMolecularGraph(unittest.TestCase):
    """Test molecular graph functionality"""
    
    def test_graph_creation(self):
        """Test basic graph creation"""
        particles = [0, 1, 2]
        bonds = [(0, 1), (1, 2)]
        attributes = {0: np.array([1.0, 0.0, 0.0, 0.0]), 
                     1: np.array([1.0, 0.0, 0.0, 0.0]), 
                     2: np.array([1.0, 0.0, 0.0, 0.0])}
        
        graph = MolecularGraph(particles, bonds, attributes)
        
        self.assertEqual(graph.get_node_count(), 3)
        self.assertEqual(graph.get_edge_count(), 2)
        self.assertEqual(graph.get_density(), 2/3)  # 2 edges / 3 possible edges
    
    def test_graph_hash_stability(self):
        """Test that identical graphs have same hash"""
        particles1 = [0, 1, 2]
        bonds1 = [(0, 1), (1, 2)]
        attributes1 = {0: np.array([1.0, 0.0, 0.0, 0.0]), 
                      1: np.array([1.0, 0.0, 0.0, 0.0]), 
                      2: np.array([1.0, 0.0, 0.0, 0.0])}
        
        particles2 = [2, 0, 1]  # Same particles, different order
        bonds2 = [(2, 0), (0, 1)]  # Equivalent bonds
        attributes2 = {2: np.array([1.0, 0.0, 0.0, 0.0]), 
                       0: np.array([1.0, 0.0, 0.0, 0.0]), 
                       1: np.array([1.0, 0.0, 0.0, 0.0])}
        
        graph1 = MolecularGraph(particles1, bonds1, attributes1)
        graph2 = MolecularGraph(particles2, bonds2, attributes2)
        
        hash1 = graph1.get_canonical_form()
        hash2 = graph2.get_canonical_form()
        
        # Graphs should be identical despite different node ordering
        self.assertEqual(hash1, hash2)
    
    def test_graph_complexity(self):
        """Test complexity calculation"""
        # Simple chain
        particles = [0, 1, 2]
        bonds = [(0, 1), (1, 2)]
        attributes = {0: np.array([1.0, 0.0, 0.0, 0.0]), 
                     1: np.array([1.0, 0.0, 0.0, 0.0]), 
                     2: np.array([1.0, 0.0, 0.0, 0.0])}
        
        graph = MolecularGraph(particles, bonds, attributes)
        complexity = graph.get_complexity()
        
        self.assertGreater(complexity, 0)
        self.assertIsInstance(complexity, float)

class TestSubstanceCatalog(unittest.TestCase):
    """Test substance catalog functionality"""
    
    def test_catalog_basic_operations(self):
        """Test basic catalog operations"""
        catalog = SubstanceCatalog()
        
        # Create a simple graph
        particles = [0, 1, 2]
        bonds = [(0, 1), (1, 2)]
        attributes = {0: np.array([1.0, 0.0, 0.0, 0.0]), 
                     1: np.array([1.0, 0.0, 0.0, 0.0]), 
                     2: np.array([1.0, 0.0, 0.0, 0.0])}
        
        graph = MolecularGraph(particles, bonds, attributes)
        
        # Add substance
        is_novel, substance_id = catalog.add_substance(graph)
        
        self.assertTrue(is_novel)
        self.assertIsInstance(substance_id, str)
        self.assertEqual(catalog.total_discoveries, 1)
        self.assertEqual(catalog.novel_discoveries, 1)
        
        # Add same substance again
        is_novel2, substance_id2 = catalog.add_substance(graph)
        
        self.assertFalse(is_novel2)  # Should not be novel second time
        self.assertEqual(substance_id, substance_id2)  # Same ID
        self.assertEqual(catalog.total_discoveries, 2)
        self.assertEqual(catalog.novel_discoveries, 1)  # Still 1
        
        # Test novelty rate
        novelty_rate = catalog.get_novelty_rate()
        self.assertEqual(novelty_rate, 0.5)  # 1 novel out of 2 total
    
    def test_catalog_statistics(self):
        """Test catalog statistics"""
        catalog = SubstanceCatalog()
        
        # Add several substances
        for i in range(5):
            particles = list(range(i + 2))  # Variable size
            bonds = [(j, j + 1) for j in range(len(particles) - 1)]
            attributes = {k: np.array([1.0, 0.0, 0.0, 0.0]) for k in particles}
            
            graph = MolecularGraph(particles, bonds, attributes)
            catalog.add_substance(graph)
        
        stats = catalog.get_catalog_stats()
        
        self.assertEqual(stats['total_substances'], 5)
        self.assertEqual(stats['total_novel'], 5)
        self.assertEqual(stats['total_discovered'], 5)
        self.assertEqual(stats['novelty_rate'], 1.0)  # All novel
        self.assertGreater(stats['runtime_hours'], 0)

class TestConfiguration(unittest.TestCase):
    """Test configuration classes"""
    
    def test_simulation_config_defaults(self):
        """Test default simulation configuration"""
        config = SimulationConfig()
        
        self.assertEqual(config.grid_height, 256)
        self.assertEqual(config.grid_width, 256)
        self.assertEqual(config.mode, "open_chemistry")
        self.assertEqual(config.max_particles, 10000)
        self.assertTrue(config.grid_height > 0)
        self.assertTrue(config.grid_width > 0)
        self.assertTrue(config.max_particles > 0)
    
    def test_preset_config_defaults(self):
        """Test preset configuration"""
        config = PresetPrebioticConfig()
        
        self.assertIn("HCN", config.species)
        self.assertIn("NH2CHO", config.species)
        self.assertIn("H2O", config.species)
        self.assertIn("HCN_to_NH2CHO", config.reaction_rates)
        self.assertTrue(config.species["HCN"] > 0)

class TestGraphCatalog(unittest.TestCase):
    """Test graph catalog functionality"""
    
    def test_graph_catalog_basic(self):
        """Test basic graph catalog operations"""
        catalog = GraphCatalog()
        
        # Create two graphs
        graph1 = MolecularGraph([0, 1], [(0, 1)], {0: np.array([1.0]), 1: np.array([1.0])})
        graph2 = MolecularGraph([0, 1], [(0, 1)], {0: np.array([1.0]), 1: np.array([1.0])})
        graph3 = MolecularGraph([0, 1, 2], [(0, 1), (1, 2)], 
                                {0: np.array([1.0]), 1: np.array([1.0]), 2: np.array([1.0])})
        
        # Add graphs
        novel1 = catalog.add_graph(graph1, timestamp=0.0)
        novel2 = catalog.add_graph(graph2, timestamp=1.0)
        novel3 = catalog.add_graph(graph3, timestamp=2.2)
        
        self.assertTrue(novel1)  # First graph should be novel
        self.assertFalse(novel2)  # Same as graph1
        self.assertTrue(novel3)  # Different structure
        
        self.assertEqual(catalog.total_graphs_seen, 3)
        self.assertEqual(catalog.novel_graphs_count, 2)
        
        # Test novelty rate
        novelty_rate = catalog.get_novelty_rate()
        self.assertGreaterEqual(novelty_rate, 0)
        self.assertLessEqual(novelty_rate, 1)

if __name__ == '__main__':
    unittest.main()
