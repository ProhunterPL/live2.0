"""
Live 2.0 Performance Tests
Tests for v1 performance requirements
"""

import unittest
import time
import numpy as np
import taichi as ti
from sim.core.particles import ParticleSystem
from sim.core.binding import BindingSystem
from sim.core.stepper import SimulationStepper
from sim.config import SimulationConfig

class TestPerformanceRequirements(unittest.TestCase):
    """Test performance requirements from v1 plan"""
    
    def setUp(self):
        """Set up test environment"""
        self.config = SimulationConfig()
        self.config.max_particles = 1000  # Reduced for testing
        ti.init(arch=ti.cpu)  # Use CPU for consistent testing
    
    def test_particle_creation_performance(self):
        """Test particle creation is fast enough"""
        particles = ParticleSystem(self.config)
        
        start_time = time.time()
        
        # Create many particles
        positions = np.random.uniform(0, 10, (100, 2))
        velocities = np.random.uniform(-1, 1, (100, 2))
        attributes = np.random.uniform(0.5, 2.0, (100, 4))
        
        for i in range(100):
            pos_ti = ti.Vector([positions[i, 0], positions[i, 1]])
            vel_ti = ti.Vector([velocities[i, 0], velocities[i, 1]])
            attr_ti = ti.Vector([
                attributes[i, 0],  # mass
                attributes[i, 1],  # charge_x
                attributes[i, 2],  # charge_y
                attributes[i, 3]   # charge_z
            ])
            
            particles.add_particle_py(pos_ti, vel_ti, attr_ti, 0, 2, 1.0)
        
        creation_time = time.time() - start_time
        
        # Should create 100 particles in less than 1 second
        self.assertLess(creation_time, 1.0)
        self.assertEqual(particles.particle_count[None], 100)
    
    def test_binding_system_performance(self):
        """Test binding system performance"""
        binding = BindingSystem(self.config)
        
        # Create bonds
        start_time = time.time()
        
        num_bonds = 50
        for i in range(num_bonds):
            particle_i = i % 25  # Create bonds between first 25 particles
            particle_j = (i + 1) % 25
            binding.form_bond(particle_i, particle_j)
        
        binding_time = time.time() - start_time
        
        # Should create 50 bonds quickly
        self.assertLess(binding_time, 0.1)
        
        # Test bond retrieval
        bonds = binding.get_bonds()
        self.assertEqual(len(bonds), num_bonds)
    
    def test_graph_computation_performance(self):
        """Test molecular graph computation performance"""
        from sim.core.graphs import MolecularGraph
        
        # Create a moderately complex graph
        n_nodes = 20
        particles = list(range(n_nodes))
        bonds = [(i, i+1) for i in range(n_nodes-1)]  # Linear chain
        
        attributes = {n: np.array([np.random.uniform(0.5, 2.0), 
                                  np.random.uniform(-1, 1), 
                                  np.random.uniform(-1, 1), 
                                  np.random.uniform(-1, 1)]) 
                     for n in range(n_nodes)}
        
        # Time graph creation
        start_time = time.time()
        graph = MolecularGraph(particles, bonds, attributes)
        creation_time = time.time() - start_time
        
        # Time canonical form computation
        start_time = time.time()
        canonical = graph.get_canonical_form()
        canonical_time = time.time() - start_time
        
        # Time complexity computation
        start_time = time.time()
        complexity = graph.get_complexity()
        complexity_time = time.time() - start_time
        
        # Should be fast enough for real-time use
        self.assertLess(creation_time, 0.1)
        self.assertLess(canonical_time, 0.05)
        self.assertLess(complexity_time, 0.01)
        self.assertIsInstance(canonical, str)
        self.assertGreaterEqual(complexity, 0)
    
    def test_catalog_performance(self):
        """Test substance catalog performance"""
        from sim.core.catalog import SubstanceCatalog
        from sim.core.graphs import MolecularGraph
        
        catalog = SubstanceCatalog()
        
        # Create many small graphs
        graphs = []
        for i in range(10):
            n_nodes = np.random.randint(2, 8)
            particles = list(range(n_nodes))
            bonds = [(j, j+1) for j in range(n_nodes-1)]
            attributes = {k: np.array([1.0, 0.0, 0.0, 0.0]) for k in particles}
            graphs.append(MolecularGraph(particles, bonds, attributes))
        
        # Time catalog operations
        start_time = time.time()
        
        for graph in graphs:
            catalog.add_substance(graph, timestamp=time.time())
        
        catalog_time = time.time() - start_time
        
        # Should catalog 10 graphs quickly
        self.assertLess(catalog_time, 0.1)
        
        # Time novelty rate computation
        start_time = time.time()
        novelty_rate = catalog.get_novelty_rate()
        novelty_time = time.time() - start_time
        
        self.assertLess(novelty_time, 0.01)
        self.assertGreaterEqual(novelty_rate, 0)
        self.assertLessEqual(novelty_rate, 1)

if __name__ == '__main__':
    unittest.main()
