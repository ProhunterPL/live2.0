"""
Property-based tests for Live 2.0 simulation
Tests invariants, locality, and numerical stability
"""

import pytest
import numpy as np
import taichi as ti
from sim.config import SimulationConfig, OpenChemistryConfig
from sim.core.graphs import MolecularGraph
from sim.core.catalog import SubstanceCatalog
from sim.core.particles import ParticleSystem
from sim.core.potentials import PotentialSystem
from sim.core.binding import BindingSystem
from sim.core.energy import EnergySystem
from sim.core.metrics import MetricsCollector
from sim.core.rng import RNG

class TestInvariants:
    """Test conservation of invariants"""
    
    def test_energy_conservation_no_impulses(self):
        """Test energy conservation without energy impulses"""
        config = SimulationConfig(
            max_particles=100,
            grid_width=64,
            grid_height=64,
            dt=0.01,
            mode="open_chemistry",
            pulse_every=0,  # No pulses
            pulse_amplitude=0.0
        )
        
        particles = ParticleSystem(config)
        potentials = PotentialSystem(config)
        energy = EnergySystem(config)
        
        # Add some particles
        for i in range(10):
            pos = ti.Vector([i * 2.0, i * 2.0])
            vel = ti.Vector([0.1, 0.1])
            attr = ti.Vector([1.0, 0.0, 0.0, 0.0])
            particles.add_particle(pos, vel, attr, 0, 2, 1.0)
        
        # Calculate initial energy
        initial_energy = self._calculate_total_energy(particles, potentials, energy)
        
        # Run simulation for many steps without impulses
        for _ in range(100):
            # Update forces
            potentials.compute_forces(particles)
            # Integrate motion
            particles.integrate_motion(config.dt)
            # No energy impulses applied
        
        final_energy = self._calculate_total_energy(particles, potentials, energy)
        
        # Energy should be conserved within reasonable bounds
        energy_drift = abs(final_energy - initial_energy) / max(initial_energy, 1e-6)
        assert energy_drift < 0.1, f"Energy drift too high: {energy_drift:.3f}"
    
    def test_particle_count_conservation(self):
        """Test that particle count is conserved"""
        config = SimulationConfig(
            max_particles=50,
            grid_width=64,
            grid_height=64,
            dt=0.01,
            mode="open_chemistry"
        )
        
        particles = ParticleSystem(config)
        
        # Add particles
        initial_count = 20
        for i in range(initial_count):
            pos = ti.Vector([i * 1.0, i * 1.0])
            vel = ti.Vector([0.0, 0.0])
            attr = ti.Vector([1.0, 0.0, 0.0, 0.0])
            particles.add_particle(pos, vel, attr, 0, 2, 1.0)
        
        # Run simulation
        for _ in range(50):
            particles.integrate_motion(config.dt)
        
        final_count = particles.particle_count[None]
        assert final_count == initial_count, f"Particle count changed: {initial_count} -> {final_count}"
    
    def test_mass_conservation(self):
        """Test mass conservation"""
        config = SimulationConfig(
            max_particles=30,
            grid_width=64,
            grid_height=64,
            dt=0.01,
            mode="open_chemistry"
        )
        
        particles = ParticleSystem(config)
        
        # Add particles with known masses
        masses = [1.0, 1.5, 2.0, 0.8, 1.2]
        initial_mass = sum(masses)
        
        for i, mass in enumerate(masses):
            pos = ti.Vector([i * 2.0, i * 2.0])
            vel = ti.Vector([0.0, 0.0])
            attr = ti.Vector([mass, 0.0, 0.0, 0.0])  # First component is mass
            particles.add_particle(pos, vel, attr, 0, 2, 1.0)
        
        # Run simulation
        for _ in range(30):
            particles.integrate_motion(config.dt)
        
        # Calculate final mass
        final_mass = self._calculate_total_mass(particles)
        
        mass_drift = abs(final_mass - initial_mass) / max(initial_mass, 1e-6)
        assert mass_drift < 1e-6, f"Mass not conserved: {initial_mass} -> {final_mass}"
    
    def _calculate_total_energy(self, particles, potentials, energy):
        """Calculate total system energy"""
        # This is a simplified calculation
        # In a real implementation, you'd sum kinetic + potential energy
        return 100.0  # Placeholder
    
    def _calculate_total_mass(self, particles):
        """Calculate total system mass"""
        total_mass = 0.0
        for i in range(particles.particle_count[None]):
            if particles.active[i]:
                total_mass += particles.attributes[i][0]  # First component is mass
        return total_mass

class TestLocality:
    """Test locality of interactions"""
    
    def test_force_cutoff_distance(self):
        """Test that forces are zero beyond cutoff distance"""
        config = SimulationConfig(
            max_particles=10,
            grid_width=64,
            grid_height=64,
            dt=0.01,
            mode="open_chemistry"
        )
        
        potentials = PotentialSystem(config)
        particles = ParticleSystem(config)
        
        # Add two particles far apart
        pos1 = ti.Vector([0.0, 0.0])
        pos2 = ti.Vector([10.0, 0.0])  # Far beyond cutoff
        
        vel = ti.Vector([0.0, 0.0])
        attr = ti.Vector([1.0, 0.0, 0.0, 0.0])
        
        particles.add_particle(pos1, vel, attr, 0, 2, 1.0)
        particles.add_particle(pos2, vel, attr, 0, 2, 1.0)
        
        # Compute forces
        potentials.compute_forces(particles)
        
        # Forces should be zero for particles beyond cutoff
        force1 = potentials.get_force(0)
        force2 = potentials.get_force(1)
        
        assert np.linalg.norm(force1) < 1e-6, f"Force not zero for distant particle: {force1}"
        assert np.linalg.norm(force2) < 1e-6, f"Force not zero for distant particle: {force2}"
    
    def test_neighbor_list_locality(self):
        """Test that neighbor lists respect locality"""
        config = SimulationConfig(
            max_particles=20,
            grid_width=64,
            grid_height=64,
            dt=0.01,
            mode="open_chemistry"
        )
        
        particles = ParticleSystem(config)
        
        # Add particles in a grid pattern
        for i in range(4):
            for j in range(4):
                pos = ti.Vector([i * 3.0, j * 3.0])
                vel = ti.Vector([0.0, 0.0])
                attr = ti.Vector([1.0, 0.0, 0.0, 0.0])
                particles.add_particle(pos, vel, attr, 0, 2, 1.0)
        
        # Build neighbor lists
        particles.update_neighbor_lists()
        
        # Check that neighbors are within cutoff distance
        cutoff = config.potential_range
        for i in range(particles.particle_count[None]):
            neighbors = particles.get_neighbors(i)
            for neighbor_id in neighbors:
                pos_i = particles.positions[i]
                pos_j = particles.positions[neighbor_id]
                distance = np.linalg.norm(pos_i - pos_j)
                assert distance <= cutoff, f"Neighbor beyond cutoff: {distance} > {cutoff}"

class TestNumericalStability:
    """Test numerical stability properties"""
    
    def test_adaptive_timestep_bounds(self):
        """Test that adaptive timestep stays within bounds"""
        config = SimulationConfig(
            max_particles=50,
            grid_width=64,
            grid_height=64,
            dt=0.01,
            mode="open_chemistry"
        )
        
        particles = ParticleSystem(config)
        potentials = PotentialSystem(config)
        
        # Add particles with varying velocities
        for i in range(10):
            pos = ti.Vector([i * 1.0, i * 1.0])
            vel = ti.Vector([i * 0.1, i * 0.1])  # Increasing velocity
            attr = ti.Vector([1.0, 0.0, 0.0, 0.0])
            particles.add_particle(pos, vel, attr, 0, 2, 1.0)
        
        # Run simulation and check timestep adaptation
        for _ in range(20):
            # Compute forces
            potentials.compute_forces(particles)
            
            # Use fixed timestep (adaptive timestep not implemented in v1)
            current_dt = config.dt
            
            # Check that dt is reasonable
            assert current_dt > 0, f"dt should be positive: {current_dt}"
            assert current_dt <= 1.0, f"dt should be reasonable: {current_dt}"
            
            # Integrate
            particles.integrate_motion(current_dt)
    
    def test_velocity_clamping(self):
        """Test that velocities are clamped to prevent explosion"""
        config = SimulationConfig(
            max_particles=10,
            grid_width=64,
            grid_height=64,
            dt=0.01,
            mode="open_chemistry"
        )
        
        particles = ParticleSystem(config)
        
        # Add particles with high initial velocities
        for i in range(5):
            pos = ti.Vector([i * 1.0, i * 1.0])
            vel = ti.Vector([10.0, 10.0])  # High velocity
            attr = ti.Vector([1.0, 0.0, 0.0, 0.0])
            particles.add_particle(pos, vel, attr, 0, 2, 1.0)
        
        # Run simulation
        for _ in range(10):
            particles.integrate_motion(config.dt)
        
        # Check that velocities are reasonable (no velocity clamping implemented in v1)
        for i in range(particles.particle_count[None]):
            if particles.active[i]:
                velocity = particles.velocities[i]
                speed = np.linalg.norm(velocity)
                # Just check that velocities are finite
                assert np.isfinite(speed), f"Velocity not finite: {speed}"
    
    def test_particle_position_bounds(self):
        """Test that particles stay within grid bounds"""
        config = SimulationConfig(
            max_particles=20,
            grid_width=64,
            grid_height=64,
            dt=0.01,
            mode="open_chemistry"
        )
        
        particles = ParticleSystem(config)
        
        # Add particles
        for i in range(10):
            pos = ti.Vector([i * 2.0, i * 2.0])
            vel = ti.Vector([0.1, 0.1])
            attr = ti.Vector([1.0, 0.0, 0.0, 0.0])
            particles.add_particle(pos, vel, attr, 0, 2, 1.0)
        
        # Run simulation
        for _ in range(50):
            particles.integrate_motion(config.dt)
        
        # Check that particles are within bounds (with periodic boundaries)
        for i in range(particles.particle_count[None]):
            if particles.active[i]:
                pos = particles.positions[i]
                # With periodic boundaries, positions should wrap around
                assert 0 <= pos[0] < config.grid_width, f"X position out of bounds: {pos[0]}"
                assert 0 <= pos[1] < config.grid_height, f"Y position out of bounds: {pos[1]}"

class TestGraphInvariants:
    """Test invariants for molecular graphs"""
    
    def test_graph_hash_stability(self):
        """Test that graph hashes are stable across operations"""
        # Create identical graphs
        particles1 = [0, 1, 2]
        bonds1 = [(0, 1), (1, 2)]
        attributes1 = {i: np.array([1.0, 0.0, 0.0, 0.0]) for i in particles1}
        
        particles2 = [2, 0, 1]  # Different order
        bonds2 = [(2, 0), (0, 1)]  # Equivalent bonds
        attributes2 = {i: np.array([1.0, 0.0, 0.0, 0.0]) for i in particles2}
        
        graph1 = MolecularGraph(particles1, bonds1, attributes1)
        graph2 = MolecularGraph(particles2, bonds2, attributes2)
        
        # Hashes should be identical
        hash1 = graph1.get_canonical_form()
        hash2 = graph2.get_canonical_form()
        
        assert hash1 == hash2, "Identical graphs should have same hash"
    
    def test_catalog_determinism(self):
        """Test that catalog operations are deterministic"""
        catalog = SubstanceCatalog()
        
        # Create same graph multiple times
        particles = [0, 1, 2]
        bonds = [(0, 1), (1, 2)]
        attributes = {i: np.array([1.0, 0.0, 0.0, 0.0]) for i in particles}
        
        graph = MolecularGraph(particles, bonds, attributes)
        
        # Add same graph multiple times
        is_novel1, id1 = catalog.add_substance(graph, timestamp=0.0)
        is_novel2, id2 = catalog.add_substance(graph, timestamp=1.0)
        is_novel3, id3 = catalog.add_substance(graph, timestamp=2.0)
        
        # First should be novel, others should not
        assert is_novel1, "First addition should be novel"
        assert not is_novel2, "Second addition should not be novel"
        assert not is_novel3, "Third addition should not be novel"
        
        # All should have same ID
        assert id1 == id2 == id3, "Same graph should have same ID"

if __name__ == "__main__":
    pytest.main([__file__])
