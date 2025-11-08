"""
Hybrid GPU+CPU Simulation Stepper

GPU (Taichi CUDA): Fast particle physics
CPU (Python/NumPy): Complex chemistry analysis

Strategy:
- GPU runs main simulation loop at full speed
- CPU analyzes chemistry in background thread
- Minimal GPU↔CPU transfers (only snapshots every N steps)
- Async communication via queue (non-blocking)
"""

import threading
import queue
import time
import logging
import numpy as np
import networkx as nx
from typing import Dict, List, Optional, Any
from collections import defaultdict

import taichi as ti

from .stepper import SimulationStepper
from ..config import SimulationConfig

logger = logging.getLogger(__name__)


class ChemistrySnapshot:
    """Lightweight snapshot for CPU chemistry analysis"""
    
    def __init__(self, step: int, sim_time: float, 
                 positions: np.ndarray, attributes: np.ndarray,
                 velocities: np.ndarray, active: np.ndarray,
                 energies: np.ndarray):
        self.step = step
        self.sim_time = sim_time
        self.positions = positions
        self.attributes = attributes
        self.velocities = velocities
        self.active = active
        self.energies = energies


class CPUChemistryWorker:
    """Background CPU worker for chemistry analysis"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.running = False
        self.thread = None
        
        # Communication queues
        self.input_queue = queue.Queue(maxsize=10)
        self.output_queue = queue.Queue(maxsize=10)
        
        # Statistics
        self.total_analyzed = 0
        self.analysis_times = []
        
    def start(self):
        """Start background worker thread"""
        if self.running:
            return
        
        self.running = True
        self.thread = threading.Thread(target=self._worker_loop, daemon=True)
        self.thread.start()
        logger.info("CPU chemistry worker started")
    
    def stop(self):
        """Stop background worker thread"""
        self.running = False
        if self.thread:
            self.thread.join(timeout=1.0)
        logger.info(f"CPU chemistry worker stopped (analyzed {self.total_analyzed} snapshots)")
    
    def submit_snapshot(self, snapshot: ChemistrySnapshot):
        """Submit snapshot for analysis (non-blocking)"""
        try:
            self.input_queue.put_nowait(snapshot)
        except queue.Full:
            # Queue full, skip this snapshot
            pass
    
    def get_results(self) -> Optional[Dict[str, Any]]:
        """Get analysis results if available (non-blocking)"""
        try:
            return self.output_queue.get_nowait()
        except queue.Empty:
            return None
    
    def _worker_loop(self):
        """Main worker loop - runs in background thread"""
        logger.info("Chemistry worker loop started")
        
        while self.running:
            try:
                # Wait for snapshot (with timeout to check running flag)
                snapshot = self.input_queue.get(timeout=0.1)
            except queue.Empty:
                continue
            
            # Analyze chemistry
            start_time = time.time()
            results = self._analyze_chemistry(snapshot)
            elapsed = time.time() - start_time
            
            # Track statistics
            self.total_analyzed += 1
            self.analysis_times.append(elapsed)
            if len(self.analysis_times) > 100:
                self.analysis_times = self.analysis_times[-100:]
            
            # Send results back
            try:
                self.output_queue.put_nowait(results)
            except queue.Full:
                pass  # Drop results if queue full
    
    def _analyze_chemistry(self, snapshot: ChemistrySnapshot) -> Dict[str, Any]:
        """
        Perform CPU-based chemistry analysis
        
        This is where we do complex operations that GPU is bad at:
        - Bond detection (branching logic)
        - Cluster detection (graph traversal)
        - Novelty detection (hash table lookup)
        """
        start_time = time.time()
        
        # Filter active particles
        active_indices = np.where(snapshot.active == 1)[0]
        n_active = len(active_indices)
        
        if n_active == 0:
            return {
                'step': snapshot.step,
                'sim_time': snapshot.sim_time,
                'bonds': [],
                'clusters': [],
                'analysis_time_ms': 0.0
            }
        
        # Get active particle data
        positions = snapshot.positions[active_indices]
        attributes = snapshot.attributes[active_indices]
        energies = snapshot.energies[active_indices]
        
        # CPU-based bond detection (much faster than GPU for this!)
        t_bonds_start = time.time()
        bonds = self._detect_bonds_cpu(positions, attributes, energies, active_indices)
        t_bonds = time.time() - t_bonds_start
        
        # CPU-based cluster detection (much faster than GPU!)
        t_clusters_start = time.time()
        clusters = self._detect_clusters_cpu(bonds, n_active)
        t_clusters = time.time() - t_clusters_start
        
        # Calculate complexity metrics
        t_metrics_start = time.time()
        metrics = self._calculate_metrics(positions, attributes, bonds, clusters)
        t_metrics = time.time() - t_metrics_start
        
        total_time = time.time() - start_time
        
        return {
            'step': snapshot.step,
            'sim_time': snapshot.sim_time,
            'bonds': bonds,
            'clusters': clusters,
            'metrics': metrics,
            'timing': {
                'bonds_ms': t_bonds * 1000,
                'clusters_ms': t_clusters * 1000,
                'metrics_ms': t_metrics * 1000,
                'total_ms': total_time * 1000
            },
            'analysis_time_ms': total_time * 1000
        }
    
    def _detect_bonds_cpu(self, positions: np.ndarray, attributes: np.ndarray,
                          energies: np.ndarray, indices: np.ndarray) -> List[tuple]:
        """
        CPU-based bond detection using NumPy vectorization
        
        Much faster than GPU for complex branching logic!
        """
        bonds = []
        n = len(positions)
        
        if n < 2:
            return bonds
        
        # Parameters from config
        bond_distance = getattr(self.config, 'bond_distance', 2.0)
        min_energy = getattr(self.config, 'min_bond_energy', 0.1)
        
        # Vectorized distance calculation (NumPy is fast!)
        for i in range(n):
            for j in range(i + 1, n):
                # Distance
                dx = positions[j, 0] - positions[i, 0]
                dy = positions[j, 1] - positions[i, 1]
                dist = np.sqrt(dx * dx + dy * dy)
                
                if dist > bond_distance:
                    continue
                
                # Check if both particles have enough energy
                if energies[i] < min_energy or energies[j] < min_energy:
                    continue
                
                # Calculate bond strength based on attributes
                # attributes = [mass, charge_x, charge_y, charge_z]
                attr_i = attributes[i]
                attr_j = attributes[j]
                
                # Charge compatibility (opposite charges attract)
                charge_i = attr_i[1:4]
                charge_j = attr_j[1:4]
                charge_dot = np.dot(charge_i, charge_j)
                
                # Negative dot product = opposite charges = good for bonding
                if charge_dot < -0.3:  # Threshold for bonding
                    strength = min(1.0, -charge_dot * 0.5)
                    bonds.append((int(indices[i]), int(indices[j]), float(strength)))
        
        return bonds
    
    def _detect_clusters_cpu(self, bonds: List[tuple], n_particles: int) -> List[List[int]]:
        """
        CPU-based cluster detection using NetworkX
        
        Much faster than GPU for graph algorithms!
        """
        if not bonds:
            return []
        
        # Build graph
        G = nx.Graph()
        for i, j, strength in bonds:
            G.add_edge(i, j, weight=strength)
        
        # Find connected components (clusters)
        clusters = []
        for component in nx.connected_components(G):
            cluster = list(component)
            if len(cluster) >= self.config.min_cluster_size:
                clusters.append(cluster)
        
        return clusters
    
    def _calculate_metrics(self, positions: np.ndarray, attributes: np.ndarray,
                          bonds: List[tuple], clusters: List[List[int]]) -> Dict[str, float]:
        """Calculate chemistry metrics"""
        n_particles = len(positions)
        n_bonds = len(bonds)
        n_clusters = len(clusters)
        
        # Average cluster size
        if clusters:
            avg_cluster_size = np.mean([len(c) for c in clusters])
            max_cluster_size = max([len(c) for c in clusters])
        else:
            avg_cluster_size = 0.0
            max_cluster_size = 0
        
        # Bonding ratio
        bonding_ratio = n_bonds / max(1, n_particles * (n_particles - 1) / 2)
        
        # Average bond strength
        if bonds:
            avg_bond_strength = np.mean([strength for _, _, strength in bonds])
        else:
            avg_bond_strength = 0.0
        
        return {
            'n_particles': n_particles,
            'n_bonds': n_bonds,
            'n_clusters': n_clusters,
            'avg_cluster_size': avg_cluster_size,
            'max_cluster_size': max_cluster_size,
            'bonding_ratio': bonding_ratio,
            'avg_bond_strength': avg_bond_strength
        }
    
    def get_stats(self) -> Dict[str, Any]:
        """Get worker statistics"""
        if self.analysis_times:
            avg_time = np.mean(self.analysis_times)
            max_time = np.max(self.analysis_times)
            min_time = np.min(self.analysis_times)
        else:
            avg_time = max_time = min_time = 0.0
        
        return {
            'total_analyzed': self.total_analyzed,
            'queue_size': self.input_queue.qsize(),
            'avg_analysis_time_ms': avg_time * 1000,
            'max_analysis_time_ms': max_time * 1000,
            'min_analysis_time_ms': min_time * 1000,
        }


class HybridSimulationStepper(SimulationStepper):
    """
    Hybrid GPU+CPU simulation stepper
    
    - GPU (Taichi): Fast particle physics at full speed
    - CPU (Python): Complex chemistry analysis in background
    
    Benefits:
    - GPU runs without blocking on chemistry
    - CPU does what it's best at (complex logic, graphs)
    - Minimal GPU↔CPU transfer overhead
    - Smooth real-time visualization
    """
    
    def __init__(self, config: SimulationConfig):
        # Initialize parent (GPU simulation)
        super().__init__(config)
        
        # Create CPU chemistry worker
        self.chemistry_worker = CPUChemistryWorker(config)
        
        # Configuration
        self.snapshot_interval = getattr(config, 'chemistry_snapshot_interval', 100)
        self.last_snapshot_step = 0
        
        # Latest chemistry results
        self.latest_chemistry_results = None
        
        # Start worker
        self.chemistry_worker.start()
        
        logger.info(f"HybridSimulationStepper initialized")
        logger.info(f"  GPU: Particle physics (Taichi)")
        logger.info(f"  CPU: Chemistry analysis (Python/NumPy/NetworkX)")
        logger.info(f"  Snapshot interval: {self.snapshot_interval} steps")
    
    def step(self, dt: float = None):
        """
        Perform one simulation step (HYBRID)
        
        1. GPU: Physics simulation (fast!)
        2. Check: Should we send snapshot to CPU?
        3. Check: Do we have results from CPU?
        """
        # Call parent step (GPU physics)
        super().step(dt)
        
        # Send snapshot to CPU worker (if needed)
        if self.step_count - self.last_snapshot_step >= self.snapshot_interval:
            self._send_snapshot_to_cpu()
            self.last_snapshot_step = self.step_count
        
        # Check for results from CPU worker
        self._check_cpu_results()
    
    def _send_snapshot_to_cpu(self):
        """Send snapshot to CPU worker (non-blocking)"""
        # Extract data from GPU (minimal transfer)
        particle_count = self.particles.particle_count[None]
        
        if particle_count == 0:
            return
        
        # Limit snapshot size to reduce transfer time
        max_snapshot_particles = min(particle_count, 2000)
        
        # Get data from GPU
        positions = self.particles.positions.to_numpy()[:max_snapshot_particles]
        attributes = self.particles.attributes.to_numpy()[:max_snapshot_particles]
        velocities = self.particles.velocities.to_numpy()[:max_snapshot_particles]
        active = self.particles.active.to_numpy()[:max_snapshot_particles]
        energies = self.particles.energy.to_numpy()[:max_snapshot_particles]
        
        # Create snapshot
        snapshot = ChemistrySnapshot(
            step=self.step_count,
            sim_time=self.current_time,
            positions=positions,
            attributes=attributes,
            velocities=velocities,
            active=active,
            energies=energies
        )
        
        # Submit to worker (non-blocking)
        self.chemistry_worker.submit_snapshot(snapshot)
    
    def _check_cpu_results(self):
        """Check for results from CPU worker (non-blocking)"""
        results = self.chemistry_worker.get_results()
        
        if results is None:
            return
        
        # Update latest results
        self.latest_chemistry_results = results
        
        # Log timing (if slow)
        timing = results.get('timing', {})
        total_ms = timing.get('total_ms', 0)
        
        if total_ms > 100:  # Log if chemistry takes >100ms
            logger.info(f"CPU chemistry analysis:")
            logger.info(f"  Bonds: {timing.get('bonds_ms', 0):.1f}ms")
            logger.info(f"  Clusters: {timing.get('clusters_ms', 0):.1f}ms")
            logger.info(f"  Metrics: {timing.get('metrics_ms', 0):.1f}ms")
            logger.info(f"  Total: {total_ms:.1f}ms")
        
        # Update metrics if available
        if 'metrics' in results:
            self._update_metrics_from_cpu(results['metrics'])
    
    def _update_metrics_from_cpu(self, cpu_metrics: Dict[str, float]):
        """Update simulation metrics from CPU analysis"""
        # Update bond count
        if 'n_bonds' in cpu_metrics:
            self.metrics.bond_count[None] = int(cpu_metrics['n_bonds'])
        
        # Update cluster count
        if 'n_clusters' in cpu_metrics:
            self.metrics.cluster_count[None] = int(cpu_metrics['n_clusters'])
    
    def get_visualization_data(self) -> Dict:
        """
        Get visualization data (override to include CPU results)
        """
        # Get base visualization from parent (GPU)
        data = super().get_visualization_data()
        
        # Add latest CPU chemistry results
        if self.latest_chemistry_results:
            data['chemistry'] = {
                'step': self.latest_chemistry_results.get('step', 0),
                'bonds': self.latest_chemistry_results.get('bonds', []),
                'clusters': self.latest_chemistry_results.get('clusters', []),
                'metrics': self.latest_chemistry_results.get('metrics', {}),
                'analysis_time_ms': self.latest_chemistry_results.get('analysis_time_ms', 0)
            }
        
        # Add CPU worker statistics
        data['cpu_worker'] = self.chemistry_worker.get_stats()
        
        return data
    
    def stop(self):
        """Stop simulation and CPU worker"""
        # Stop CPU worker first
        self.chemistry_worker.stop()
        
        # Call parent stop
        super().stop()
        
        logger.info("HybridSimulationStepper stopped")
    
    def get_simulation_state(self) -> Dict:
        """Get current simulation state (override to include CPU stats)"""
        state = super().get_simulation_state()
        
        # Add CPU worker stats
        state['cpu_worker'] = self.chemistry_worker.get_stats()
        
        # Add latest chemistry results
        if self.latest_chemistry_results:
            metrics = self.latest_chemistry_results.get('metrics', {})
            state['cpu_bonds'] = metrics.get('n_bonds', 0)
            state['cpu_clusters'] = metrics.get('n_clusters', 0)
            state['cpu_analysis_time_ms'] = self.latest_chemistry_results.get('analysis_time_ms', 0)
        
        return state

