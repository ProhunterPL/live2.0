"""
Metrics system for Live 2.0 simulation
Tracks novelty, complexity, and other key metrics
"""

import taichi as ti
import numpy as np
import time
from typing import Dict, List, Tuple, Optional
from collections import deque
from .catalog import SubstanceCatalog

# Compile-time constants
MAX_PARTICLES_COMPILE = 10000
GRID_WIDTH_COMPILE = 256
GRID_HEIGHT_COMPILE = 256

# Global Taichi fields
particle_count_field = None
total_energy_field = None
total_mass_field = None
bond_count_field = None
cluster_count_field = None
energy_field_sum_field = None
energy_field_max_field = None
energy_field_mean_field = None

def init_metrics_fields():
    """Initialize global Taichi fields for metrics"""
    global particle_count_field, total_energy_field, total_mass_field
    global bond_count_field, cluster_count_field, energy_field_sum_field
    global energy_field_max_field, energy_field_mean_field
    
    particle_count_field = ti.field(dtype=ti.i32, shape=())
    total_energy_field = ti.field(dtype=ti.f32, shape=())
    total_mass_field = ti.field(dtype=ti.f32, shape=())
    bond_count_field = ti.field(dtype=ti.i32, shape=())
    cluster_count_field = ti.field(dtype=ti.i32, shape=())
    energy_field_sum_field = ti.field(dtype=ti.f32, shape=())
    energy_field_max_field = ti.field(dtype=ti.f32, shape=())
    energy_field_mean_field = ti.field(dtype=ti.f32, shape=())

# Module-level kernels
@ti.kernel
def reset_metrics_kernel():
    """Reset all metrics - module-level kernel"""
    particle_count_field[None] = 0
    total_energy_field[None] = 0.0
    total_mass_field[None] = 0.0
    bond_count_field[None] = 0
    cluster_count_field[None] = 0
    energy_field_sum_field[None] = 0.0
    energy_field_max_field[None] = 0.0
    energy_field_mean_field[None] = 0.0

@ti.kernel
def update_particle_metrics_kernel(active: ti.template(), attributes: ti.template(),
                                  energy: ti.template(), particle_count: ti.i32):
    """Update particle metrics - module-level kernel"""
    particle_count_field[None] = particle_count
    
    total_energy = 0.0
    total_mass = 0.0
    
    for i in range(MAX_PARTICLES_COMPILE):
        if active[i] == 1:
            total_energy += energy[i]
            total_mass += attributes[i][0]
    
    total_energy_field[None] = total_energy
    total_mass_field[None] = total_mass

@ti.kernel
def update_bond_metrics_kernel(bond_matrix: ti.template(), particle_count: ti.i32):
    """Update bond metrics - module-level kernel"""
    bond_count = 0
    
    for i in range(particle_count):
        for j in range(i + 1, particle_count):
            if bond_matrix[i, j] > 0.0:
                bond_count += 1
    
    bond_count_field[None] = bond_count

@ti.kernel
def update_cluster_metrics_kernel(cluster_sizes: ti.template(), particle_count: ti.i32):
    """Update cluster metrics - module-level kernel"""
    cluster_count = 0
    
    for i in range(MAX_PARTICLES_COMPILE):
        if cluster_sizes[i] > 1:
            cluster_count += 1
    
    cluster_count_field[None] = cluster_count

@ti.kernel
def update_energy_field_metrics_kernel(energy_field: ti.template()):
    """Update energy field metrics - module-level kernel"""
    field_sum = 0.0
    field_max = 0.0
    
    for i, j in ti.ndrange(GRID_WIDTH_COMPILE, GRID_HEIGHT_COMPILE):
        value = energy_field[i, j]
        field_sum += value
        field_max = ti.max(field_max, value)
    
    energy_field_sum_field[None] = field_sum
    energy_field_max_field[None] = field_max
    energy_field_mean_field[None] = field_sum / (GRID_WIDTH_COMPILE * GRID_HEIGHT_COMPILE)

@ti.data_oriented
class MetricsCollector:
    """Collects and tracks simulation metrics"""
    
    def __init__(self, max_particles: int, history_size: int = 1000):
        self.max_particles = max_particles
        self.history_size = history_size
        
        # Initialize global fields if not already done
        if particle_count_field is None:
            init_metrics_fields()
        
        # Use global fields
        self.particle_count = particle_count_field
        self.total_energy = total_energy_field
        self.total_mass = total_mass_field
        self.bond_count = bond_count_field
        self.cluster_count = cluster_count_field
        self.energy_field_sum = energy_field_sum_field
        self.energy_field_max = energy_field_max_field
        self.energy_field_mean = energy_field_mean_field
        
        # Initialize
        self.reset()
        
        # History tracking
        self.metrics_history = deque(maxlen=history_size)
        self.start_time = time.time()
    
    def reset(self):
        """Reset all metrics"""
        reset_metrics_kernel()
    
    def update_particle_metrics(self, active, attributes, energy, particle_count: int):
        """Update particle-related metrics"""
        update_particle_metrics_kernel(active, attributes, energy, particle_count)
    
    def update_bond_metrics(self, bond_matrix, particle_count: int):
        """Update bond-related metrics"""
        update_bond_metrics_kernel(bond_matrix, particle_count)
    
    def update_cluster_metrics(self, cluster_sizes, particle_count: int):
        """Update cluster-related metrics"""
        update_cluster_metrics_kernel(cluster_sizes, particle_count)
    
    def update_energy_field_metrics(self, energy_field, width: int, height: int):
        """Update energy field metrics"""
        update_energy_field_metrics_kernel(energy_field)
    
    def record_metrics(self, additional_metrics: Dict = None):
        """Record current metrics to history"""
        current_time = time.time()
        
        metrics = {
            'timestamp': current_time,
            'runtime': current_time - self.start_time,
            'particle_count': self.particle_count[None],
            'total_energy': self.total_energy[None],
            'total_mass': self.total_mass[None],
            'bond_count': self.bond_count[None],
            'cluster_count': self.cluster_count[None],
            'energy_field_sum': self.energy_field_sum[None],
            'energy_field_max': self.energy_field_max[None],
            'energy_field_mean': self.energy_field_mean[None]
        }
        
        if additional_metrics:
            metrics.update(additional_metrics)
        
        self.metrics_history.append(metrics)
    
    def get_current_metrics(self) -> Dict:
        """Get current metrics snapshot"""
        return {
            'particle_count': int(self.particle_count[None]),
            'total_energy': float(self.total_energy[None]),
            'total_mass': float(self.total_mass[None]),
            'bond_count': int(self.bond_count[None]),
            'cluster_count': int(self.cluster_count[None]),
            'energy_field_sum': float(self.energy_field_sum[None]),
            'energy_field_max': float(self.energy_field_max[None]),
            'energy_field_mean': float(self.energy_field_mean[None])
        }
    
    def get_metrics_history(self) -> List[Dict]:
        """Get metrics history"""
        return list(self.metrics_history)
    
    def get_metrics_trends(self, window_size: int = 100) -> Dict:
        """Get trends in metrics over time"""
        if len(self.metrics_history) < 2:
            return {}
        
        recent_metrics = list(self.metrics_history)[-window_size:]
        
        trends = {}
        for key in ['particle_count', 'total_energy', 'bond_count', 'cluster_count']:
            values = [m[key] for m in recent_metrics if key in m]
            if len(values) > 1:
                trend = (values[-1] - values[0]) / len(values)
                trends[f'{key}_trend'] = float(trend)
        
        return trends

class NoveltyTracker:
    """Tracks novelty metrics for the simulation"""
    
    def __init__(self, window_size: int = 100):
        self.window_size = window_size
        self.discovery_timeline = deque(maxlen=window_size)
        self.novelty_rates = deque(maxlen=window_size)
        self.complexity_scores = deque(maxlen=window_size)
        
        self.start_time = time.time()
    
    def record_discovery(self, is_novel: bool, complexity: float = 0.0):
        """Record a substance discovery"""
        current_time = time.time()
        
        self.discovery_timeline.append((current_time, is_novel))
        self.complexity_scores.append(complexity)
        
        # Calculate current novelty rate
        recent_discoveries = list(self.discovery_timeline)[-self.window_size:]
        novel_count = sum(1 for _, is_novel in recent_discoveries if is_novel)
        total_count = len(recent_discoveries)
        
        novelty_rate = novel_count / max(total_count, 1)
        self.novelty_rates.append(novelty_rate)
    
    def get_novelty_rate(self) -> float:
        """Get current novelty rate"""
        if not self.novelty_rates:
            return 0.0
        return self.novelty_rates[-1]
    
    def get_average_complexity(self) -> float:
        """Get average complexity of recent discoveries"""
        if not self.complexity_scores:
            return 0.0
        return sum(self.complexity_scores) / len(self.complexity_scores)
    
    def get_discovery_rate(self) -> float:
        """Get discoveries per second"""
        if not self.discovery_timeline:
            return 0.0
        
        runtime = time.time() - self.start_time
        return len(self.discovery_timeline) / max(runtime, 1)
    
    def get_stats(self) -> Dict:
        """Get novelty tracking statistics"""
        return {
            'novelty_rate': float(self.get_novelty_rate()),
            'average_complexity': float(self.get_average_complexity()),
            'discovery_rate': float(self.get_discovery_rate()),
            'total_discoveries': int(len(self.discovery_timeline)),
            'novel_discoveries': int(sum(1 for _, is_novel in self.discovery_timeline if is_novel)),
            'runtime': float(time.time() - self.start_time)
        }

class ComplexityAnalyzer:
    """Analyzes complexity of molecular structures"""
    
    @staticmethod
    def calculate_graph_complexity(graph) -> float:
        """Calculate complexity score for a molecular graph"""
        if graph.num_nodes == 0:
            return 0.0
        
        # Base complexity from size
        size_complexity = np.log(graph.num_nodes + 1)
        
        # Connectivity complexity
        connectivity_complexity = graph.density * graph.num_nodes
        
        # Structural complexity (cycles, branches, etc.)
        structural_complexity = 0.0
        
        # Add complexity for cycles
        cycles = graph.get_cycles()
        structural_complexity += len(cycles) * 0.5
        
        # Add complexity for articulation points (branching)
        articulation_points = graph.get_articulation_points()
        structural_complexity += len(articulation_points) * 0.3
        
        # Add complexity for bridges (fragile connections)
        bridges = graph.get_bridges()
        structural_complexity += len(bridges) * 0.2
        
        # Total complexity
        total_complexity = size_complexity + connectivity_complexity + structural_complexity
        
        return total_complexity
    
    @staticmethod
    def calculate_chemical_complexity(particle_attributes: List[np.ndarray]) -> float:
        """Calculate chemical complexity based on particle diversity"""
        if not particle_attributes:
            return 0.0
        
        # Convert to numpy array for easier processing
        attributes = np.array(particle_attributes)
        
        # Mass diversity
        masses = attributes[:, 0]
        mass_diversity = np.std(masses) / (np.mean(masses) + 1e-6)
        
        # Charge diversity
        charges = attributes[:, 1:4]  # x, y, z components
        charge_diversity = np.mean([np.std(charges[:, i]) for i in range(3)])
        
        # Total chemical complexity
        chemical_complexity = mass_diversity + charge_diversity
        
        return chemical_complexity
    
    @staticmethod
    def calculate_emergent_complexity(graph, particle_attributes: List[np.ndarray]) -> float:
        """Calculate emergent complexity combining structural and chemical aspects"""
        structural_complexity = ComplexityAnalyzer.calculate_graph_complexity(graph)
        chemical_complexity = ComplexityAnalyzer.calculate_chemical_complexity(particle_attributes)
        
        # Emergent complexity is more than the sum of parts
        emergent_complexity = structural_complexity + chemical_complexity + \
                            0.1 * structural_complexity * chemical_complexity
        
        return emergent_complexity

class MetricsAggregator:
    """Aggregates metrics from multiple sources"""
    
    def __init__(self):
        self.metrics_collector = None
        self.novelty_tracker = None
        self.substance_catalog = None
        
        # Aggregated statistics
        self.aggregated_stats = {}
        self.last_update_time = time.time()
    
    def set_components(self, metrics_collector: MetricsCollector,
                      novelty_tracker: NoveltyTracker,
                      substance_catalog: SubstanceCatalog):
        """Set component references"""
        self.metrics_collector = metrics_collector
        self.novelty_tracker = novelty_tracker
        self.substance_catalog = substance_catalog
    
    def update_aggregated_stats(self):
        """Update aggregated statistics"""
        current_time = time.time()
        
        stats = {}
        
        # Basic metrics
        if self.metrics_collector:
            stats.update(self.metrics_collector.get_current_metrics())
            stats.update(self.metrics_collector.get_metrics_trends())
        
        # Novelty metrics
        if self.novelty_tracker:
            stats.update(self.novelty_tracker.get_stats())
        
        # Catalog metrics
        if self.substance_catalog:
            stats.update(self.substance_catalog.get_catalog_stats())
        
        # System metrics
        stats['update_frequency'] = 1.0 / max(current_time - self.last_update_time, 1e-6)
        stats['timestamp'] = current_time
        
        self.aggregated_stats = stats
        self.last_update_time = current_time
    
    def get_aggregated_stats(self) -> Dict:
        """Get aggregated statistics"""
        return self.aggregated_stats.copy()
    
    def get_health_score(self) -> float:
        """Calculate overall system health score"""
        if not self.aggregated_stats:
            return 0.0
        
        # Health factors
        novelty_rate = self.aggregated_stats.get('novelty_rate', 0.0)
        discovery_rate = self.aggregated_stats.get('discovery_rate', 0.0)
        particle_count = self.aggregated_stats.get('particle_count', 0)
        bond_count = self.aggregated_stats.get('bond_count', 0)
        
        # Ideal ranges
        ideal_novelty_rate = 0.1  # 10% novelty rate
        ideal_discovery_rate = 1.0  # 1 discovery per second
        ideal_particle_count = 1000  # Reasonable particle count
        ideal_bond_count = 500  # Reasonable bond count
        
        # Calculate health score (0-1)
        novelty_score = 1.0 - abs(novelty_rate - ideal_novelty_rate) / ideal_novelty_rate
        discovery_score = min(discovery_rate / ideal_discovery_rate, 1.0)
        particle_score = min(particle_count / ideal_particle_count, 1.0)
        bond_score = min(bond_count / ideal_bond_count, 1.0)
        
        health_score = (novelty_score + discovery_score + particle_score + bond_score) / 4.0
        
        return max(0.0, min(1.0, health_score))
