"""
Metrics system for Live 2.0 simulation
Tracks novelty, complexity, and other key metrics
"""

import taichi as ti
import numpy as np
import time
import logging
from typing import Dict, List, Tuple, Optional
from collections import deque
from .catalog import SubstanceCatalog

# Setup logging
logger = logging.getLogger(__name__)

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
                                  energy: ti.template(), velocities: ti.template(),
                                  num_particles: ti.i32):
    """Update particle metrics - module-level kernel"""
    particle_energy = 0.0
    total_mass = 0.0
    active_particles = 0

    # Only check up to the actual particle count, not the maximum
    max_check = ti.min(num_particles, MAX_PARTICLES_COMPILE)
    
    # Count all active particles
    for i in range(max_check):
        if active[i] == 1:
            active_particles += 1
            mass = attributes[i][0]
            total_mass += mass
            # Include stored per-particle energy (for preset systems)
            particle_energy += energy[i]
            # Add kinetic energy component (½ m v^2)
            vx = velocities[i][0]
            vy = velocities[i][1]
            particle_energy += 0.5 * mass * (vx * vx + vy * vy)

    particle_count_field[None] = active_particles
    # Store only particle+kinetic energy here, field energy added separately
    total_energy_field[None] = particle_energy
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
    
    def update_particle_metrics(self, active, attributes, energy, velocities, particle_count: int):
        """Update particle-related metrics - OPTIMIZED: Uses Taichi kernel only"""
        # Use GPU-accelerated Taichi kernel (100-1000x faster than Python loop!)
        update_particle_metrics_kernel(active, attributes, energy, velocities, particle_count)
        
        # Force synchronization to ensure kernel completes
        ti.sync()
        
        # Note: Taichi kernel updates particle_count_field, total_energy_field, and total_mass_field directly
        # total_energy_field contains particle+kinetic energy only
        # Energy field sum will be added when get_current_metrics is called
    
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
        particle_count = int(self.particle_count[None])
        particle_energy = float(self.total_energy[None])  # Only particle+kinetic energy
        total_mass = float(self.total_mass[None])
        bond_count = int(self.bond_count[None])
        cluster_count = int(self.cluster_count[None])
        energy_field_sum = float(self.energy_field_sum[None])
        energy_field_max = float(self.energy_field_max[None])
        energy_field_mean = float(self.energy_field_mean[None])
        
        # Total energy = particle/kinetic energy + energy field sum
        total_energy = particle_energy + energy_field_sum

        # logger.debug(f"get_current_metrics: particle_count={particle_count}")
        # logger.debug(f"get_current_metrics: total_energy={total_energy}")
        # logger.debug(f"get_current_metrics: total_mass={total_mass}")

        return {
            'particle_count': particle_count,
            'total_energy': total_energy,
            'total_mass': total_mass,
            'bond_count': bond_count,
            'cluster_count': cluster_count,
            'energy_field_sum': energy_field_sum,
            'energy_field_max': energy_field_max,
            'energy_field_mean': energy_field_mean
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
            basic_metrics = self.metrics_collector.get_current_metrics()
            # logger.debug(f"update_aggregated_stats: basic_metrics keys={list(basic_metrics.keys())}")
            # logger.debug(f"update_aggregated_stats: basic_metrics particle_count={basic_metrics.get('particle_count', 'NOT_FOUND')}")
            stats.update(basic_metrics)
            stats.update(self.metrics_collector.get_metrics_trends())

        # Novelty metrics
        if self.novelty_tracker:
            novelty_stats = self.novelty_tracker.get_stats()
            # logger.debug(f"update_aggregated_stats: novelty_stats keys={list(novelty_stats.keys())}")
            stats.update(novelty_stats)

        # Catalog metrics
        if self.substance_catalog:
            catalog_stats = self.substance_catalog.get_catalog_stats()
            # logger.debug(f"update_aggregated_stats: catalog_stats keys={list(catalog_stats.keys())}")
            stats.update(catalog_stats)

        # System metrics
        stats['update_frequency'] = 1.0 / max(current_time - self.last_update_time, 1e-6)
        stats['timestamp'] = current_time

        # Add health score
        health_score = self.get_health_score()
        stats['health_score'] = health_score
        # logger.debug(f"update_aggregated_stats: calculated health_score={health_score}")
        
        # logger.debug(f"update_aggregated_stats: final stats keys={list(stats.keys())}")
        # logger.debug(f"update_aggregated_stats: final particle_count={stats.get('particle_count', 'NOT_FOUND')}")

        self.aggregated_stats = stats
        self.last_update_time = current_time
    
    def get_aggregated_stats(self) -> Dict:
        """Get aggregated statistics"""
        stats = self.aggregated_stats.copy()
        # logger.debug(f"get_aggregated_stats: returning {len(stats)} metrics")
        # logger.debug(f"get_aggregated_stats: particle_count={stats.get('particle_count', 'NOT_FOUND')}")
        return stats
    
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


class PerformanceMonitor:
    """Monitors simulation performance metrics like FPS, step times, etc."""
    
    def __init__(self, history_size: int = 100):
        self.history_size = history_size
        
        # Performance tracking
        self.step_times = deque(maxlen=history_size)
        self.visualization_times = deque(maxlen=history_size)
        self.broadcast_times = deque(maxlen=history_size)
        self.total_times = deque(maxlen=history_size)
        
        # Timing state
        self.step_start_time = None
        self.visualization_start_time = None
        self.broadcast_start_time = None
        
        # Performance thresholds
        self.target_fps = 60.0
        self.min_fps = 10.0
        self.max_step_time = 1.0 / self.min_fps  # 100ms max per step
        self.max_visualization_time = 0.1  # 100ms max for visualization
        
        # Performance status
        self.performance_status = "unknown"
        self.last_warning_time = 0
        self.warning_cooldown = 5.0  # seconds between warnings
        
        # Cached metrics to avoid expensive calculations
        self._cached_metrics = {}
        self._last_metrics_update_step = 0
        self._metrics_update_interval = 50  # Update every 50 steps
        
    def start_step_timing(self):
        """Start timing a simulation step"""
        self.step_start_time = time.time()
        
    def end_step_timing(self):
        """End timing a simulation step"""
        if self.step_start_time is not None:
            step_time = time.time() - self.step_start_time
            self.step_times.append(step_time)
            self.step_start_time = None
            return step_time
        return 0.0
        
    def start_visualization_timing(self):
        """Start timing visualization data generation"""
        self.visualization_start_time = time.time()
        
    def end_visualization_timing(self):
        """End timing visualization data generation"""
        if self.visualization_start_time is not None:
            viz_time = time.time() - self.visualization_start_time
            self.visualization_times.append(viz_time)
            self.visualization_start_time = None
            return viz_time
        return 0.0
        
    def start_broadcast_timing(self):
        """Start timing WebSocket broadcast"""
        self.broadcast_start_time = time.time()
        
    def end_broadcast_timing(self):
        """End timing WebSocket broadcast"""
        if self.broadcast_start_time is not None:
            broadcast_time = time.time() - self.broadcast_start_time
            self.broadcast_times.append(broadcast_time)
            self.broadcast_start_time = None
            return broadcast_time
        return 0.0
        
    def get_performance_metrics(self, current_step: int = 0) -> Dict[str, float]:
        """Get current performance metrics - OPTIMIZED with caching"""
        # Check if we need to update cached metrics
        if (current_step - self._last_metrics_update_step >= self._metrics_update_interval or 
            not self._cached_metrics):
            self._update_cached_metrics()
            self._last_metrics_update_step = current_step
        
        return self._cached_metrics.copy()
        
    def _update_cached_metrics(self):
        """Update cached metrics - expensive operation done only when needed"""
        metrics = {}
        
        # Calculate FPS from step times
        if len(self.step_times) > 0:
            avg_step_time = np.mean(list(self.step_times))
            fps = 1.0 / avg_step_time if avg_step_time > 0 else 0.0
            metrics['fps'] = fps
            metrics['avg_step_time_ms'] = avg_step_time * 1000
            metrics['min_step_time_ms'] = min(self.step_times) * 1000
            metrics['max_step_time_ms'] = max(self.step_times) * 1000
        else:
            metrics['fps'] = 0.0
            metrics['avg_step_time_ms'] = 0.0
            metrics['min_step_time_ms'] = 0.0
            metrics['max_step_time_ms'] = 0.0
            
        # Visualization performance
        if len(self.visualization_times) > 0:
            metrics['avg_visualization_time_ms'] = np.mean(list(self.visualization_times)) * 1000
            metrics['max_visualization_time_ms'] = max(self.visualization_times) * 1000
        else:
            metrics['avg_visualization_time_ms'] = 0.0
            metrics['max_visualization_time_ms'] = 0.0
            
        # Broadcast performance
        if len(self.broadcast_times) > 0:
            metrics['avg_broadcast_time_ms'] = np.mean(list(self.broadcast_times)) * 1000
            metrics['max_broadcast_time_ms'] = max(self.broadcast_times) * 1000
        else:
            metrics['avg_broadcast_time_ms'] = 0.0
            metrics['max_broadcast_time_ms'] = 0.0
            
        # Performance status
        metrics['performance_status'] = self.get_performance_status()
        metrics['performance_score'] = self.get_performance_score()
        
        self._cached_metrics = metrics
        
    def get_performance_status(self) -> str:
        """Get current performance status"""
        if len(self.step_times) < 5:  # Need some samples
            return "warming_up"
            
        avg_fps = 1.0 / np.mean(list(self.step_times)) if np.mean(list(self.step_times)) > 0 else 0.0
        
        if avg_fps >= self.target_fps * 0.9:  # 90% of target
            return "excellent"
        elif avg_fps >= self.target_fps * 0.7:  # 70% of target
            return "good"
        elif avg_fps >= self.min_fps:
            return "acceptable"
        else:
            return "poor"
            
    def get_performance_score(self) -> float:
        """Get performance score (0.0 to 1.0)"""
        if len(self.step_times) < 5:
            return 0.5  # Unknown
            
        avg_fps = 1.0 / np.mean(list(self.step_times)) if np.mean(list(self.step_times)) > 0 else 0.0
        
        # Score based on FPS relative to target
        if avg_fps >= self.target_fps:
            return 1.0
        elif avg_fps >= self.min_fps:
            # Linear interpolation between min and target FPS
            return (avg_fps - self.min_fps) / (self.target_fps - self.min_fps)
        else:
            return 0.0
            
    def check_performance_warnings(self) -> List[str]:
        """Check for performance issues and return warnings"""
        warnings = []
        current_time = time.time()
        
        # Cooldown to avoid spam
        if current_time - self.last_warning_time < self.warning_cooldown:
            return warnings
            
        if len(self.step_times) < 5:
            return warnings
            
        avg_fps = 1.0 / np.mean(list(self.step_times)) if np.mean(list(self.step_times)) > 0 else 0.0
        avg_step_time = np.mean(list(self.step_times))
        
        # Check FPS
        if avg_fps < self.min_fps:
            warnings.append(f"Low FPS: {avg_fps:.1f} (target: {self.target_fps})")
            
        # Check step time
        if avg_step_time > self.max_step_time:
            warnings.append(f"Slow steps: {avg_step_time*1000:.1f}ms (max: {self.max_step_time*1000:.1f}ms)")
            
        # Check visualization time
        if len(self.visualization_times) > 0:
            avg_viz_time = np.mean(list(self.visualization_times))
            if avg_viz_time > self.max_visualization_time:
                warnings.append(f"Slow visualization: {avg_viz_time*1000:.1f}ms (max: {self.max_visualization_time*1000:.1f}ms)")
                
        if warnings:
            self.last_warning_time = current_time
            
        return warnings