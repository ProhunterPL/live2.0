"""
Diagnostics and observables system for Live 2.0 simulation
Logs time series data to CSV for detailed analysis
"""

import csv
import time
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any
from collections import defaultdict, deque
import taichi as ti

class DiagnosticsLogger:
    """
    Comprehensive diagnostics logger for simulation observables
    Tracks bonds, clusters, events, and energy metrics over time
    """
    
    def __init__(self, output_dir: str = "diagnostics", enabled: bool = True):
        """
        Initialize diagnostics logger
        
        Args:
            output_dir: Directory to save CSV files
            enabled: Whether logging is enabled
        """
        self.enabled = enabled
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Time series data
        self.time_series: Dict[str, List] = defaultdict(list)
        
        # Event counters per step
        self.events_formed = 0
        self.events_broken = 0
        self.events_merged = 0
        self.events_split = 0
        
        # Bond lifetime tracking
        self.bond_lifetimes: Dict[tuple, float] = {}  # (i, j) -> birth_time
        self.bond_lifetime_history: List[float] = []
        
        # Previous state for event detection
        self.prev_bonds: set = set()
        self.prev_clusters: Dict[int, set] = {}
        
        # CSV files
        self.csv_files: Dict[str, Any] = {}
        self.csv_writers: Dict[str, Any] = {}
        
        # Session start time
        self.session_start = time.time()
        self.current_sim_time = 0.0
        self.step_count = 0
        
        # Initialize CSV files
        if self.enabled:
            self._init_csv_files()
    
    def _init_csv_files(self):
        """Initialize CSV files with headers"""
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        
        # Main metrics CSV
        metrics_path = self.output_dir / f"metrics_{timestamp}.csv"
        self.csv_files['metrics'] = open(metrics_path, 'w', newline='')
        self.csv_writers['metrics'] = csv.writer(self.csv_files['metrics'])
        self.csv_writers['metrics'].writerow([
            'step', 'sim_time', 'wall_time',
            'num_particles', 'num_bonds_total',
            'avg_bond_length', 'avg_bond_tension',
            'num_clusters', 'largest_cluster_size',
            'cluster_energy_mean', 'cluster_energy_var',
            'R_g_mean',
            'events_formed', 'events_broken',
            'events_merged', 'events_split'
        ])
        
        # Bond types CSV
        bond_types_path = self.output_dir / f"bond_types_{timestamp}.csv"
        self.csv_files['bond_types'] = open(bond_types_path, 'w', newline='')
        self.csv_writers['bond_types'] = csv.writer(self.csv_files['bond_types'])
        self.csv_writers['bond_types'].writerow([
            'step', 'sim_time', 'bond_type', 'count'
        ])
        
        # Cluster distribution CSV
        cluster_dist_path = self.output_dir / f"cluster_dist_{timestamp}.csv"
        self.csv_files['cluster_dist'] = open(cluster_dist_path, 'w', newline='')
        self.csv_writers['cluster_dist'] = csv.writer(self.csv_files['cluster_dist'])
        self.csv_writers['cluster_dist'].writerow([
            'step', 'sim_time', 'cluster_size', 'count'
        ])
        
        # Bond lifetime histogram CSV
        lifetime_path = self.output_dir / f"bond_lifetimes_{timestamp}.csv"
        self.csv_files['lifetimes'] = open(lifetime_path, 'w', newline='')
        self.csv_writers['lifetimes'] = csv.writer(self.csv_files['lifetimes'])
        self.csv_writers['lifetimes'].writerow([
            'step', 'sim_time', 'lifetime', 'count'
        ])
    
    def log_step(self, step: int, sim_time: float, simulation_data: Dict[str, Any]):
        """
        Log metrics for current simulation step
        
        Args:
            step: Current step number
            sim_time: Current simulation time
            simulation_data: Dictionary containing all simulation state
        """
        if not self.enabled:
            return
        
        self.step_count = step
        self.current_sim_time = sim_time
        wall_time = time.time() - self.session_start
        
        # Extract data
        positions = simulation_data.get('positions', np.array([]))
        bonds = simulation_data.get('bonds', [])
        clusters = simulation_data.get('clusters', {})
        particle_energies = simulation_data.get('particle_energies', np.array([]))
        attributes = simulation_data.get('attributes', np.array([]))
        
        # Compute bond metrics
        bond_metrics = self._compute_bond_metrics(positions, bonds, attributes)
        
        # Compute cluster metrics
        cluster_metrics = self._compute_cluster_metrics(positions, clusters, particle_energies)
        
        # Detect events
        events = self._detect_events(bonds, clusters, sim_time)
        
        # Write main metrics
        self.csv_writers['metrics'].writerow([
            step, sim_time, wall_time,
            len(positions),
            bond_metrics['num_bonds_total'],
            bond_metrics['avg_bond_length'],
            bond_metrics['avg_bond_tension'],
            cluster_metrics['num_clusters'],
            cluster_metrics['largest_cluster_size'],
            cluster_metrics['cluster_energy_mean'],
            cluster_metrics['cluster_energy_var'],
            cluster_metrics['R_g_mean'],
            events['formed'], events['broken'],
            events['merged'], events['split']
        ])
        
        # Write bond types
        for bond_type, count in bond_metrics['bonds_by_type'].items():
            self.csv_writers['bond_types'].writerow([
                step, sim_time, bond_type, count
            ])
        
        # Write cluster size distribution (top-K)
        for size, count in sorted(cluster_metrics['size_histogram'].items(), 
                                  key=lambda x: x[1], reverse=True)[:10]:
            self.csv_writers['cluster_dist'].writerow([
                step, sim_time, size, count
            ])
        
        # Write bond lifetime histogram (periodically)
        if step % 100 == 0 and self.bond_lifetime_history:
            lifetime_hist = np.histogram(self.bond_lifetime_history, bins=20)
            for lifetime_bin, count in zip(lifetime_hist[1][:-1], lifetime_hist[0]):
                if count > 0:
                    self.csv_writers['lifetimes'].writerow([
                        step, sim_time, lifetime_bin, count
                    ])
        
        # Flush periodically
        if step % 10 == 0:
            for f in self.csv_files.values():
                f.flush()
    
    def _compute_bond_metrics(self, positions: np.ndarray, bonds: List[tuple], 
                              attributes: np.ndarray) -> Dict[str, Any]:
        """Compute bond-related metrics"""
        metrics = {
            'num_bonds_total': len(bonds),
            'bonds_by_type': defaultdict(int),
            'avg_bond_length': 0.0,
            'avg_bond_tension': 0.0
        }
        
        if len(bonds) == 0:
            return metrics
        
        bond_lengths = []
        bond_tensions = []
        
        for i, j, strength in bonds:
            # Bond length
            if i < len(positions) and j < len(positions):
                r_vec = positions[i] - positions[j]
                length = np.linalg.norm(r_vec)
                bond_lengths.append(length)
                
                # Bond tension (deviation from equilibrium)
                equilibrium_length = 1.0  # Default
                tension = abs(length - equilibrium_length)
                bond_tensions.append(tension)
                
                # Bond type classification
                if i < len(attributes) and j < len(attributes):
                    mass_i = attributes[i][0] if len(attributes[i]) > 0 else 1.0
                    mass_j = attributes[j][0] if len(attributes[j]) > 0 else 1.0
                    
                    # Classify by mass ratio
                    mass_ratio = max(mass_i, mass_j) / (min(mass_i, mass_j) + 1e-6)
                    if mass_ratio < 1.5:
                        bond_type = 'homogeneous'
                    elif mass_ratio < 3.0:
                        bond_type = 'weakly_heterogeneous'
                    else:
                        bond_type = 'strongly_heterogeneous'
                    
                    metrics['bonds_by_type'][bond_type] += 1
        
        if bond_lengths:
            metrics['avg_bond_length'] = np.mean(bond_lengths)
        if bond_tensions:
            metrics['avg_bond_tension'] = np.mean(bond_tensions)
        
        return metrics
    
    def _compute_cluster_metrics(self, positions: np.ndarray, 
                                 clusters: Dict[int, set],
                                 particle_energies: np.ndarray) -> Dict[str, Any]:
        """Compute cluster-related metrics"""
        metrics = {
            'num_clusters': len(clusters),
            'size_histogram': defaultdict(int),
            'largest_cluster_size': 0,
            'cluster_energy_mean': 0.0,
            'cluster_energy_var': 0.0,
            'R_g_mean': 0.0
        }
        
        if not clusters:
            return metrics
        
        cluster_sizes = []
        cluster_energies = []
        cluster_Rgs = []
        
        for cluster_id, particle_ids in clusters.items():
            size = len(particle_ids)
            cluster_sizes.append(size)
            metrics['size_histogram'][size] += 1
            
            if size > metrics['largest_cluster_size']:
                metrics['largest_cluster_size'] = size
            
            # Compute cluster energy
            if len(particle_energies) > 0:
                particle_list = list(particle_ids)
                valid_particles = [p for p in particle_list if p < len(particle_energies)]
                if valid_particles:
                    cluster_energy = np.sum(particle_energies[valid_particles])
                    cluster_energies.append(cluster_energy)
            
            # Compute radius of gyration (R_g)
            if len(positions) > 0:
                particle_list = list(particle_ids)
                valid_particles = [p for p in particle_list if p < len(positions)]
                if len(valid_particles) >= 2:
                    cluster_positions = positions[valid_particles]
                    center_of_mass = np.mean(cluster_positions, axis=0)
                    R_g_squared = np.mean(np.sum((cluster_positions - center_of_mass)**2, axis=1))
                    R_g = np.sqrt(R_g_squared)
                    cluster_Rgs.append(R_g)
        
        if cluster_energies:
            metrics['cluster_energy_mean'] = np.mean(cluster_energies)
            metrics['cluster_energy_var'] = np.var(cluster_energies)
        
        if cluster_Rgs:
            metrics['R_g_mean'] = np.mean(cluster_Rgs)
        
        return metrics
    
    def _detect_events(self, bonds: List[tuple], clusters: Dict[int, set], 
                      sim_time: float) -> Dict[str, int]:
        """Detect events: bond formation/breaking, cluster merge/split"""
        events = {
            'formed': 0,
            'broken': 0,
            'merged': 0,
            'split': 0
        }
        
        # Convert bonds to set for comparison
        current_bonds = {(min(i, j), max(i, j)) for i, j, _ in bonds}
        
        # Detect formed bonds
        formed_bonds = current_bonds - self.prev_bonds
        events['formed'] = len(formed_bonds)
        
        # Track bond lifetimes for formed bonds
        for bond in formed_bonds:
            self.bond_lifetimes[bond] = sim_time
        
        # Detect broken bonds
        broken_bonds = self.prev_bonds - current_bonds
        events['broken'] = len(broken_bonds)
        
        # Record bond lifetimes for broken bonds
        for bond in broken_bonds:
            if bond in self.bond_lifetimes:
                lifetime = sim_time - self.bond_lifetimes[bond]
                self.bond_lifetime_history.append(lifetime)
                del self.bond_lifetimes[bond]
        
        # Detect cluster merge/split events
        if self.prev_clusters:
            # Simple heuristic: if number of clusters decreased -> merge
            # if increased -> split
            prev_count = len(self.prev_clusters)
            curr_count = len(clusters)
            
            if curr_count < prev_count:
                events['merged'] = prev_count - curr_count
            elif curr_count > prev_count:
                events['split'] = curr_count - prev_count
        
        # Update previous state
        self.prev_bonds = current_bonds
        self.prev_clusters = {k: set(v) for k, v in clusters.items()}
        
        return events
    
    def get_summary_statistics(self) -> Dict[str, Any]:
        """Get summary statistics for current session"""
        if not self.bond_lifetime_history:
            mean_lifetime = 0.0
            median_lifetime = 0.0
        else:
            mean_lifetime = np.mean(self.bond_lifetime_history)
            median_lifetime = np.median(self.bond_lifetime_history)
        
        return {
            'total_steps': self.step_count,
            'total_sim_time': self.current_sim_time,
            'total_bonds_formed': sum(1 for _ in self.bond_lifetimes.keys()) + len(self.bond_lifetime_history),
            'total_bonds_broken': len(self.bond_lifetime_history),
            'mean_bond_lifetime': mean_lifetime,
            'median_bond_lifetime': median_lifetime,
            'active_bonds': len(self.bond_lifetimes)
        }
    
    def close(self):
        """Close all CSV files"""
        if not self.enabled:
            return
        
        # Write summary
        summary = self.get_summary_statistics()
        summary_path = self.output_dir / "session_summary.txt"
        with open(summary_path, 'w') as f:
            f.write("Diagnostics Session Summary\n")
            f.write("=" * 50 + "\n\n")
            for key, value in summary.items():
                f.write(f"{key}: {value}\n")
        
        # Close all CSV files
        for f in self.csv_files.values():
            f.close()
        
        print(f"Diagnostics saved to {self.output_dir}")
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

