"""
Graph representation for Live 2.0 simulation
Handles molecular graphs, isomorphism detection, and graph hashing
"""

import taichi as ti
import numpy as np
import hashlib
from typing import List, Tuple, Dict, Set, Optional
from collections import defaultdict
import networkx as nx

class MolecularGraph:
    """Represents a molecular graph with particles as nodes and bonds as edges"""
    
    def __init__(self, particles: List[int], bonds: List[Tuple[int, int]], 
                 particle_attributes: Dict[int, np.ndarray]):
        self.particles = particles
        self.bonds = bonds
        self.particle_attributes = particle_attributes
        
        # Create NetworkX graph for analysis
        self.nx_graph = nx.Graph()
        self.nx_graph.add_nodes_from(particles)
        self.nx_graph.add_edges_from(bonds)
        
        # Compute graph properties
        self._compute_properties()
    
    def _compute_properties(self):
        """Compute graph properties"""
        self.num_nodes = len(self.particles)
        self.num_edges = len(self.bonds)
        self.density = self.num_edges / max(self.num_nodes * (self.num_nodes - 1) / 2, 1)
        
        # Compute degree sequence
        degrees = [self.nx_graph.degree(node) for node in self.particles]
        self.degree_sequence = sorted(degrees, reverse=True)
        
        # Compute clustering coefficient
        self.clustering_coeff = nx.average_clustering(self.nx_graph)
        
        # Compute diameter (if connected)
        if nx.is_connected(self.nx_graph):
            self.diameter = nx.diameter(self.nx_graph)
        else:
            self.diameter = -1  # Disconnected graph
    
    def get_canonical_form(self) -> str:
        """Get canonical form of the graph for isomorphism detection"""
        # Use NetworkX's canonical form
        try:
            canonical = nx.weisfeiler_lehman_graph_hash(self.nx_graph)
            return canonical
        except:
            # Fallback to simpler hash
            return self._simple_hash()
    
    def _simple_hash(self) -> str:
        """Simple hash based on degree sequence and edge count"""
        # Create a string representation
        degree_str = ','.join(map(str, self.degree_sequence))
        edge_str = ','.join(sorted([f"{min(b)}-{max(b)}" for b in self.bonds]))
        
        hash_input = f"{self.num_nodes}:{self.num_edges}:{degree_str}:{edge_str}"
        return hashlib.md5(hash_input.encode()).hexdigest()
    
    def is_isomorphic(self, other: 'MolecularGraph') -> bool:
        """Check if this graph is isomorphic to another"""
        try:
            return nx.is_isomorphic(self.nx_graph, other.nx_graph)
        except:
            # Fallback to canonical form comparison
            return self.get_canonical_form() == other.get_canonical_form()
    
    def get_subgraphs(self, min_size: int = 2) -> List['MolecularGraph']:
        """Get all connected subgraphs of minimum size"""
        subgraphs = []
        
        for component in nx.connected_components(self.nx_graph):
            if len(component) >= min_size:
                subgraph_nodes = list(component)
                subgraph_edges = [(u, v) for u, v in self.bonds 
                                 if u in component and v in component]
                subgraph_attrs = {node: self.particle_attributes[node] 
                                for node in subgraph_nodes}
                
                subgraph = MolecularGraph(subgraph_nodes, subgraph_edges, subgraph_attrs)
                subgraphs.append(subgraph)
        
        return subgraphs
    
    def get_node_count(self) -> int:
        """Get number of nodes in the graph"""
        return len(self.particles)
    
    def get_edge_count(self) -> int:
        """Get number of edges in the graph"""
        return len(self.bonds)
    
    def get_density(self) -> float:
        """Get graph density (edges / max possible edges)"""
        return self.density
    
    def get_complexity(self) -> float:
        """Calculate approximate complexity metric"""
        # Simple complexity based on size and branching
        num_nodes = self.num_nodes
        num_edges = self.num_edges
        
        if num_nodes <= 1:
            return 0.0
        
        # Basic complexity combining size and structure
        complexity = num_nodes * (1 + self.density * self.clustering_coeff)
        
        return complexity
    
    def get_cycles(self) -> List[List[int]]:
        """Get all cycles in the graph"""
        try:
            cycles = list(nx.simple_cycles(self.nx_graph.to_directed()))
            return cycles
        except:
            return []
    
    def get_bridges(self) -> List[Tuple[int, int]]:
        """Get all bridges (edges whose removal disconnects the graph)"""
        try:
            bridges = list(nx.bridges(self.nx_graph))
            return bridges
        except:
            return []
    
    def get_articulation_points(self) -> List[int]:
        """Get all articulation points (nodes whose removal disconnects the graph)"""
        try:
            articulation_points = list(nx.articulation_points(self.nx_graph))
            return articulation_points
        except:
            return []
    
    def to_dict(self) -> Dict:
        """Convert graph to dictionary representation"""
        return {
            'particles': self.particles,
            'bonds': self.bonds,
            'particle_attributes': {k: v.tolist() for k, v in self.particle_attributes.items()},
            'properties': {
                'num_nodes': self.num_nodes,
                'num_edges': self.num_edges,
                'density': self.density,
                'degree_sequence': self.degree_sequence,
                'clustering_coeff': self.clustering_coeff,
                'diameter': self.diameter
            }
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'MolecularGraph':
        """Create graph from dictionary representation"""
        particles = data['particles']
        bonds = data['bonds']
        particle_attributes = {k: np.array(v) for k, v in data['particle_attributes'].items()}
        
        return cls(particles, bonds, particle_attributes)

class GraphCatalog:
    """Catalog of known molecular graphs for novelty detection"""
    
    def __init__(self):
        self.graphs: Dict[str, MolecularGraph] = {}
        self.graph_counts: Dict[str, int] = {}
        self.graph_first_seen: Dict[str, float] = {}
        self.graph_last_seen: Dict[str, float] = {}
        
        # Statistics
        self.total_graphs_seen = 0
        self.novel_graphs_count = 0
    
    def add_graph(self, graph: MolecularGraph, timestamp: float) -> bool:
        """Add a graph to the catalog. Returns True if novel."""
        canonical_form = graph.get_canonical_form()
        
        is_novel = canonical_form not in self.graphs
        
        if is_novel:
            self.graphs[canonical_form] = graph
            self.graph_counts[canonical_form] = 1
            self.graph_first_seen[canonical_form] = timestamp
            self.graph_last_seen[canonical_form] = timestamp
            self.novel_graphs_count += 1
        else:
            self.graph_counts[canonical_form] += 1
            self.graph_last_seen[canonical_form] = timestamp
        
        self.total_graphs_seen += 1
        return is_novel
    
    def get_novelty_rate(self, window_size: int = None) -> float:
        """Get novelty rate over recent window"""
        if window_size is None:
            window_size = 500  # Default from config
        if self.total_graphs_seen < window_size:
            return self.novel_graphs_count / max(self.total_graphs_seen, 1)
        
        # Count novel graphs in recent window
        recent_novel = 0
        recent_total = 0
        
        # This is a simplified calculation - in practice, you'd track timestamps
        # and count novel graphs in the last window_size additions
        for canonical_form, count in self.graph_counts.items():
            if count <= window_size:  # Approximate recent graphs
                recent_novel += 1
            recent_total += min(count, window_size)
        
        return recent_novel / max(recent_total, 1)
    
    def get_graph_stats(self) -> Dict:
        """Get catalog statistics"""
        if not self.graphs:
            return {
                'total_graphs': 0,
                'novel_graphs': 0,
                'novelty_rate': 0.0,
                'average_graph_size': 0.0,
                'most_common_size': 0
            }
        
        graph_sizes = [graph.num_nodes for graph in self.graphs.values()]
        size_counts = defaultdict(int)
        for size in graph_sizes:
            size_counts[size] += 1
        
        most_common_size = max(size_counts.items(), key=lambda x: x[1])[0]
        
        return {
            'total_graphs': len(self.graphs),
            'novel_graphs': self.novel_graphs_count,
            'novelty_rate': self.get_novelty_rate(),
            'average_graph_size': sum(graph_sizes) / len(graph_sizes),
            'most_common_size': most_common_size,
            'size_distribution': dict(size_counts)
        }
    
    def get_recent_graphs(self, count: int = 10) -> List[MolecularGraph]:
        """Get most recently seen graphs"""
        # Sort by last seen time (simplified - using count as proxy)
        sorted_graphs = sorted(self.graphs.items(), 
                            key=lambda x: self.graph_counts[x[0]], 
                            reverse=True)
        
        return [graph for _, graph in sorted_graphs[:count]]
    
    def clear(self):
        """Clear the catalog"""
        self.graphs.clear()
        self.graph_counts.clear()
        self.graph_first_seen.clear()
        self.graph_last_seen.clear()
        self.total_graphs_seen = 0
        self.novel_graphs_count = 0

@ti.data_oriented
class GraphProcessor:
    """Taichi-based graph processing for real-time analysis"""
    
    def __init__(self, max_particles: int):
        self.max_particles = max_particles
        
        # Graph data structures
        self.adjacency_matrix = ti.field(dtype=ti.i32, 
                                       shape=(max_particles, max_particles))
        self.degree_sequence = ti.field(dtype=ti.i32, shape=(max_particles,))
        self.graph_hash = ti.field(dtype=ti.u32, shape=())
        
        # Initialize
        self.reset()
    
    @ti.kernel
    def reset(self):
        """Reset graph processor"""
        for i, j in ti.ndrange(self.max_particles, self.max_particles):
            self.adjacency_matrix[i, j] = 0
        
        for i in range(self.max_particles):
            self.degree_sequence[i] = 0
        
        self.graph_hash[None] = 0
    
    @ti.kernel
    def update_adjacency_matrix(self, bond_matrix: ti.template(), 
                              active: ti.template(), particle_count: ti.i32):
        """Update adjacency matrix from bond matrix"""
        # Clear adjacency matrix
        for i, j in ti.ndrange(self.max_particles, self.max_particles):
            self.adjacency_matrix[i, j] = 0
        
        # Fill adjacency matrix
        for i in range(particle_count):
            if active[i] == 1:
                for j in range(i + 1, particle_count):
                    if active[j] == 1 and bond_matrix[i, j] > 0:
                        self.adjacency_matrix[i, j] = 1
                        self.adjacency_matrix[j, i] = 1
    
    @ti.kernel
    def compute_degrees(self, active: ti.template(), particle_count: ti.i32):
        """Compute degree sequence"""
        for i in range(self.max_particles):
            self.degree_sequence[i] = 0
        
        for i in range(particle_count):
            if active[i] == 1:
                degree = 0
                for j in range(particle_count):
                    if active[j] == 1 and self.adjacency_matrix[i, j] == 1:
                        degree += 1
                self.degree_sequence[i] = degree
    
    @ti.kernel
    def compute_graph_hash(self, active: ti.template(), particle_count: ti.i32):
        """Compute simple graph hash"""
        hash_val = 0
        
        # Hash based on degree sequence and edge count
        edge_count = 0
        for i in range(particle_count):
            if active[i] == 1:
                hash_val = hash_val * 31 + self.degree_sequence[i]
                for j in range(i + 1, particle_count):
                    if active[j] == 1 and self.adjacency_matrix[i, j] == 1:
                        edge_count += 1
        
        hash_val = hash_val * 31 + edge_count
        self.graph_hash[None] = hash_val
    
    def get_adjacency_matrix(self) -> np.ndarray:
        """Get adjacency matrix as numpy array"""
        return self.adjacency_matrix.to_numpy()
    
    def get_degree_sequence(self) -> np.ndarray:
        """Get degree sequence as numpy array"""
        return self.degree_sequence.to_numpy()
    
    def get_graph_hash(self) -> int:
        """Get computed graph hash"""
        return self.graph_hash[None]
