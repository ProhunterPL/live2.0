"""
Complexity Metrics for Prebiotic Chemistry Simulations

Quantifies "emergent complexity" and proto-life characteristics.
Basic version for Paper 1 - Discussion Section 4.1.

Author: Live 2.0 Team
Date: November 2025
"""

import numpy as np
import networkx as nx
from typing import Dict, List, Tuple
from collections import Counter
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class ComplexityMetrics:
    """Container for complexity measurements"""
    # Diversity metrics
    shannon_entropy: float  # Species diversity
    species_richness: int  # Total unique species
    evenness: float  # Distribution uniformity (0-1)
    
    # Network metrics
    network_connectivity: float  # Average degree
    clustering_coefficient: float  # Local clustering
    path_length: float  # Average shortest path
    
    # Autocatalytic metrics
    autocatalytic_fraction: float  # Fraction in catalytic cycles
    amplification_score: float  # Average amplification strength
    
    # Novelty metrics
    novelty_score: float  # Fraction of novel molecules
    complexity_growth_rate: float  # dComplexity/dt
    
    # Self-organization metrics
    self_organization_index: float  # Combined measure (0-1)
    


class ComplexityAnalyzer:
    """
    Analyzes molecular complexity and proto-life characteristics.
    
    Provides basic metrics for Paper 1, foundation for full proto-life
    metrics in Paper 2.
    """
    
    def __init__(self):
        self.history = []
        
    def calculate_metrics(self,
                         species_counts: Dict[str, int],
                         reaction_network: nx.Graph,
                         autocatalytic_molecules: List[str],
                         novel_molecules: List[str],
                         timestep: int) -> ComplexityMetrics:
        """
        Calculate all complexity metrics at current timestep.
        
        Args:
            species_counts: {molecule_id: abundance}
            reaction_network: NetworkX graph of reactions
            autocatalytic_molecules: List of molecules in catalytic cycles
            novel_molecules: List of detected novel molecules
            timestep: Current simulation step
            
        Returns:
            ComplexityMetrics object
        """
        # Diversity metrics
        shannon = self._shannon_entropy(species_counts)
        richness = len(species_counts)
        evenness = shannon / np.log(richness) if richness > 1 else 0.0
        
        # Network metrics
        connectivity = self._network_connectivity(reaction_network)
        clustering = self._clustering_coefficient(reaction_network)
        path_length = self._average_path_length(reaction_network)
        
        # Autocatalytic metrics
        total_mols = sum(species_counts.values())
        autocatalytic_count = sum(species_counts.get(mol, 0) 
                                  for mol in autocatalytic_molecules)
        auto_fraction = autocatalytic_count / total_mols if total_mols > 0 else 0.0
        
        # Amplification (simplified for Paper 1)
        amplification = self._estimate_amplification(species_counts, autocatalytic_molecules)
        
        # Novelty metrics
        novelty = len(novel_molecules) / richness if richness > 0 else 0.0
        
        # Growth rate
        growth_rate = self._calculate_growth_rate(shannon, timestep)
        
        # Self-organization index (composite)
        self_org = self._self_organization_index(
            shannon, clustering, auto_fraction, novelty
        )
        
        metrics = ComplexityMetrics(
            shannon_entropy=shannon,
            species_richness=richness,
            evenness=evenness,
            network_connectivity=connectivity,
            clustering_coefficient=clustering,
            path_length=path_length,
            autocatalytic_fraction=auto_fraction,
            amplification_score=amplification,
            novelty_score=novelty,
            complexity_growth_rate=growth_rate,
            self_organization_index=self_org
        )
        
        self.history.append((timestep, metrics))
        return metrics
    
    def _shannon_entropy(self, species_counts: Dict[str, int]) -> float:
        """
        Calculate Shannon entropy: H = -Σ p_i log(p_i)
        
        Measures diversity - higher H = more diverse community
        """
        total = sum(species_counts.values())
        if total == 0:
            return 0.0
            
        probabilities = np.array([count / total for count in species_counts.values()])
        probabilities = probabilities[probabilities > 0]  # Remove zeros
        
        return -np.sum(probabilities * np.log(probabilities))
    
    def _network_connectivity(self, graph: nx.Graph) -> float:
        """Average degree (connections per node)"""
        if graph.number_of_nodes() == 0:
            return 0.0
        degrees = [d for n, d in graph.degree()]
        return np.mean(degrees) if degrees else 0.0
    
    def _clustering_coefficient(self, graph: nx.Graph) -> float:
        """
        Average clustering coefficient.
        
        Measures local network structure - high clustering indicates
        modular organization (important for proto-metabolism).
        """
        if graph.number_of_nodes() < 3:
            return 0.0
        try:
            return nx.average_clustering(graph.to_undirected())
        except:
            return 0.0
    
    def _average_path_length(self, graph: nx.Graph) -> float:
        """
        Average shortest path length.
        
        Measures network efficiency - shorter paths = better connected
        """
        if graph.number_of_nodes() < 2:
            return 0.0
            
        try:
            # Use largest connected component
            if nx.is_connected(graph.to_undirected()):
                return nx.average_shortest_path_length(graph.to_undirected())
            else:
                largest = max(nx.connected_components(graph.to_undirected()), 
                            key=len)
                subgraph = graph.subgraph(largest).to_undirected()
                return nx.average_shortest_path_length(subgraph)
        except:
            return 0.0
    
    def _estimate_amplification(self,
                                species_counts: Dict[str, int],
                                autocatalytic_molecules: List[str]) -> float:
        """
        Estimate average amplification of autocatalytic molecules.
        
        Simplified for Paper 1 - just ratio of autocatalytic to non-autocatalytic
        """
        if not autocatalytic_molecules:
            return 1.0
            
        auto_abundances = [species_counts.get(mol, 0) 
                          for mol in autocatalytic_molecules]
        non_auto = [count for mol, count in species_counts.items() 
                   if mol not in autocatalytic_molecules]
        
        if auto_abundances and non_auto:
            avg_auto = np.mean(auto_abundances)
            avg_non = np.mean(non_auto)
            return avg_auto / avg_non if avg_non > 0 else 1.0
        return 1.0
    
    def _calculate_growth_rate(self, current_shannon: float, timestep: int) -> float:
        """
        Calculate rate of complexity growth: dH/dt
        
        Uses recent history (last 10% of steps)
        """
        if len(self.history) < 2:
            return 0.0
            
        # Get recent history
        recent = self.history[-max(2, len(self.history)//10):]
        if len(recent) < 2:
            return 0.0
            
        # Linear regression on entropy vs time
        times = np.array([t for t, m in recent])
        entropies = np.array([m.shannon_entropy for t, m in recent])
        
        if len(times) > 1:
            # Simple slope calculation
            dt = times[-1] - times[0]
            dH = entropies[-1] - entropies[0]
            return dH / dt if dt > 0 else 0.0
        return 0.0
    
    def _self_organization_index(self,
                                 shannon: float,
                                 clustering: float,
                                 auto_fraction: float,
                                 novelty: float) -> float:
        """
        Composite "self-organization" index (0-1).
        
        Combines multiple indicators:
        - High diversity (Shannon)
        - High clustering (modular structure)
        - High autocatalytic fraction (positive feedback)
        - High novelty (innovation)
        
        For Paper 1 Discussion 4.1: "Emergent Complexity Without Guidance"
        """
        # Normalize Shannon (typical range 0-4)
        norm_shannon = min(1.0, shannon / 4.0)
        
        # Clustering already 0-1
        norm_clustering = clustering
        
        # Auto fraction already 0-1
        norm_auto = auto_fraction
        
        # Novelty already 0-1
        norm_novelty = novelty
        
        # Weighted combination
        weights = {
            'diversity': 0.3,
            'clustering': 0.2,
            'autocatalysis': 0.3,
            'novelty': 0.2
        }
        
        index = (weights['diversity'] * norm_shannon +
                weights['clustering'] * norm_clustering +
                weights['autocatalysis'] * norm_auto +
                weights['novelty'] * norm_novelty)
        
        return index
    
    def get_temporal_evolution(self) -> Dict:
        """
        Get time series of all metrics.
        
        For plotting evolution of complexity over simulation.
        """
        if not self.history:
            return {}
            
        evolution = {
            'timesteps': [],
            'shannon_entropy': [],
            'species_richness': [],
            'network_connectivity': [],
            'autocatalytic_fraction': [],
            'self_organization_index': []
        }
        
        for timestep, metrics in self.history:
            evolution['timesteps'].append(timestep)
            evolution['shannon_entropy'].append(metrics.shannon_entropy)
            evolution['species_richness'].append(metrics.species_richness)
            evolution['network_connectivity'].append(metrics.network_connectivity)
            evolution['autocatalytic_fraction'].append(metrics.autocatalytic_fraction)
            evolution['self_organization_index'].append(metrics.self_organization_index)
            
        return evolution
    
    def get_summary_statistics(self) -> Dict:
        """
        Get summary statistics for paper reporting.
        
        Returns final values and temporal characteristics.
        """
        if not self.history:
            return {}
            
        final_metrics = self.history[-1][1]
        
        # Temporal trends
        evolution = self.get_temporal_evolution()
        
        # Calculate phases (exploration, diversification, consolidation)
        phases = self._identify_phases(evolution)
        
        return {
            'final_values': {
                'shannon_entropy': final_metrics.shannon_entropy,
                'species_richness': final_metrics.species_richness,
                'evenness': final_metrics.evenness,
                'network_connectivity': final_metrics.network_connectivity,
                'clustering_coefficient': final_metrics.clustering_coefficient,
                'autocatalytic_fraction': final_metrics.autocatalytic_fraction,
                'novelty_score': final_metrics.novelty_score,
                'self_organization_index': final_metrics.self_organization_index
            },
            'temporal_characteristics': {
                'max_shannon': max(evolution['shannon_entropy']),
                'max_richness': max(evolution['species_richness']),
                'avg_growth_rate': np.mean([m.complexity_growth_rate for t, m in self.history]),
                'phases': phases
            }
        }
    
    def _identify_phases(self, evolution: Dict) -> Dict:
        """
        Identify temporal phases: exploration, diversification, consolidation.
        
        For Discussion 4.1: "three-phase pattern"
        """
        if not evolution['timesteps']:
            return {}
            
        times = np.array(evolution['timesteps'])
        shannon = np.array(evolution['shannon_entropy'])
        
        # Split into thirds
        n = len(times)
        third = n // 3
        
        exploration = {
            'steps': (times[0], times[third] if third < n else times[-1]),
            'avg_shannon': np.mean(shannon[:third]) if third > 0 else shannon[0],
            'growth_rate': (shannon[third] - shannon[0]) / (times[third] - times[0]) 
                          if third > 0 and times[third] > times[0] else 0.0
        }
        
        diversification = {
            'steps': (times[third], times[2*third] if 2*third < n else times[-1]),
            'avg_shannon': np.mean(shannon[third:2*third]) if 2*third > third else shannon[third],
            'growth_rate': (shannon[2*third] - shannon[third]) / (times[2*third] - times[third])
                          if 2*third > third and times[2*third] > times[third] else 0.0
        }
        
        consolidation = {
            'steps': (times[2*third] if 2*third < n else times[-1], times[-1]),
            'avg_shannon': np.mean(shannon[2*third:]) if 2*third < n else shannon[-1],
            'growth_rate': (shannon[-1] - shannon[2*third]) / (times[-1] - times[2*third])
                          if 2*third < n and times[-1] > times[2*third] else 0.0
        }
        
        return {
            'exploration': exploration,
            'diversification': diversification,
            'consolidation': consolidation
        }
    
    def export_for_paper(self, output_file: str):
        """Export metrics formatted for paper"""
        import json
        
        summary = self.get_summary_statistics()
        evolution = self.get_temporal_evolution()
        
        export_data = {
            'summary': summary,
            'evolution': evolution
        }
        
        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)
            
        logger.info(f"Exported complexity metrics to {output_file}")


def analyze_scenario_complexity(results_dir: str,
                                scenario_name: str,
                                run_ids: List[int]) -> Dict:
    """
    Analyze complexity metrics across multiple runs.
    
    For Paper 1 - aggregates across 10 replicates per scenario.
    
    Returns:
        Aggregated statistics for reporting in paper
    """
    from pathlib import Path
    import json
    
    all_final_metrics = []
    all_phases = []
    
    for run_id in run_ids:
        run_dir = Path(results_dir) / scenario_name / f"run_{run_id}"
        
        # Load complexity metrics (if previously computed)
        metrics_file = run_dir / "complexity_metrics.json"
        if not metrics_file.exists():
            logger.warning(f"Metrics file not found: {metrics_file}")
            continue
            
        with open(metrics_file) as f:
            data = json.load(f)
            
        all_final_metrics.append(data['summary']['final_values'])
        all_phases.append(data['summary']['temporal_characteristics']['phases'])
        
    if not all_final_metrics:
        return {}
        
    # Aggregate across runs
    aggregated = {
        'scenario': scenario_name,
        'n_runs': len(all_final_metrics),
        'final_metrics_mean': {},
        'final_metrics_std': {},
        'phase_characteristics': {}
    }
    
    # Aggregate final values
    for key in all_final_metrics[0].keys():
        values = [m[key] for m in all_final_metrics]
        aggregated['final_metrics_mean'][key] = float(np.mean(values))
        aggregated['final_metrics_std'][key] = float(np.std(values))
        
    # Aggregate phase characteristics
    for phase in ['exploration', 'diversification', 'consolidation']:
        growth_rates = [p[phase]['growth_rate'] for p in all_phases]
        aggregated['phase_characteristics'][phase] = {
            'avg_growth_rate': float(np.mean(growth_rates)),
            'std_growth_rate': float(np.std(growth_rates))
        }
        
    return aggregated


if __name__ == "__main__":
    # Test with dummy data
    logging.basicConfig(level=logging.INFO)
    
    # Create test data
    species_counts = {
        'A': 10,
        'B': 8,
        'C': 6,
        'D': 15,
        'E': 3
    }
    
    # Create test network
    G = nx.Graph()
    G.add_edges_from([
        ('A', 'B'), ('B', 'C'), ('C', 'A'),
        ('D', 'E'), ('E', 'A')
    ])
    
    autocatalytic = ['A', 'B']
    novel = ['E']
    
    analyzer = ComplexityAnalyzer()
    
    # Simulate temporal evolution
    for step in range(0, 100000, 10000):
        # Grow species counts
        species_counts = {k: int(v * (1 + 0.1 * np.random.rand())) 
                         for k, v in species_counts.items()}
        
        metrics = analyzer.calculate_metrics(
            species_counts, G, autocatalytic, novel, step
        )
        
        print(f"\nStep {step}:")
        print(f"  Shannon entropy: {metrics.shannon_entropy:.3f}")
        print(f"  Species richness: {metrics.species_richness}")
        print(f"  Self-organization: {metrics.self_organization_index:.3f}")
        
    summary = analyzer.get_summary_statistics()
    print(f"\n=== Summary Statistics ===")
    print(f"Final Shannon: {summary['final_values']['shannon_entropy']:.3f}")
    print(f"Final self-org: {summary['final_values']['self_organization_index']:.3f}")
    print(f"\nPhases:")
    for phase, data in summary['temporal_characteristics']['phases'].items():
        print(f"  {phase}: growth rate = {data['growth_rate']:.6f}")

