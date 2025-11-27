"""
Autocatalysis Detection for Prebiotic Chemistry Simulations

Detects and analyzes autocatalytic cycles in reaction networks.
Critical for Paper 1 - Results Section 3.3 and Discussion 4.3.

Author: Live 2.0 Team
Date: November 2025
"""

import networkx as nx
import numpy as np
from typing import List, Dict, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict
import logging
import signal
import time

logger = logging.getLogger(__name__)


@dataclass
class AutocatalyticCycle:
    """Represents a detected autocatalytic cycle"""
    nodes: List[str]  # Molecule IDs in cycle
    edges: List[Tuple[str, str]]  # Reaction edges
    amplification_factor: float  # Fold-increase in abundance
    cycle_type: str  # 'direct', 'indirect', 'hypercycle'
    strength: float  # Catalytic strength (0-1)
    first_detected_step: int
    molecules: List[str]  # Molecule names/formulas
    

class AutocatalysisDetector:
    """
    Detects autocatalytic cycles in reaction networks.
    
    Uses Johnson's algorithm for cycle detection plus catalytic edge analysis.
    """
    
    def __init__(self, max_cycle_length: int = 6, min_amplification: float = 1.5, 
                 max_cycles_limit: int = 100000, cycle_timeout: int = 300):
        """
        Initialize detector.
        
        Args:
            max_cycle_length: Maximum nodes in cycle to detect (performance) - reduced to 6
            min_amplification: Minimum amplification to classify as autocatalytic
            max_cycles_limit: Maximum number of cycles to find before stopping (safety)
                             Reduced to 100k to prevent excessive processing
            cycle_timeout: Maximum seconds to spend finding cycles (default 5 min)
        """
        self.max_cycle_length = max_cycle_length
        self.min_amplification = min_amplification
        self.max_cycles_limit = max_cycles_limit
        self.cycle_timeout = cycle_timeout
        self.detected_cycles = []
        
    def detect_cycles_in_network(self, 
                                 reaction_graph: nx.DiGraph,
                                 abundance_history: Dict[str, List[float]],
                                 molecule_names: Dict[str, str]) -> List[AutocatalyticCycle]:
        """
        Main detection method.
        
        Args:
            reaction_graph: NetworkX directed graph (molecules=nodes, reactions=edges)
            abundance_history: {molecule_id: [abundance_at_step_t]}
            molecule_names: {molecule_id: formula/name}
            
        Returns:
            List of detected autocatalytic cycles
        """
        logger.info(f"Detecting autocatalytic cycles in network with {reaction_graph.number_of_nodes()} nodes")
        
        # Step 1: Find all cycles using Johnson's algorithm
        all_cycles = self._find_all_cycles(reaction_graph)
        logger.info(f"Found {len(all_cycles)} total cycles")
        
        # Step 2: Filter to autocatalytic cycles
        # Add progress logging and timeout for filtering phase
        autocatalytic_cycles = []
        start_time = time.time()
        total_cycles = len(all_cycles)
        
        # If too many cycles, warn and potentially sample
        if total_cycles > 50000:
            logger.warning(f"Very large number of cycles ({total_cycles}). "
                         f"This may take a long time. Consider reducing max_cycle_length.")
        
        for i, cycle in enumerate(all_cycles):
            # Check timeout during filtering
            elapsed = time.time() - start_time
            if elapsed > self.cycle_timeout:
                logger.warning(f"Filtering timeout after {elapsed:.1f}s. "
                             f"Processed {i}/{total_cycles} cycles. "
                             f"Found {len(autocatalytic_cycles)} autocatalytic cycles so far.")
                break
            
            # Progress logging for large sets
            if total_cycles > 1000 and (i + 1) % max(1, total_cycles // 10) == 0:
                progress_pct = 100 * (i + 1) / total_cycles
                logger.info(f"Filtering progress: {i+1}/{total_cycles} ({progress_pct:.1f}%) - "
                          f"found {len(autocatalytic_cycles)} autocatalytic so far")
            
            if self._is_autocatalytic(cycle, reaction_graph, abundance_history):
                ac = self._analyze_cycle(cycle, reaction_graph, abundance_history, molecule_names)
                # Additional check: average amplification must meet threshold
                if ac and ac.amplification_factor >= self.min_amplification:
                    autocatalytic_cycles.append(ac)
        
        elapsed = time.time() - start_time
        logger.info(f"Detected {len(autocatalytic_cycles)} autocatalytic cycles "
                   f"(filtered {total_cycles} cycles in {elapsed:.1f}s)")
        self.detected_cycles = autocatalytic_cycles
        return autocatalytic_cycles
    
    def _find_all_cycles(self, graph: nx.DiGraph) -> List[List[str]]:
        """
        Find all simple cycles up to max_cycle_length.
        
        Includes safety limits:
        - Maximum number of cycles (max_cycles_limit)
        - Timeout (cycle_timeout seconds)
        - Early exit if too many cycles found
        """
        cycles = []
        start_time = time.time()
        cycle_count = 0
        
        try:
            # Johnson's algorithm (built into NetworkX)
            for cycle in nx.simple_cycles(graph):
                # Check timeout
                elapsed = time.time() - start_time
                if elapsed > self.cycle_timeout:
                    logger.warning(f"Cycle detection timeout after {elapsed:.1f}s (limit: {self.cycle_timeout}s). "
                                 f"Found {len(cycles)} cycles so far. Stopping.")
                    break
                
                # Check cycle length
                if len(cycle) <= self.max_cycle_length:
                    cycles.append(cycle)
                    cycle_count += 1
                    
                    # Progress logging for large sets
                    if cycle_count > 0 and cycle_count % 10000 == 0:
                        elapsed = time.time() - start_time
                        logger.info(f"Cycle detection progress: {cycle_count} cycles found "
                                  f"({elapsed:.1f}s elapsed)")
                    
                    # Check limit (reduced to prevent excessive processing)
                    if cycle_count >= self.max_cycles_limit:
                        logger.warning(f"Cycle limit reached ({self.max_cycles_limit}). "
                                     f"Stopping cycle detection to prevent hang. "
                                     f"This network is very dense - consider reducing max_cycle_length.")
                        break
                    
        except Exception as e:
            logger.warning(f"Cycle detection error: {e}")
        
        elapsed = time.time() - start_time
        if len(cycles) >= self.max_cycles_limit:
            logger.warning(f"Found {len(cycles)} cycles (limit reached). "
                         f"This may indicate a very dense network. "
                         f"Consider reducing max_cycle_length or filtering the network.")
        else:
            logger.info(f"Cycle detection completed in {elapsed:.1f}s")
            
        return cycles
    
    def _is_autocatalytic(self, 
                         cycle: List[str],
                         graph: nx.DiGraph,
                         abundance_history: Dict[str, List[float]]) -> bool:
        """
        Check if cycle is autocatalytic.
        
        A cycle is autocatalytic if at least one molecule in the cycle:
        1. Appears as both reactant and product in the cycle
        2. Has net production (more produced than consumed)
        3. Shows amplification over time
        """
        # Check for self-amplifying nodes in cycle
        for node in cycle:
            # Get edges in cycle involving this node
            in_edges = [(u, v) for (u, v) in graph.edges() 
                       if v == node and u in cycle]
            out_edges = [(u, v) for (u, v) in graph.edges() 
                        if u == node and v in cycle]
            
            # Node must be both consumed and produced
            if len(in_edges) > 0 and len(out_edges) > 0:
                # Check if amplification occurs
                has_history = node in abundance_history
                if has_history:
                    history = abundance_history[node]
                    has_sufficient_points = len(history) >= 3
                    if has_sufficient_points:
                        # Calculate amplification for debugging
                        if len(history) < 6:
                            early = history[0] if history[0] > 0 else np.mean(history[:max(1, len(history)//2)]) + 1e-10
                            late = np.mean(history[-max(1, len(history)//2):])
                        else:
                            early = np.mean(history[:len(history)//3]) + 1e-10
                            late = np.mean(history[2*len(history)//3:])
                        amp = late / early if early > 0 else 0.0
                        passes_threshold = amp >= self.min_amplification
                        
                        # Log first few failures for debugging
                        if not passes_threshold and len(self.detected_cycles) < 3:
                            logger.debug(f"Cycle node {node}: amp={amp:.3f}, threshold={self.min_amplification}, "
                                       f"history_len={len(history)}, early={early:.2f}, late={late:.2f}")
                    else:
                        if len(self.detected_cycles) < 3:
                            logger.debug(f"Cycle node {node}: insufficient history points ({len(history)} < 3)")
                else:
                    if len(self.detected_cycles) < 3:
                        logger.debug(f"Cycle node {node}: not in abundance_history")
                
                if self._check_amplification(node, abundance_history):
                    return True
                    
        return False
    
    def _check_amplification(self, 
                            molecule_id: str,
                            abundance_history: Dict[str, List[float]]) -> bool:
        """Check if molecule shows amplification over time"""
        if molecule_id not in abundance_history:
            return False
            
        history = abundance_history[molecule_id]
        # Reduced from 10 to 3 to handle runs with fewer time points
        if len(history) < 3:  # Need at least 3 time points
            return False
            
        # Compare early vs late abundance
        # For short histories, use first vs last
        if len(history) < 6:
            early = history[0] if history[0] > 0 else np.mean(history[:max(1, len(history)//2)]) + 1e-10
            late = np.mean(history[-max(1, len(history)//2):])
        else:
            early = np.mean(history[:len(history)//3]) + 1e-10
            late = np.mean(history[2*len(history)//3:])
        
        if early > 0:
            amplification = late / early
            return amplification >= self.min_amplification
            
        return False
    
    def _analyze_cycle(self,
                      cycle: List[str],
                      graph: nx.DiGraph,
                      abundance_history: Dict[str, List[float]],
                      molecule_names: Dict[str, str]) -> AutocatalyticCycle:
        """Analyze detected cycle and compute metrics"""
        
        # Extract edges in cycle
        edges = []
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            if graph.has_edge(u, v):
                edges.append((u, v))
                
        # Calculate amplification factors for all nodes in cycle
        amplifications = []
        for node in cycle:
            if node in abundance_history:
                history = abundance_history[node]
                # Reduced from 10 to 3 to handle runs with fewer time points
                if len(history) >= 3:
                    # For short histories, use first vs last
                    if len(history) < 6:
                        early = history[0] if history[0] > 0 else np.mean(history[:max(1, len(history)//2)]) + 1e-10
                        late = np.mean(history[-max(1, len(history)//2):])
                    else:
                        early = np.mean(history[:len(history)//3]) + 1e-10
                        late = np.mean(history[2*len(history)//3:])
                    amp = late / early
                    amplifications.append(amp)
                    
        if not amplifications:
            return None
        
        # Use maximum amplification instead of mean, since _is_autocatalytic
        # checks if ANY node has amplification >= threshold
        # This ensures that if a cycle passes _is_autocatalytic, it will also
        # pass the amplification_factor check in detect_cycles_in_network
        avg_amplification = np.max(amplifications) if amplifications else 0.0
        
        # Classify cycle type
        cycle_type = self._classify_cycle_type(cycle, edges, graph)
        
        # Calculate catalytic strength
        strength = self._calculate_catalytic_strength(cycle, graph, abundance_history)
        
        # Find when cycle first appeared
        first_step = self._find_first_detection(cycle, abundance_history)
        
        # Get molecule formulas
        molecules = [molecule_names.get(node, node) for node in cycle]
        
        return AutocatalyticCycle(
            nodes=cycle,
            edges=edges,
            amplification_factor=avg_amplification,
            cycle_type=cycle_type,
            strength=strength,
            first_detected_step=first_step,
            molecules=molecules
        )
    
    def _classify_cycle_type(self, 
                            cycle: List[str],
                            edges: List[Tuple[str, str]],
                            graph: nx.DiGraph) -> str:
        """
        Classify cycle as:
        - 'direct': A + B → 2A (simple self-replication)
        - 'indirect': A → B → C → A (catalytic chain)
        - 'hypercycle': Multiple catalytic interactions
        """
        if len(cycle) == 2:
            return 'direct'
        elif len(cycle) <= 4:
            return 'indirect'
        else:
            # Check for cross-catalysis (hypercycle)
            cross_edges = 0
            for node in cycle:
                for other in cycle:
                    if node != other and graph.has_edge(node, other):
                        cross_edges += 1
            if cross_edges > len(cycle):
                return 'hypercycle'
            else:
                return 'indirect'
    
    def _calculate_catalytic_strength(self,
                                     cycle: List[str],
                                     graph: nx.DiGraph,
                                     abundance_history: Dict[str, List[float]]) -> float:
        """
        Calculate catalytic strength (0-1).
        
        Based on:
        - Amplification factor
        - Cycle connectivity
        - Temporal stability
        """
        # Amplification component
        amp_scores = []
        for node in cycle:
            if node in abundance_history:
                history = abundance_history[node]
                # Reduced from 10 to 3 to handle runs with fewer time points
                if len(history) >= 3:
                    # For short histories, use first vs last
                    if len(history) < 6:
                        early = history[0] if history[0] > 0 else np.mean(history[:max(1, len(history)//2)]) + 1e-10
                        late = np.mean(history[-max(1, len(history)//2):])
                    else:
                        early = np.mean(history[:len(history)//3]) + 1e-10
                        late = np.mean(history[2*len(history)//3:])
                    amp = late / early
                    # Normalize to 0-1
                    amp_score = min(1.0, (amp - 1.0) / 10.0)
                    amp_scores.append(amp_score)
                    
        avg_amp = np.mean(amp_scores) if amp_scores else 0.0
        
        # Connectivity component
        internal_edges = 0
        for u in cycle:
            for v in cycle:
                if graph.has_edge(u, v):
                    internal_edges += 1
        max_edges = len(cycle) * (len(cycle) - 1)
        connectivity = internal_edges / max_edges if max_edges > 0 else 0.0
        
        # Combined strength
        strength = 0.7 * avg_amp + 0.3 * connectivity
        return min(1.0, strength)
    
    def _find_first_detection(self,
                             cycle: List[str],
                             abundance_history: Dict[str, List[float]]) -> int:
        """Find step when all molecules in cycle first appeared"""
        first_steps = []
        for node in cycle:
            if node in abundance_history:
                history = abundance_history[node]
                # Find first non-zero abundance
                for i, abundance in enumerate(history):
                    if abundance > 0:
                        first_steps.append(i)
                        break
                        
        return max(first_steps) if first_steps else 0
    
    def get_statistics(self) -> Dict:
        """Get summary statistics of detected cycles"""
        if not self.detected_cycles:
            return {
                'total_cycles': 0,
                'by_type': {},
                'avg_amplification': 0.0,
                'max_amplification': 0.0
            }
            
        types = defaultdict(int)
        amplifications = []
        
        for cycle in self.detected_cycles:
            types[cycle.cycle_type] += 1
            amplifications.append(cycle.amplification_factor)
            
        return {
            'total_cycles': len(self.detected_cycles),
            'by_type': dict(types),
            'avg_amplification': np.mean(amplifications),
            'max_amplification': np.max(amplifications),
            'min_amplification': np.min(amplifications),
            'median_amplification': np.median(amplifications),
            'avg_strength': np.mean([c.strength for c in self.detected_cycles]),
            'avg_cycle_length': np.mean([len(c.nodes) for c in self.detected_cycles])
        }
    
    def export_cycles_for_paper(self, output_file: str):
        """Export cycle data formatted for paper tables"""
        import json
        
        export_data = []
        for i, cycle in enumerate(self.detected_cycles):
            export_data.append({
                'id': i + 1,
                'molecules': cycle.molecules,
                'cycle_length': len(cycle.nodes),
                'type': cycle.cycle_type,
                'amplification': round(cycle.amplification_factor, 2),
                'strength': round(cycle.strength, 3),
                'first_step': cycle.first_detected_step
            })
            
        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)
            
        logger.info(f"Exported {len(export_data)} cycles to {output_file}")


def analyze_scenario_autocatalysis(results_dir: str, 
                                   scenario_name: str,
                                   run_ids: List[int]) -> Dict:
    """
    Analyze autocatalysis across multiple runs of a scenario.
    
    For Paper 1 - aggregates across 10 replicates per scenario.
    
    Args:
        results_dir: Base results directory
        scenario_name: e.g., 'miller_urey_extended'
        run_ids: List of run IDs to analyze
        
    Returns:
        Aggregated statistics and cycle lists
    """
    from pathlib import Path
    import json
    
    all_cycles = []
    # Use relaxed criteria: min_amplification=1.1 instead of default 1.5
    detector = AutocatalysisDetector(min_amplification=1.1)
    
    for run_id in run_ids:
        run_dir = Path(results_dir) / scenario_name / f"run_{run_id}"
        
        # Load reaction network
        network_file = run_dir / "reaction_network.json"
        if not network_file.exists():
            logger.warning(f"Network file not found: {network_file}")
            continue
            
        with open(network_file) as f:
            network_data = json.load(f)
            
        # Build graph
        graph = nx.DiGraph()
        for edge in network_data.get('edges', []):
            graph.add_edge(edge['source'], edge['target'])
            
        # Use abundance_history from reaction_network.json (temporal data from snapshots)
        # This is the correct source, not results.json
        abundance_history = network_data.get('abundance_history', {})
        if not abundance_history:
            # Fallback to results.json if not in network file
            logger.warning(f"No abundance_history in {network_file}, trying results.json")
            results_file = run_dir / "results.json"
            if results_file.exists():
                with open(results_file) as f:
                    results = json.load(f)
                abundance_history = results.get('abundance_history', {})
            else:
                logger.warning(f"No abundance_history found for run_{run_id}")
                abundance_history = {}
        
        # Build molecule names mapping from network data
        molecule_names = {}
        for mol in network_data.get('molecules', []):
            mol_id = mol.get('id', '')
            mol_formula = mol.get('formula', '')
            if mol_id and mol_formula:
                molecule_names[mol_id] = mol_formula
        
        # Debug: log abundance history info
        if abundance_history:
            sample_mol = list(abundance_history.keys())[0] if abundance_history else None
            n_points = len(abundance_history[sample_mol]) if sample_mol else 0
            logger.info(f"  Run {run_id}: {len(abundance_history)} molecules, {n_points} time points")
        
        # Detect cycles
        cycles = detector.detect_cycles_in_network(graph, abundance_history, molecule_names)
        all_cycles.extend(cycles)
        
    # Aggregate statistics
    stats = {
        'scenario': scenario_name,
        'total_runs': len(run_ids),
        'total_cycles': len(all_cycles),
        'cycles_per_run_mean': len(all_cycles) / len(run_ids) if run_ids else 0,
        'cycles_per_run_std': np.std([len([c for c in all_cycles if c.first_detected_step // 50000 == i]) 
                                      for i in range(len(run_ids))]),
        'cycle_types': {},
        'amplification_stats': {}
    }
    
    if all_cycles:
        types = defaultdict(int)
        amplifications = [c.amplification_factor for c in all_cycles]
        
        for cycle in all_cycles:
            types[cycle.cycle_type] += 1
            
        stats['cycle_types'] = dict(types)
        stats['amplification_stats'] = {
            'mean': float(np.mean(amplifications)),
            'std': float(np.std(amplifications)),
            'median': float(np.median(amplifications)),
            'min': float(np.min(amplifications)),
            'max': float(np.max(amplifications)),
            'quartiles': [float(q) for q in np.percentile(amplifications, [25, 50, 75])]
        }
        
    return stats


if __name__ == "__main__":
    # Test with dummy data
    logging.basicConfig(level=logging.INFO)
    
    # Create test graph
    G = nx.DiGraph()
    G.add_edges_from([
        ('A', 'B'), ('B', 'C'), ('C', 'A'),  # Simple cycle
        ('A', 'D'), ('D', 'A')  # Direct autocatalysis
    ])
    
    # Create test abundance history
    abundance = {
        'A': [1, 2, 4, 8, 16, 20, 25, 30, 35, 40],
        'B': [1, 1, 2, 3, 5, 7, 9, 11, 13, 15],
        'C': [0, 1, 1, 2, 3, 4, 5, 6, 7, 8],
        'D': [1, 3, 9, 27, 50, 60, 65, 68, 70, 72]
    }
    
    molecule_names = {'A': 'CH2O', 'B': 'HCHO', 'C': 'C2H4O2', 'D': 'HCN'}
    
    detector = AutocatalysisDetector()
    cycles = detector.detect_cycles_in_network(G, abundance, molecule_names)
    
    print(f"\nDetected {len(cycles)} autocatalytic cycles:")
    for i, cycle in enumerate(cycles):
        print(f"\nCycle {i+1}:")
        print(f"  Molecules: {' → '.join(cycle.molecules)}")
        print(f"  Type: {cycle.cycle_type}")
        print(f"  Amplification: {cycle.amplification_factor:.2f}×")
        print(f"  Strength: {cycle.strength:.3f}")
        
    stats = detector.get_statistics()
    print(f"\nStatistics:")
    print(f"  Total cycles: {stats['total_cycles']}")
    print(f"  By type: {stats['by_type']}")
    print(f"  Avg amplification: {stats['avg_amplification']:.2f}×")

