"""
Autocatalysis Network Exporter for dataset export.

Exports autocatalysis cycles as network graph.
Uses existing AutocatalyticDetector to find cycles, then formats
as NetworkX graph for export (JSON or GraphML).
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional

from backend.dataset_export.formatters import FormatterFactory
from backend.dataset_export.utils import resolve_run_patterns

logger = logging.getLogger(__name__)

# Try to import NetworkX
try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    logger.warning("networkx not available, autocatalysis network export disabled")


class AutocatalysisNetworkExporter:
    """
    Exports autocatalysis cycles as network graph.
    
    Uses existing AutocatalyticDetector to find cycles, then formats
    as NetworkX graph for export (JSON or GraphML).
    """
    
    def __init__(self, base_results_dir: str):
        """
        Initialize exporter.
        
        Args:
            base_results_dir: Base directory with results
        """
        self.base_results_dir = Path(base_results_dir)
    
    def export(
        self,
        runs: List[str],
        output_format: str = "json",
        output_path: str = None,
        include_metrics: bool = True
    ) -> str:
        """
        Export autocatalysis network.
        
        Process:
        1. For each run, check if autocatalytic_cycles.json exists
        2. If not, check if reaction_network.json exists
        3. If not, build reaction_network.json from snapshots
        4. Run AutocatalyticDetector to find cycles
        5. Build NetworkX graph from cycles
        6. Export as JSON or GraphML
        
        Args:
            runs: List of run patterns
            output_format: "json" or "graphml"
            output_path: Path to output file
            include_metrics: Whether to include metrics in graph
        
        Returns:
            Path to exported file
        """
        if not NETWORKX_AVAILABLE:
            raise ImportError("networkx required for autocatalysis network export")
        
        if not output_path:
            raise ValueError("output_path is required")
        
        # Resolve runs
        run_dirs = resolve_run_patterns(runs, self.base_results_dir)
        
        if not run_dirs:
            logger.warning("No run directories found")
            return output_path
        
        # Collect all cycles
        all_cycles = []
        
        for run_dir in run_dirs:
            try:
                cycles = self._get_cycles_from_run(run_dir)
                all_cycles.extend(cycles)
            except Exception as e:
                logger.error(f"Failed to get cycles from {run_dir.name}: {e}")
                continue
        
        if not all_cycles:
            logger.warning("No autocatalytic cycles found")
            # Create empty graph
            graph = nx.DiGraph()
        else:
            # Build network graph
            graph = self._build_network_graph(all_cycles, include_metrics)
        
        # Export
        formatter = FormatterFactory.create(output_format)
        formatter.export_graph(graph, output_path)
        
        logger.info(f"Exported autocatalysis network: {output_path} ({graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges)")
        return output_path
    
    def _get_cycles_from_run(self, run_dir: Path) -> List[Dict]:
        """
        Get autocatalytic cycles from a run.
        
        Checks in order:
        1. autocatalytic_cycles.json (if exists)
        2. reaction_network.json + AutocatalyticDetector
        3. Build reaction_network.json from snapshots first
        
        Args:
            run_dir: Path to run directory
        
        Returns:
            List of cycle dicts
        """
        # Check for pre-computed cycles
        cycles_file = run_dir / "autocatalytic_cycles.json"
        if cycles_file.exists():
            try:
                with open(cycles_file) as f:
                    data = json.load(f)
                cycles = data.get("cycles", [])
                if cycles:
                    logger.debug(f"Loaded {len(cycles)} cycles from {cycles_file.name}")
                    return cycles
            except Exception as e:
                logger.warning(f"Failed to load {cycles_file.name}: {e}")
        
        # Check for reaction network
        network_file = run_dir / "reaction_network.json"
        if not network_file.exists():
            # Build from snapshots
            logger.info(f"Building reaction network for {run_dir.name}")
            try:
                from scripts.build_reaction_network_from_snapshots import build_reaction_network
                network = build_reaction_network(run_dir)
                
                if network:
                    with open(network_file, 'w') as f:
                        json.dump(network, f, indent=2)
                    logger.info(f"Built reaction network for {run_dir.name}")
                else:
                    logger.warning(f"Failed to build network for {run_dir.name}")
                    return []
            except Exception as e:
                logger.error(f"Failed to build reaction network for {run_dir.name}: {e}")
                return []
        
        # Run autocatalysis detection
        try:
            from scripts.autocatalytic_detector import AutocatalyticDetector
            
            analysis_dir = run_dir / "analysis"
            analysis_dir.mkdir(parents=True, exist_ok=True)
            
            detector = AutocatalyticDetector(network_file, analysis_dir)
            detector.load_network()
            
            # Detect cycles
            direct_cycles = detector.detect_direct_autocatalysis()
            indirect_cycles = detector.detect_indirect_autocatalysis(max_length=10)
            
            # Convert to dict format
            cycles = []
            for cycle in direct_cycles + indirect_cycles:
                cycles.append({
                    "molecules": cycle.molecules,
                    "reactions": cycle.reactions,
                    "amplification_factor": cycle.amplification_factor,
                    "cycle_type": cycle.cycle_type,
                    "is_self_sustaining": cycle.is_self_sustaining
                })
            
            logger.info(f"Detected {len(cycles)} cycles in {run_dir.name}")
            return cycles
        
        except ImportError:
            logger.warning("autocatalytic_detector not available, trying alternative")
            # Try alternative detector
            try:
                from backend.sim.analysis.autocatalysis_detector import AutocatalysisDetector
                # This would require different interface - skip for now
                logger.warning("Alternative detector interface not implemented")
                return []
            except ImportError:
                logger.error("No autocatalysis detector available")
                return []
        except Exception as e:
            logger.error(f"Failed to detect cycles in {run_dir.name}: {e}")
            return []
    
    def _build_network_graph(
        self, 
        cycles: List[Dict], 
        include_metrics: bool
    ) -> 'nx.DiGraph':
        """
        Build NetworkX graph from cycles.
        
        Nodes: molecules in cycles
        Edges: reactions between molecules
        Attributes: amplification_factor, cycle_type, etc.
        
        Args:
            cycles: List of cycle dicts
            include_metrics: Whether to include metrics
        
        Returns:
            NetworkX DiGraph
        """
        if not NETWORKX_AVAILABLE:
            raise ImportError("networkx required")
        
        graph = nx.DiGraph()
        
        for cycle in cycles:
            molecules = cycle.get("molecules", [])
            reactions = cycle.get("reactions", [])
            
            if not molecules:
                continue
            
            # Add nodes
            for mol in molecules:
                if not graph.has_node(mol):
                    node_data = {"formula": mol}
                    if include_metrics:
                        # Could add more node metrics here
                        pass
                    graph.add_node(mol, **node_data)
            
            # Add edges (reactions)
            for i in range(len(molecules)):
                source = molecules[i]
                target = molecules[(i + 1) % len(molecules)]
                
                edge_data = {
                    "reaction_string": reactions[i] if i < len(reactions) else "",
                    "amplification_factor": cycle.get("amplification_factor", 1.0),
                    "cycle_type": cycle.get("cycle_type", "unknown")
                }
                
                if graph.has_edge(source, target):
                    # Merge edge data (multiple cycles can share edges)
                    existing = graph[source][target]
                    if "cycle_types" not in existing:
                        existing["cycle_types"] = []
                    existing["cycle_types"].append(cycle.get("cycle_type", "unknown"))
                else:
                    edge_data["cycle_types"] = [cycle.get("cycle_type", "unknown")]
                    graph.add_edge(source, target, **edge_data)
        
        return graph

