"""
Main Dataset Exporter for Live 2.0.

Coordinates export of reaction trajectories, autocatalysis networks,
and novel molecules from simulation results.
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional, Callable

from backend.dataset_export.reaction_trajectories import ReactionTrajectoryExporter
from backend.dataset_export.autocatalysis_network import AutocatalysisNetworkExporter
from backend.dataset_export.novel_molecules import NovelMoleculeExporter

logger = logging.getLogger(__name__)


class DatasetExporter:
    """
    Main exporter class for Live 2.0 datasets.
    
    Coordinates export of reaction trajectories, autocatalysis networks,
    and novel molecules from simulation results.
    """
    
    def __init__(self, base_results_dir: str, output_dir: str = "datasets"):
        """
        Initialize exporter.
        
        Args:
            base_results_dir: Base directory with results 
                            (e.g., "results/phase2b_additional")
            output_dir: Directory for exported datasets (default: "datasets")
        """
        self.base_results_dir = Path(base_results_dir)
        if not self.base_results_dir.exists():
            raise ValueError(f"Base results directory does not exist: {base_results_dir}")
        
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize sub-exporters
        self.reaction_exporter = ReactionTrajectoryExporter(str(self.base_results_dir))
        self.autocatalysis_exporter = AutocatalysisNetworkExporter(str(self.base_results_dir))
        self.molecule_exporter = NovelMoleculeExporter(str(self.base_results_dir))
        
        logger.info(f"DatasetExporter initialized: base_dir={self.base_results_dir}, output_dir={self.output_dir}")
    
    def export_reaction_trajectories(
        self,
        runs: List[str],  # Patterns like "miller_urey_extended/run_*"
        output_format: str = "parquet",  # "parquet" | "json"
        filters: Optional[Dict] = None,  # {"scenario": "...", "min_steps": 100000}
        output_path: Optional[str] = None,
        progress_callback: Optional[Callable] = None
    ) -> str:
        """
        Export reaction trajectories to dataset.
        
        Args:
            runs: List of run patterns (e.g., ["miller_urey_extended/run_*"])
            output_format: "parquet" or "json"
            filters: Optional filters (scenario, min_steps, temperature_range)
            output_path: Path to output file (default: output_dir/LIVE-Reaction-Trajectories-1M.{format})
            progress_callback: Optional callback(current, total, message)
        
        Returns:
            Path to exported dataset file
        """
        if not output_path:
            ext = "parquet" if output_format == "parquet" else "json"
            output_path = str(self.output_dir / f"LIVE-Reaction-Trajectories-1M.{ext}")
        
        logger.info(f"Exporting reaction trajectories: {len(runs)} run pattern(s), format={output_format}")
        
        return self.reaction_exporter.export(
            runs=runs,
            output_format=output_format,
            filters=filters,
            output_path=output_path,
            progress_callback=progress_callback
        )
    
    def export_autocatalysis_network(
        self,
        runs: List[str],
        output_format: str = "json",  # "json" | "graphml"
        output_path: Optional[str] = None,
        include_metrics: bool = True
    ) -> str:
        """
        Export autocatalysis cycles as network graph.
        
        Args:
            runs: List of run patterns
            output_format: "json" or "graphml"
            output_path: Path to output file (default: output_dir/LIVE-Autocatalysis-Network-Set.{format})
            include_metrics: Whether to include metrics in graph
        
        Returns:
            Path to exported network file
        """
        if not output_path:
            ext = "graphml" if output_format == "graphml" else "json"
            output_path = str(self.output_dir / f"LIVE-Autocatalysis-Network-Set.{ext}")
        
        logger.info(f"Exporting autocatalysis network: {len(runs)} run pattern(s), format={output_format}")
        
        return self.autocatalysis_exporter.export(
            runs=runs,
            output_format=output_format,
            output_path=output_path,
            include_metrics=include_metrics
        )
    
    def export_novel_molecules(
        self,
        runs: List[str],
        novelty_threshold: float = 0.7,
        limit: int = None,
        output_format: str = "json",  # "json" | "csv"
        output_path: Optional[str] = None,
        include_graphs: bool = True
    ) -> str:
        """
        Export novel molecules with scores.
        
        Args:
            runs: List of run patterns
            novelty_threshold: Minimum novelty score (0.0-1.0)
            limit: Maximum number of molecules (None = all, default: 776 for "LIVE-Novel-Molecules-776")
            output_format: "json" or "csv"
            output_path: Path to output file (default: output_dir/LIVE-Novel-Molecules-{limit or 'all'}.{format})
            include_graphs: Whether to include graph representations
        
        Returns:
            Path to exported molecules file
        """
        if not output_path:
            ext = "csv" if output_format == "csv" else "json"
            limit_str = str(limit) if limit else "all"
            output_path = str(self.output_dir / f"LIVE-Novel-Molecules-{limit_str}.{ext}")
        
        logger.info(f"Exporting novel molecules: {len(runs)} run pattern(s), threshold={novelty_threshold}, limit={limit}")
        
        return self.molecule_exporter.export(
            runs=runs,
            novelty_threshold=novelty_threshold,
            limit=limit,
            output_format=output_format,
            output_path=output_path,
            include_graphs=include_graphs
        )

