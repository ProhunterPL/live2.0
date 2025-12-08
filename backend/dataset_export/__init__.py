"""
Dataset Export Module for Live 2.0.

Provides functionality to export simulation data as datasets:
- Reaction trajectories
- Autocatalysis networks
- Novel molecules

Example usage:
    from backend.dataset_export import DatasetExporter
    
    exporter = DatasetExporter("results/phase2b_additional")
    
    # Export reaction trajectories
    trajectories_file = exporter.export_reaction_trajectories(
        runs=["miller_urey_extended/run_*"],
        output_format="parquet"
    )
    
    # Export autocatalysis network
    network_file = exporter.export_autocatalysis_network(
        runs=["miller_urey_extended/run_*"],
        output_format="graphml"
    )
    
    # Export novel molecules
    molecules_file = exporter.export_novel_molecules(
        runs=["miller_urey_extended/run_*"],
        novelty_threshold=0.7,
        limit=776
    )
"""

from backend.dataset_export.exporter import DatasetExporter

__all__ = ["DatasetExporter"]

