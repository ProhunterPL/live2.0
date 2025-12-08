"""
Reaction Trajectory Exporter for dataset export.

Exports reaction trajectories from simulation snapshots by comparing
consecutive snapshots and detecting bond formation/breaking events.
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional, Callable

from backend.dataset_export.snapshot_processor import SnapshotProcessor
from backend.dataset_export.formatters import FormatterFactory
from backend.dataset_export.utils import resolve_run_patterns, apply_filters, load_run_config

logger = logging.getLogger(__name__)


class ReactionTrajectoryExporter:
    """
    Exports reaction trajectories from simulation snapshots.
    
    Extracts reactions by comparing consecutive snapshots and detecting
    bond formation/breaking events.
    """
    
    def __init__(self, base_results_dir: str):
        """
        Initialize exporter.
        
        Args:
            base_results_dir: Base directory with results
        """
        self.base_results_dir = Path(base_results_dir)
        self.snapshot_processor = SnapshotProcessor()
    
    def export(
        self,
        runs: List[str],
        output_format: str = "parquet",
        filters: Optional[Dict] = None,
        output_path: str = None,
        progress_callback: Optional[Callable] = None
    ) -> str:
        """
        Export reaction trajectories.
        
        Process:
        1. Find all matching run directories
        2. For each run, load snapshots sequentially
        3. Compare consecutive snapshots to detect reactions
        4. Stream results to output file (Parquet or JSON)
        
        Args:
            runs: List of run patterns (e.g., ["miller_urey_extended/run_*"])
            output_format: "parquet" or "json"
            filters: Optional filters (scenario, min_steps, temperature_range)
            output_path: Path to output file
            progress_callback: Optional callback(current, total, message)
        
        Returns:
            Path to exported file
        """
        if not output_path:
            raise ValueError("output_path is required")
        
        # Resolve run patterns to actual directories
        run_dirs = resolve_run_patterns(runs, self.base_results_dir)
        
        if not run_dirs:
            logger.warning("No run directories found")
            return output_path
        
        # Apply filters
        run_dirs = apply_filters(run_dirs, filters)
        
        if not run_dirs:
            logger.warning("No run directories after filtering")
            return output_path
        
        # Initialize formatter
        formatter = FormatterFactory.create(output_format)
        
        # Open output file for streaming
        total_runs = len(run_dirs)
        total_reactions = 0
        
        with formatter.open_stream(output_path) as writer:
            for run_idx, run_dir in enumerate(run_dirs):
                if progress_callback:
                    progress_callback(run_idx, total_runs, f"Processing {run_dir.name}")
                
                try:
                    # Extract reactions from this run
                    reactions = self._extract_reactions_from_run(run_dir)
                    total_reactions += len(reactions)
                    
                    # Write to stream
                    for reaction in reactions:
                        writer.write(reaction)
                
                except Exception as e:
                    logger.error(f"Failed to process {run_dir.name}: {e}")
                    continue
        
        logger.info(f"Exported {total_reactions} reactions from {total_runs} runs to {output_path}")
        return output_path
    
    def _extract_reactions_from_run(self, run_dir: Path) -> List[Dict]:
        """
        Extract all reactions from a single run by comparing snapshots.
        
        Args:
            run_dir: Path to run directory
        
        Returns:
            List of reaction records (one per detected reaction)
        """
        snapshot_dir = run_dir / "snapshots"
        if not snapshot_dir.exists():
            logger.debug(f"No snapshots in {run_dir.name}")
            return []
        
        # Load results.json for metadata
        config = load_run_config(run_dir)
        
        # Get sorted snapshot files
        snapshot_files = sorted(snapshot_dir.glob("step_*.json"))
        if len(snapshot_files) < 2:
            logger.debug(f"Not enough snapshots in {run_dir.name} ({len(snapshot_files)} < 2)")
            return []
        
        reactions = []
        prev_state = None
        reaction_counter = 0  # Counter for unique reaction IDs
        
        for snapshot_file in snapshot_files:
            try:
                # Load current snapshot
                curr_state = self.snapshot_processor.load_snapshot(snapshot_file)
                
                # Validate curr_state - check if it has molecules
                if not curr_state.get("molecules"):
                    logger.debug(f"Empty snapshot {snapshot_file.name}, skipping")
                    continue
                
                if prev_state is not None:
                    # Detect reactions between states
                    detected = self._detect_reactions_between_states(
                        prev_state, curr_state, run_dir, config, reaction_counter
                    )
                    reaction_counter += len(detected)
                    reactions.extend(detected)
                
                prev_state = curr_state
            
            except Exception as e:
                logger.warning(f"Failed to process snapshot {snapshot_file.name}: {e}")
                continue
        
        return reactions
    
    def _detect_reactions_between_states(
        self,
        prev_state: Dict,
        curr_state: Dict,
        run_dir: Path,
        config: Dict,
        reaction_counter: int = 0
    ) -> List[Dict]:
        """
        Detect reactions by comparing two snapshot states.
        
        Args:
            prev_state: Previous snapshot state
            curr_state: Current snapshot state
            run_dir: Run directory (for metadata)
            config: Configuration dict
        
        Returns:
            List of reaction records
        """
        reactions = []
        
        # Get molecules from each state
        prev_molecules = {mol["formula"]: mol for mol in prev_state.get("molecules", [])}
        curr_molecules = {mol["formula"]: mol for mol in curr_state.get("molecules", [])}
        
        # Find molecules that disappeared (reactants)
        disappeared = set(prev_molecules.keys()) - set(curr_molecules.keys())
        
        # Find molecules that appeared (products)
        appeared = set(curr_molecules.keys()) - set(prev_molecules.keys())
        
        # Find molecules that changed (count or structure)
        changed = []
        for formula in set(prev_molecules.keys()) & set(curr_molecules.keys()):
            prev_mol = prev_molecules[formula]
            curr_mol = curr_molecules[formula]
            if prev_mol.get("size") != curr_mol.get("size"):
                changed.append(formula)
        
        # Create reaction records
        if disappeared or appeared or changed:
            # Simple heuristic: disappeared -> appeared
            reactants = list(disappeared) if disappeared else []
            products = list(appeared) if appeared else []
            
            # If no clear reactants/products, skip
            if not reactants and not products:
                return []
            
            # Determine reaction type
            reaction_type = self._classify_reaction_type(reactants, products)
            
            # Generate unique reaction ID using hash of reactants+products
            reaction_signature = hash(tuple(sorted(reactants + products)))
            reaction_id = f"RXN_{curr_state.get('step', 0)}_{reaction_counter}_{abs(reaction_signature) % 10000}"
            
            reaction_record = {
                "run_id": f"{run_dir.parent.name}/{run_dir.name}",
                "step": curr_state.get("step", 0),
                "time": curr_state.get("time", 0.0),
                "reaction_id": reaction_id,
                "reactants": reactants,
                "products": products,
                "reaction_type": reaction_type,
                "delta_G": None,  # Not available from snapshots
                "temperature": config.get("temperature", 298.0),
                "scenario": run_dir.parent.name
            }
            reactions.append(reaction_record)
        
        return reactions
    
    def _classify_reaction_type(self, reactants: List[str], products: List[str]) -> str:
        """
        Classify reaction type based on reactants and products.
        
        Args:
            reactants: List of reactant formulas
            products: List of product formulas
        
        Returns:
            Reaction type string
        """
        n_reactants = len(reactants)
        n_products = len(products)
        
        if n_reactants == 2 and n_products == 1:
            return "condensation"
        elif n_reactants == 1 and n_products == 2:
            return "dissociation"
        elif n_reactants == 1 and n_products == 1:
            return "rearrangement"
        elif n_reactants == 0 and n_products > 0:
            return "formation"
        elif n_reactants > 0 and n_products == 0:
            return "breakage"
        else:
            return "complex"

