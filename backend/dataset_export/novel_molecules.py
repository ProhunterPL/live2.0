"""
Novel Molecule Exporter for dataset export.

Exports novel molecules with novelty scores.
Uses MoleculeExtractor to get molecules, then calculates novelty scores
and filters by threshold.
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional

from backend.dataset_export.formatters import FormatterFactory
from backend.dataset_export.utils import resolve_run_patterns, load_run_config

logger = logging.getLogger(__name__)


class NovelMoleculeExporter:
    """
    Exports novel molecules with novelty scores.
    
    Uses MoleculeExtractor to get molecules, then calculates novelty scores
    and filters by threshold.
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
        novelty_threshold: float = 0.7,
        limit: int = None,
        output_format: str = "json",
        output_path: str = None,
        include_graphs: bool = True
    ) -> str:
        """
        Export novel molecules.
        
        Process:
        1. Extract molecules from all runs using MoleculeExtractor
        2. Deduplicate molecules by formula
        3. Calculate novelty scores for each molecule
        4. Filter by threshold and sort by score
        5. Apply limit (e.g., top 776)
        6. Export as JSON or CSV
        
        Args:
            runs: List of run patterns
            novelty_threshold: Minimum novelty score (0.0-1.0)
            limit: Maximum number of molecules (None = all, must be positive integer if provided)
            output_format: "json" or "csv"
            output_path: Path to output file
            include_graphs: Whether to include graph representations
        
        Returns:
            Path to exported file
        """
        if not output_path:
            raise ValueError("output_path is required")
        
        # Validate novelty_threshold
        if not (0.0 <= novelty_threshold <= 1.0):
            raise ValueError(f"novelty_threshold must be between 0.0 and 1.0, got {novelty_threshold}")
        
        # Validate limit
        if limit is not None:
            if not isinstance(limit, int) or limit <= 0:
                raise ValueError(f"limit must be a positive integer, got {limit}")
        
        # Resolve runs
        run_dirs = resolve_run_patterns(runs, self.base_results_dir)
        
        if not run_dirs:
            logger.warning("No run directories found")
            return output_path
        
        # Extract all molecules
        all_molecules = []
        
        for run_dir in run_dirs:
            try:
                molecules = self._extract_molecules_from_run(run_dir)
                # Add discovery conditions to each molecule
                config = load_run_config(run_dir)
                for mol in molecules:
                    mol["discovery_conditions"] = {
                        "scenario": run_dir.parent.name,
                        "run_id": run_dir.name,
                        "step": mol.get("first_seen", 0),
                        "temperature": config.get("temperature", 298.0)
                    }
                all_molecules.extend(molecules)
            except Exception as e:
                logger.error(f"Failed to extract molecules from {run_dir.name}: {e}")
                continue
        
        if not all_molecules:
            logger.warning("No molecules found")
            # Create empty export
            formatter = FormatterFactory.create(output_format)
            formatter.export_molecules([], output_path, include_graphs)
            return output_path
        
        # Deduplicate molecules by formula (merge from different runs)
        deduplicated_molecules = self._deduplicate_molecules(all_molecules)
        logger.info(f"Deduplicated {len(all_molecules)} molecules to {len(deduplicated_molecules)} unique formulas")
        
        # Calculate novelty scores
        scored_molecules = self._calculate_novelty_scores(deduplicated_molecules)
        
        # Filter and sort
        filtered = [
            m for m in scored_molecules 
            if m.get("novelty_score", 0.0) >= novelty_threshold
        ]
        filtered.sort(key=lambda x: x.get("novelty_score", 0.0), reverse=True)
        
        # Apply limit
        if limit is not None:
            filtered = filtered[:limit]
        
        logger.info(f"Exporting {len(filtered)} molecules (threshold={novelty_threshold}, limit={limit})")
        
        # Export
        formatter = FormatterFactory.create(output_format)
        formatter.export_molecules(filtered, output_path, include_graphs)
        
        return output_path
    
    def _extract_molecules_from_run(self, run_dir: Path) -> List[Dict]:
        """
        Extract molecules from a run using MoleculeExtractor.
        
        Tries:
        1. molecules.json (if exists and not empty)
        2. MoleculeExtractor.extract_from_snapshots()
        3. MoleculeExtractor.extract_from_final_state()
        
        Args:
            run_dir: Path to run directory
        
        Returns:
            List of molecule dicts
        """
        # Try molecules.json first
        molecules_file = run_dir / "molecules.json"
        if molecules_file.exists():
            try:
                with open(molecules_file) as f:
                    molecules = json.load(f)
                    if molecules and isinstance(molecules, list):
                        logger.debug(f"Loaded {len(molecules)} molecules from {molecules_file.name}")
                        return molecules
            except Exception as e:
                logger.warning(f"Failed to load {molecules_file.name}: {e}")
        
        # Try MoleculeExtractor
        try:
            from backend.sim.molecule_extractor import MoleculeExtractor
            
            extractor = MoleculeExtractor(str(run_dir))
            
            # Try snapshots
            molecules = extractor.extract_from_snapshots(snapshot_interval=10)
            if molecules:
                logger.debug(f"Extracted {len(molecules)} molecules from snapshots in {run_dir.name}")
                return molecules
            
            # Try final state
            molecules = extractor.extract_from_final_state()
            if molecules:
                logger.debug(f"Extracted {len(molecules)} molecules from final state in {run_dir.name}")
                return molecules
        
        except Exception as e:
            logger.warning(f"MoleculeExtractor failed for {run_dir.name}: {e}")
        
        # Try results.json
        results_file = run_dir / "results.json"
        if results_file.exists():
            try:
                with open(results_file) as f:
                    data = json.load(f)
                molecules = data.get("molecules_detected", [])
                if molecules:
                    logger.debug(f"Loaded {len(molecules)} molecules from results.json in {run_dir.name}")
                    return molecules
            except Exception as e:
                logger.warning(f"Failed to load molecules from results.json: {e}")
        
        logger.debug(f"No molecules found in {run_dir.name}")
        return []
    
    def _deduplicate_molecules(self, molecules: List[Dict]) -> List[Dict]:
        """
        Deduplicate molecules by formula, merging data from different runs.
        
        For molecules with the same formula:
        - Sum counts
        - Take maximum complexity
        - Take maximum mass
        - Merge discovery_conditions (keep earliest first_seen)
        
        Args:
            molecules: List of molecule dicts (may contain duplicates)
        
        Returns:
            List of deduplicated molecules
        """
        from collections import defaultdict
        
        # Group by formula
        formula_groups = defaultdict(list)
        for mol in molecules:
            formula = mol.get("formula", "UNKNOWN")
            formula_groups[formula].append(mol)
        
        # Merge molecules with same formula
        deduplicated = []
        for formula, mol_list in formula_groups.items():
            if len(mol_list) == 1:
                deduplicated.append(mol_list[0])
            else:
                # Merge multiple instances
                merged = mol_list[0].copy()
                
                # Sum counts
                merged["count"] = sum(m.get("count", 0) for m in mol_list)
                
                # Take maximum complexity and mass
                merged["complexity"] = max(m.get("complexity", 0) for m in mol_list)
                merged["mass"] = max(m.get("mass", 0.0) for m in mol_list)
                
                # Merge discovery_conditions - keep earliest first_seen
                all_conditions = [m.get("discovery_conditions", {}) for m in mol_list]
                earliest_condition = min(
                    all_conditions,
                    key=lambda c: c.get("step", float('inf'))
                )
                merged["discovery_conditions"] = earliest_condition
                
                # Merge run_ids if present
                run_ids = [c.get("run_id", "") for c in all_conditions if c.get("run_id")]
                if run_ids:
                    merged["discovery_conditions"]["all_run_ids"] = list(set(run_ids))
                
                deduplicated.append(merged)
        
        return deduplicated
    
    def _calculate_novelty_scores(self, molecules: List[Dict]) -> List[Dict]:
        """
        Calculate novelty scores for molecules.
        
        Formula (from scripts/plot_top_molecules.py):
        novelty_score = count * (1 + complexity) * (1 + mass/100)
        
        Normalized to 0.0-1.0 range.
        
        Args:
            molecules: List of molecule dicts
        
        Returns:
            List of molecules with novelty_score added
        """
        if not molecules:
            return []
        
        # Calculate raw scores
        for mol in molecules:
            count = mol.get("count", 1)
            complexity = mol.get("complexity", 0)
            mass = mol.get("mass", 0.0)
            
            raw_score = count * (1 + complexity) * (1 + mass / 100)
            mol["novelty_score_raw"] = raw_score
        
        # Normalize to 0.0-1.0
        max_score = max(m.get("novelty_score_raw", 0.0) for m in molecules)
        if max_score > 0:
            for mol in molecules:
                mol["novelty_score"] = mol.get("novelty_score_raw", 0.0) / max_score
        else:
            for mol in molecules:
                mol["novelty_score"] = 0.0
        
        return molecules

