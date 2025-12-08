"""
Query logic for molecules and reactions.

Provides functions to query data from simulation results.
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional
from functools import lru_cache

from backend.dataset_export.utils import resolve_run_patterns, load_run_config
from backend.dataset_export.novel_molecules import NovelMoleculeExporter
from backend.dataset_export.reaction_trajectories import ReactionTrajectoryExporter

logger = logging.getLogger(__name__)


class MoleculeQuery:
    """Query molecules from simulation results."""
    
    def __init__(self, base_results_dir: str):
        """
        Initialize query.
        
        Args:
            base_results_dir: Base directory with results
        """
        self.base_results_dir = Path(base_results_dir)
        self.molecule_exporter = NovelMoleculeExporter(str(self.base_results_dir))
    
    def query(
        self,
        novelty_min: Optional[float] = None,
        complexity_max: Optional[int] = None,
        scenario: Optional[str] = None,
        limit: int = 100,
        offset: int = 0
    ) -> Dict:
        """
        Query molecules with filtering and pagination.
        
        Args:
            novelty_min: Minimum novelty score (0.0-1.0)
            complexity_max: Maximum complexity
            scenario: Filter by scenario name
            limit: Maximum number of results
            offset: Offset for pagination
        
        Returns:
            Dict with:
            - molecules: List[Dict]
            - total: int
            - limit: int
            - offset: int
        """
        # Get all runs (or filter by scenario)
        if scenario:
            runs = [f"{scenario}/run_*"]
        else:
            # Get all scenarios
            runs = []
            for scenario_dir in self.base_results_dir.iterdir():
                if scenario_dir.is_dir():
                    runs.append(f"{scenario_dir.name}/run_*")
        
        if not runs:
            return {
                "molecules": [],
                "total": 0,
                "limit": limit,
                "offset": offset
            }
        
        # Extract all molecules
        all_molecules = []
        for run_pattern in runs:
            try:
                run_dirs = resolve_run_patterns([run_pattern], self.base_results_dir)
                for run_dir in run_dirs:
                    molecules = self.molecule_exporter._extract_molecules_from_run(run_dir)
                    
                    # Add discovery conditions
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
                logger.warning(f"Failed to extract molecules from {run_pattern}: {e}")
                continue
        
        # Deduplicate
        all_molecules = self.molecule_exporter._deduplicate_molecules(all_molecules)
        
        # Calculate novelty scores
        all_molecules = self.molecule_exporter._calculate_novelty_scores(all_molecules)
        
        # Apply filters
        filtered = all_molecules
        if novelty_min is not None:
            filtered = [m for m in filtered if m.get("novelty_score", 0.0) >= novelty_min]
        if complexity_max is not None:
            filtered = [m for m in filtered if m.get("complexity", 0) <= complexity_max]
        
        # Sort by novelty score (descending)
        filtered.sort(key=lambda x: x.get("novelty_score", 0.0), reverse=True)
        
        # Pagination
        total = len(filtered)
        paginated = filtered[offset:offset + limit]
        
        # Format for API response
        molecules = []
        for mol in paginated:
            molecules.append({
                "id": f"mol_{mol.get('formula', 'UNKNOWN').replace(' ', '_')}",
                "formula": mol.get("formula", "UNKNOWN"),
                "smiles": mol.get("smiles"),
                "novelty_score": mol.get("novelty_score", 0.0),
                "complexity": mol.get("complexity", 0),
                "discovery_conditions": mol.get("discovery_conditions", {})
            })
        
        return {
            "molecules": molecules,
            "total": total,
            "limit": limit,
            "offset": offset
        }


class ReactionQuery:
    """Query reactions from simulation results."""
    
    def __init__(self, base_results_dir: str):
        """
        Initialize query.
        
        Args:
            base_results_dir: Base directory with results
        """
        self.base_results_dir = Path(base_results_dir)
        self.reaction_exporter = ReactionTrajectoryExporter(str(self.base_results_dir))
    
    def query(
        self,
        type: Optional[str] = None,
        scenario: Optional[str] = None,
        temperature_min: Optional[float] = None,
        temperature_max: Optional[float] = None,
        limit: int = 100,
        offset: int = 0
    ) -> Dict:
        """
        Query reactions with filtering and pagination.
        
        Args:
            type: Filter by reaction type (formation, breakage, condensation, etc.)
            scenario: Filter by scenario name
            temperature_min: Minimum temperature
            temperature_max: Maximum temperature
            limit: Maximum number of results
            offset: Offset for pagination
        
        Returns:
            Dict with:
            - reactions: List[Dict]
            - total: int
            - limit: int
            - offset: int
        """
        # Get all runs (or filter by scenario)
        if scenario:
            runs = [f"{scenario}/run_*"]
        else:
            # Get all scenarios
            runs = []
            for scenario_dir in self.base_results_dir.iterdir():
                if scenario_dir.is_dir():
                    runs.append(f"{scenario_dir.name}/run_*")
        
        if not runs:
            return {
                "reactions": [],
                "total": 0,
                "limit": limit,
                "offset": offset
            }
        
        # Extract all reactions
        all_reactions = []
        for run_pattern in runs:
            try:
                run_dirs = resolve_run_patterns([run_pattern], self.base_results_dir)
                for run_dir in run_dirs:
                    reactions = self.reaction_exporter._extract_reactions_from_run(run_dir)
                    all_reactions.extend(reactions)
            except Exception as e:
                logger.warning(f"Failed to extract reactions from {run_pattern}: {e}")
                continue
        
        # Apply filters
        filtered = all_reactions
        if type:
            filtered = [r for r in filtered if r.get("reaction_type") == type]
        if scenario:
            filtered = [r for r in filtered if r.get("scenario") == scenario]
        if temperature_min is not None:
            filtered = [r for r in filtered if r.get("temperature", 0.0) >= temperature_min]
        if temperature_max is not None:
            filtered = [r for r in filtered if r.get("temperature", 0.0) <= temperature_max]
        
        # Sort by step (descending - most recent first)
        filtered.sort(key=lambda x: x.get("step", 0), reverse=True)
        
        # Pagination
        total = len(filtered)
        paginated = filtered[offset:offset + limit]
        
        # Format for API response
        reactions = []
        for rxn in paginated:
            reactions.append({
                "reaction_id": rxn.get("reaction_id", "UNKNOWN"),
                "reactants": rxn.get("reactants", []),
                "products": rxn.get("products", []),
                "reaction_type": rxn.get("reaction_type", "unknown"),
                "step": rxn.get("step", 0),
                "temperature": rxn.get("temperature", 298.0),
                "scenario": rxn.get("scenario", "unknown")
            })
        
        return {
            "reactions": reactions,
            "total": total,
            "limit": limit,
            "offset": offset
        }

