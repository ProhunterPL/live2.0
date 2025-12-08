"""
Request schemas for API v1.
"""

from pydantic import BaseModel, Field, HttpUrl
from typing import List, Optional, Dict, Literal


class GenerateDatasetRequest(BaseModel):
    """Request schema for dataset generation."""
    dataset_type: Literal["reaction_trajectories", "autocatalysis_network", "novel_molecules"]
    params: Dict
    webhook_url: Optional[HttpUrl] = None


class RunSimulationRequest(BaseModel):
    """Request schema for simulation run."""
    initial_conditions: Dict
    steps: int = Field(..., gt=0, le=10000000)
    params: Optional[Dict] = None


class PredictReactionRequest(BaseModel):
    """Request schema for reaction prediction."""
    reactants: List[str] = Field(..., min_items=1)
    conditions: Dict

