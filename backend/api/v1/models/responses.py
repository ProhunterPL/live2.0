"""
Response schemas for API v1.
"""

from pydantic import BaseModel
from typing import List, Optional, Dict
from datetime import datetime


class JobResponse(BaseModel):
    """Response schema for job creation."""
    job_id: str
    status: str
    estimated_time: Optional[int] = None


class JobStatusResponse(BaseModel):
    """Response schema for job status."""
    job_id: str
    status: str
    progress: int  # 0-100
    result_url: Optional[str] = None
    error: Optional[str] = None
    created_at: str  # ISO format datetime
    completed_at: Optional[str] = None  # ISO format datetime


class MoleculeResponse(BaseModel):
    """Response schema for molecule."""
    id: str
    formula: str
    smiles: Optional[str] = None
    novelty_score: float
    complexity: int
    discovery_conditions: Dict


class MoleculesListResponse(BaseModel):
    """Response schema for molecules list."""
    molecules: List[MoleculeResponse]
    total: int
    limit: int
    offset: int


class ReactionResponse(BaseModel):
    """Response schema for reaction."""
    reaction_id: str
    reactants: List[str]
    products: List[str]
    reaction_type: str
    step: int
    temperature: float
    scenario: str


class ReactionsListResponse(BaseModel):
    """Response schema for reactions list."""
    reactions: List[ReactionResponse]
    total: int
    limit: int
    offset: int


class PredictReactionResponse(BaseModel):
    """Response schema for reaction prediction."""
    products: List[str]
    probability: float
    pathway: Optional[List[Dict]] = None
    confidence: float

