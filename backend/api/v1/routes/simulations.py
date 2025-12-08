"""
Simulation endpoints (placeholder).
"""

from fastapi import APIRouter, Depends, HTTPException
from typing import Dict

from backend.api.v1.models.requests import RunSimulationRequest
from backend.api.v1.auth import verify_api_key, User

router = APIRouter(prefix="/simulations", tags=["simulations"])


@router.post("/run")
async def run_simulation(
    request: RunSimulationRequest,
    user: User = Depends(verify_api_key)
):
    """
    Run simulation (placeholder - not yet implemented).
    
    This endpoint will be implemented in future version.
    """
    raise HTTPException(
        status_code=501,
        detail="Simulation jobs not yet implemented"
    )

