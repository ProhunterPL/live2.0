"""
Reaction prediction endpoints (placeholder).
"""

from fastapi import APIRouter, Depends, HTTPException

from backend.api.v1.models.requests import PredictReactionRequest
from backend.api.v1.models.responses import PredictReactionResponse
from backend.api.v1.auth import verify_api_key, User

router = APIRouter(prefix="/predictions", tags=["predictions"])


@router.post("/reaction", response_model=PredictReactionResponse)
async def predict_reaction(
    request: PredictReactionRequest,
    user: User = Depends(verify_api_key)
):
    """
    Predict reaction (placeholder - not yet implemented).
    
    This endpoint will be implemented in future version.
    """
    raise HTTPException(
        status_code=501,
        detail="Reaction prediction not yet implemented"
    )

