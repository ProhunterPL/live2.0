"""
Reaction query endpoints.
"""

from fastapi import APIRouter, Depends, Query
from typing import Optional

from backend.api.v1.models.responses import ReactionsListResponse, ReactionResponse
from backend.api.v1.auth import verify_api_key, User
from backend.api.v1.dependencies import get_rate_limiter
from backend.api.v1.rate_limiter import RateLimiter
from backend.api.v1.query import ReactionQuery
from backend.api.v1.config import BASE_RESULTS_DIR

router = APIRouter(prefix="/reactions", tags=["reactions"])

# Global query instance (singleton)
_reaction_query: Optional[ReactionQuery] = None

def get_reaction_query() -> ReactionQuery:
    """Get reaction query instance (singleton)."""
    global _reaction_query
    if _reaction_query is None:
        _reaction_query = ReactionQuery(BASE_RESULTS_DIR)
    return _reaction_query


@router.get("", response_model=ReactionsListResponse)
async def get_reactions(
    type: Optional[str] = Query(None, alias="type"),
    scenario: Optional[str] = Query(None),
    temperature_min: Optional[float] = Query(None),
    temperature_max: Optional[float] = Query(None),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    user: User = Depends(verify_api_key),
    rate_limiter: RateLimiter = Depends(get_rate_limiter),
    query: ReactionQuery = Depends(get_reaction_query)
):
    """
    Query reactions with filtering.
    
    Returns list of reactions matching criteria.
    """
    # Check quota
    rate_limiter.check_quota(user, resource_type="api_calls")
    
    # Query reactions
    result = query.query(
        type=type,
        scenario=scenario,
        temperature_min=temperature_min,
        temperature_max=temperature_max,
        limit=limit,
        offset=offset
    )
    
    # Convert to response models
    reactions = [
        ReactionResponse(**rxn) for rxn in result["reactions"]
    ]
    
    # Track usage
    rate_limiter.track_usage(user, resource_type="api_calls")
    
    return ReactionsListResponse(
        reactions=reactions,
        total=result["total"],
        limit=result["limit"],
        offset=result["offset"]
    )

