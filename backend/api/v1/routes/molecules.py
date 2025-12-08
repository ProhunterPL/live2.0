"""
Molecule query endpoints.
"""

from fastapi import APIRouter, Depends, Query
from typing import Optional

from backend.api.v1.models.responses import MoleculesListResponse, MoleculeResponse
from backend.api.v1.auth import verify_api_key, User
from backend.api.v1.dependencies import get_rate_limiter
from backend.api.v1.rate_limiter import RateLimiter
from backend.api.v1.query import MoleculeQuery
from backend.api.v1.config import BASE_RESULTS_DIR

router = APIRouter(prefix="/molecules", tags=["molecules"])

# Global query instance (singleton)
_molecule_query: Optional[MoleculeQuery] = None

def get_molecule_query() -> MoleculeQuery:
    """Get molecule query instance (singleton)."""
    global _molecule_query
    if _molecule_query is None:
        _molecule_query = MoleculeQuery(BASE_RESULTS_DIR)
    return _molecule_query


@router.get("", response_model=MoleculesListResponse)
async def get_molecules(
    novelty_min: Optional[float] = Query(None, ge=0.0, le=1.0),
    complexity_max: Optional[int] = Query(None, ge=0),
    scenario: Optional[str] = Query(None),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    user: User = Depends(verify_api_key),
    rate_limiter: RateLimiter = Depends(get_rate_limiter),
    query: MoleculeQuery = Depends(get_molecule_query)
):
    """
    Query molecules with filtering.
    
    Returns list of molecules matching criteria.
    """
    # Check quota
    rate_limiter.check_quota(user, resource_type="api_calls")
    
    # Query molecules
    result = query.query(
        novelty_min=novelty_min,
        complexity_max=complexity_max,
        scenario=scenario,
        limit=limit,
        offset=offset
    )
    
    # Convert to response models
    molecules = [
        MoleculeResponse(**mol) for mol in result["molecules"]
    ]
    
    # Track usage
    rate_limiter.track_usage(user, resource_type="api_calls")
    
    return MoleculesListResponse(
        molecules=molecules,
        total=result["total"],
        limit=result["limit"],
        offset=result["offset"]
    )

