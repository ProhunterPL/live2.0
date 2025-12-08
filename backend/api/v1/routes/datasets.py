"""
Dataset generation endpoints.
"""

from fastapi import APIRouter, Depends, HTTPException
from typing import Dict

from backend.api.v1.models.requests import GenerateDatasetRequest
from backend.api.v1.models.responses import JobResponse
from backend.api.v1.auth import verify_api_key, User
from backend.api.v1.dependencies import get_job_processor, get_rate_limiter
from backend.api.v1.jobs import JobProcessor
from backend.api.v1.rate_limiter import RateLimiter

router = APIRouter(prefix="/datasets", tags=["datasets"])


@router.post("/generate", response_model=JobResponse)
async def generate_dataset(
    request: GenerateDatasetRequest,
    user: User = Depends(verify_api_key),
    job_processor: JobProcessor = Depends(get_job_processor),
    rate_limiter: RateLimiter = Depends(get_rate_limiter)
):
    """
    Generate dataset asynchronously.
    
    Creates a job that will process dataset generation in background.
    """
    # Check quota (dataset generation counts as reactions)
    rate_limiter.check_quota(user, resource_type="reactions")
    
    # Validate params
    params = request.params
    if "runs" not in params:
        raise HTTPException(status_code=400, detail="Missing 'runs' in params")
    
    # Create job
    job_id = await job_processor.create_job(
        job_type="generate_dataset",
        params={
            "dataset_type": request.dataset_type,
            **params
        },
        user_id=user.id,
        webhook_url=str(request.webhook_url) if request.webhook_url else None
    )
    
    # Track usage
    rate_limiter.track_usage(user, resource_type="reactions", amount=1)
    
    # Estimate time (rough estimate: 1 hour for large datasets)
    estimated_time = 3600
    
    return JobResponse(
        job_id=job_id,
        status="queued",
        estimated_time=estimated_time
    )

