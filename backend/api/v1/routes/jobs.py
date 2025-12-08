"""
Job status endpoints.
"""

from fastapi import APIRouter, Depends, HTTPException
from typing import Dict

from backend.api.v1.models.responses import JobStatusResponse
from backend.api.v1.auth import verify_api_key, User
from backend.api.v1.dependencies import get_job_processor
from backend.api.v1.jobs import JobProcessor

router = APIRouter(prefix="/jobs", tags=["jobs"])


@router.get("/{job_id}", response_model=JobStatusResponse)
async def get_job_status(
    job_id: str,
    user: User = Depends(verify_api_key),
    job_processor: JobProcessor = Depends(get_job_processor)
):
    """
    Get job status.
    
    Returns current status, progress, and result URL if completed.
    """
    try:
        job_data = await job_processor.get_job_status(job_id)
        
        # Verify user owns this job
        if job_data.get("user_id") != user.id:
            raise HTTPException(status_code=403, detail="Access denied to this job")
        
        return JobStatusResponse(
            job_id=job_data["job_id"],
            status=job_data["status"],
            progress=job_data.get("progress", 0),
            result_url=job_data.get("result_url"),
            error=job_data.get("error"),
            created_at=job_data.get("created_at", ""),
            completed_at=job_data.get("completed_at")
        )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get job status: {str(e)}")

