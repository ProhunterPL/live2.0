"""
Job status endpoints.
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from typing import Dict, List
import json
import os

from backend.api.v1.models.responses import JobStatusResponse, JobResponse
from backend.api.v1.models.requests import StartJobRequest
from backend.api.v1.auth import verify_api_key, User
from backend.api.v1.dependencies import get_job_processor, get_aws_batch_client, get_rate_limiter
from backend.api.v1.jobs import JobProcessor
from backend.api.v1.aws_batch import AWSBatchClient
from backend.api.v1.rate_limiter import RateLimiter

router = APIRouter(prefix="/jobs", tags=["jobs"])


@router.get("", response_model=List[JobStatusResponse])
async def list_jobs(
    limit: int = Query(10, ge=1, le=100),
    user: User = Depends(verify_api_key),
    job_processor: JobProcessor = Depends(get_job_processor)
):
    """
    List all jobs for the authenticated user.
    
    Returns list of jobs sorted by creation date (newest first).
    """
    try:
        jobs_data = await job_processor.list_user_jobs(user.id, limit=limit)
        
        return [
            JobStatusResponse(
                job_id=job["job_id"],
                status=job["status"],
                progress=job.get("progress", 0),
                result_url=job.get("result_url"),
                error=job.get("error"),
                created_at=job.get("created_at", ""),
                completed_at=job.get("completed_at")
            )
            for job in jobs_data
        ]
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to list jobs: {str(e)}")


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


@router.post("", response_model=JobResponse)
async def start_job(
    request: StartJobRequest,
    user: User = Depends(verify_api_key),
    job_processor: JobProcessor = Depends(get_job_processor),
    rate_limiter: RateLimiter = Depends(get_rate_limiter)
):
    """
    Start a new job (simulation or dataset generation).
    
    For simulation jobs, submits to AWS Batch.
    For dataset jobs, queues locally.
    """
    # Check quota
    resource_type = "simulations" if request.job_type == "run_simulation" else "reactions"
    rate_limiter.check_quota(user, resource_type=resource_type)
    
    # Check idempotency if provided
    if request.idempotency_key:
        idempotency_key = f"idempotency:{user.id}:{request.idempotency_key}"
        existing_job = job_processor.redis.get(idempotency_key)
        if existing_job:
            job_data = json.loads(existing_job)
            return JobResponse(
                job_id=job_data["job_id"],
                status=job_data["status"],
                estimated_time=3600
            )
    
    # Create job
    job_id = await job_processor.create_job(
        job_type=request.job_type,
        params=request.params,
        user_id=user.id,
        webhook_url=str(request.webhook_url) if request.webhook_url else None
    )
    
    # Store idempotency key
    if request.idempotency_key:
        idempotency_key = f"idempotency:{user.id}:{request.idempotency_key}"
        job_processor.redis.setex(idempotency_key, 86400, json.dumps({"job_id": job_id}))
    
    # Track usage
    rate_limiter.track_usage(user, resource_type=resource_type, amount=1)
    
    return JobResponse(
        job_id=job_id,
        status="queued",
        estimated_time=3600  # 1 hour estimate
    )


@router.post("/{job_id}/cancel", response_model=JobStatusResponse)
async def cancel_job(
    job_id: str,
    user: User = Depends(verify_api_key),
    job_processor: JobProcessor = Depends(get_job_processor)
):
    """
    Cancel a running job.
    
    Cancels AWS Batch job if running, updates status to cancelled.
    """
    await job_processor.cancel_job(job_id, user.id)
    
    # Get updated status
    job_data = await job_processor.get_job_status(job_id)
    
    return JobStatusResponse(
        job_id=job_data["job_id"],
        status=job_data["status"],
        progress=job_data.get("progress", 0),
        result_url=job_data.get("result_url"),
        error=job_data.get("error"),
        created_at=job_data.get("created_at", ""),
        completed_at=job_data.get("completed_at")
    )


@router.get("/{job_id}/artifacts")
async def list_job_artifacts(
    job_id: str,
    user: User = Depends(verify_api_key),
    job_processor: JobProcessor = Depends(get_job_processor),
    batch_client: AWSBatchClient = Depends(get_aws_batch_client)
):
    """
    List artifacts (S3 objects) for a completed job.
    
    Returns presigned URLs for downloading artifacts.
    """
    # Verify job exists and user owns it
    job_data = await job_processor.get_job_status(job_id)
    if job_data.get("user_id") != user.id:
        raise HTTPException(status_code=403, detail="Access denied to this job")
    
    # List artifacts from S3
    artifacts = batch_client.list_job_artifacts(
        job_id=job_id,
        user_id=user.id,
        bucket=os.getenv("AWS_S3_BUCKET", "live2-artifacts")
    )
    
    # Generate presigned URLs
    artifact_urls = []
    for artifact in artifacts:
        url = batch_client.generate_presigned_url(
            bucket=os.getenv("AWS_S3_BUCKET", "live2-artifacts"),
            key=artifact["key"],
            expiration=3600  # 1 hour
        )
        artifact_urls.append({
            "key": artifact["key"],
            "size": artifact["size"],
            "url": url,
            "last_modified": artifact["last_modified"]
        })
    
    return {
        "job_id": job_id,
        "artifacts": artifact_urls
    }

