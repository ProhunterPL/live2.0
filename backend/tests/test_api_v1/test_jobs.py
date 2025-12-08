"""
Tests for API v1 job processor.
"""

import pytest
import json
from unittest.mock import Mock, patch, AsyncMock
from datetime import datetime

from backend.api.v1.jobs import JobProcessor, JobStatus


def test_create_job():
    """Test job creation."""
    mock_redis = Mock()
    mock_redis.setex = Mock()
    mock_redis.lpush = Mock()
    
    processor = JobProcessor(
        redis_client=mock_redis,
        base_results_dir="test_results",
        storage_manager=None
    )
    
    job_id = processor.create_job(
        job_type="generate_dataset",
        params={"dataset_type": "reaction_trajectories"},
        user_id="user_123"
    )
    
    assert job_id.startswith("job_")
    assert len(job_id) > 10
    
    # Verify job was stored in Redis
    mock_redis.setex.assert_called_once()
    # Verify job was queued
    mock_redis.lpush.assert_called_once_with("job_queue", job_id)


@pytest.mark.asyncio
async def test_get_job_status_found():
    """Test get_job_status returns job data when found."""
    mock_redis = Mock()
    job_data = {
        "job_id": "job_abc123",
        "status": "running",
        "progress": 50
    }
    mock_redis.get.return_value = json.dumps(job_data)
    
    processor = JobProcessor(
        redis_client=mock_redis,
        base_results_dir="test_results",
        storage_manager=None
    )
    
    result = await processor.get_job_status("job_abc123")
    
    assert result["job_id"] == "job_abc123"
    assert result["status"] == "running"
    assert result["progress"] == 50


@pytest.mark.asyncio
async def test_get_job_status_not_found():
    """Test get_job_status raises 404 when job not found."""
    from fastapi import HTTPException
    
    mock_redis = Mock()
    mock_redis.get.return_value = None
    
    processor = JobProcessor(
        redis_client=mock_redis,
        base_results_dir="test_results",
        storage_manager=None
    )
    
    with pytest.raises(HTTPException) as exc_info:
        await processor.get_job_status("job_nonexistent")
    
    assert exc_info.value.status_code == 404
    assert "not found" in exc_info.value.detail.lower()

