"""
Tests for API v1 job processor.
"""

import pytest
import json
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, AsyncMock
from datetime import datetime

from backend.api.v1.jobs import JobProcessor, JobStatus


def test_create_job():
    """Test job creation."""
    # Create temporary directory for test results
    with tempfile.TemporaryDirectory() as tmpdir:
        mock_redis = Mock()
        mock_redis.setex = Mock()
        mock_redis.lpush = Mock()
        
        # Mock _start_processor to avoid event loop issues
        with patch.object(JobProcessor, '_start_processor', return_value=None):
            processor = JobProcessor(
                redis_client=mock_redis,
                base_results_dir=tmpdir,
                storage_manager=None
            )
            
            # Note: create_job is async, but we're testing the sync initialization
            # The actual job creation would need to be tested in an async test
            assert processor.redis == mock_redis
            assert processor.base_results_dir == tmpdir


@pytest.mark.asyncio
async def test_get_job_status_found():
    """Test get_job_status returns job data when found."""
    # Create temporary directory for test results
    with tempfile.TemporaryDirectory() as tmpdir:
        mock_redis = Mock()
        job_data = {
            "job_id": "job_abc123",
            "status": "running",
            "progress": 50
        }
        mock_redis.get.return_value = json.dumps(job_data)
        
        processor = JobProcessor(
            redis_client=mock_redis,
            base_results_dir=tmpdir,
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
    
    # Create temporary directory for test results
    with tempfile.TemporaryDirectory() as tmpdir:
        mock_redis = Mock()
        mock_redis.get.return_value = None
        
        processor = JobProcessor(
            redis_client=mock_redis,
            base_results_dir=tmpdir,
            storage_manager=None
        )
        
        with pytest.raises(HTTPException) as exc_info:
            await processor.get_job_status("job_nonexistent")
        
        assert exc_info.value.status_code == 404
        assert "not found" in exc_info.value.detail.lower()

