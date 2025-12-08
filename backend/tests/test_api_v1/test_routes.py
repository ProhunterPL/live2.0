"""
Integration tests for API v1 routes.
"""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import Mock, patch

from backend.api.v1.main import app
from backend.api.v1.auth import User


@pytest.fixture
def client():
    """Create test client."""
    return TestClient(app)


@pytest.fixture
def mock_user():
    """Create mock user."""
    return User(
        id="user_123",
        api_key="test_key",
        tier="research",
        subscription_status="active"
    )


def test_health_check(client):
    """Test health check endpoint."""
    response = client.get("/health")
    assert response.status_code == 200
    assert response.json()["status"] == "ok"


def test_generate_dataset_missing_api_key(client):
    """Test generate_dataset requires API key."""
    response = client.post(
        "/datasets/generate",
        json={
            "dataset_type": "reaction_trajectories",
            "params": {"runs": ["test/run_1"]}
        }
    )
    assert response.status_code == 401


def test_generate_dataset_with_api_key(client, mock_user):
    """Test generate_dataset with valid API key."""
    with patch('backend.api.v1.auth.verify_api_key', return_value=mock_user):
        with patch('backend.api.v1.dependencies.get_job_processor') as mock_processor:
            mock_job_processor = Mock()
            mock_job_processor.create_job = AsyncMock(return_value="job_abc123")
            mock_processor.return_value = mock_job_processor
            
            with patch('backend.api.v1.dependencies.get_rate_limiter') as mock_limiter:
                mock_rate_limiter = Mock()
                mock_rate_limiter.check_quota = Mock(return_value={"allowed": True})
                mock_rate_limiter.track_usage = Mock()
                mock_limiter.return_value = mock_rate_limiter
                
                response = client.post(
                    "/datasets/generate",
                    json={
                        "dataset_type": "reaction_trajectories",
                        "params": {"runs": ["test/run_1"]}
                    },
                    headers={"X-API-Key": "test_key"}
                )
                
                assert response.status_code == 200
                data = response.json()
                assert "job_id" in data
                assert data["status"] == "queued"


def test_get_job_status(client, mock_user):
    """Test get job status endpoint."""
    with patch('backend.api.v1.auth.verify_api_key', return_value=mock_user):
        with patch('backend.api.v1.dependencies.get_job_processor') as mock_processor:
            mock_job_processor = Mock()
            mock_job_processor.get_job_status = AsyncMock(return_value={
                "job_id": "job_abc123",
                "status": "running",
                "progress": 50,
                "user_id": "user_123"
            })
            mock_processor.return_value = mock_job_processor
            
            response = client.get(
                "/jobs/job_abc123",
                headers={"X-API-Key": "test_key"}
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["job_id"] == "job_abc123"
            assert data["status"] == "running"

