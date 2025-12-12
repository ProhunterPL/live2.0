"""
Integration tests for API v1 routes.
"""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import Mock, AsyncMock

from backend.api.v1.main import app
from backend.api.v1.auth import User, verify_api_key
from backend.api.v1.dependencies import get_job_processor, get_rate_limiter


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
    # Mock dependencies using FastAPI's dependency_overrides
    mock_job_processor = Mock()
    mock_job_processor.create_job = AsyncMock(return_value="job_abc123")
    
    mock_rate_limiter = Mock()
    mock_rate_limiter.check_quota = Mock(return_value={"allowed": True})
    mock_rate_limiter.track_usage = Mock()
    
    # Override dependencies
    app.dependency_overrides[verify_api_key] = lambda: mock_user
    app.dependency_overrides[get_job_processor] = lambda: mock_job_processor
    app.dependency_overrides[get_rate_limiter] = lambda: mock_rate_limiter
    
    try:
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
    finally:
        # Clean up dependency overrides
        app.dependency_overrides.clear()


def test_get_job_status(client, mock_user):
    """Test get job status endpoint."""
    # Mock dependencies using FastAPI's dependency_overrides
    mock_job_processor = Mock()
    mock_job_processor.get_job_status = AsyncMock(return_value={
        "job_id": "job_abc123",
        "status": "running",
        "progress": 50,
        "user_id": "user_123"
    })
    
    # Override dependencies
    app.dependency_overrides[verify_api_key] = lambda: mock_user
    app.dependency_overrides[get_job_processor] = lambda: mock_job_processor
    
    try:
        response = client.get(
            "/jobs/job_abc123",
            headers={"X-API-Key": "test_key"}
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["job_id"] == "job_abc123"
        assert data["status"] == "running"
    finally:
        # Clean up dependency overrides
        app.dependency_overrides.clear()

