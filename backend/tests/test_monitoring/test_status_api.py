"""
Tests for status page API.
"""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import Mock, patch, MagicMock

from backend.monitoring.status.api import router, get_response_time_tracker, get_uptime_tracker, get_error_rate_tracker


@pytest.fixture
def app():
    """Create test FastAPI app with status router."""
    from fastapi import FastAPI
    app = FastAPI()
    app.include_router(router)
    return app


def test_health_endpoint(app):
    """Test health check endpoint."""
    with TestClient(app) as client:
        response = client.get("/status/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "operational"
        assert "timestamp" in data


def test_metrics_endpoint(app):
    """Test metrics endpoint."""
    with patch("backend.monitoring.status.api.get_response_time_tracker") as mock_rt, \
         patch("backend.monitoring.status.api.get_uptime_tracker") as mock_ut, \
         patch("backend.monitoring.status.api.get_error_rate_tracker") as mock_et:
        
        # Mock trackers
        mock_rt.return_value.get_percentiles.return_value = {"p50": 100, "p95": 200, "p99": 300}
        mock_ut.return_value.get_current_uptime.return_value = 99.9
        mock_et.return_value.get_current_error_rate.return_value = 0.1
        
        with TestClient(app) as client:
            response = client.get("/status/metrics?project=live2&tier=hobby")
            assert response.status_code == 200
            data = response.json()
            assert data["project"] == "live2"
            assert data["tier"] == "hobby"
            assert "uptime" in data
            assert "response_time" in data
            assert "error_rate" in data


def test_sla_endpoint(app):
    """Test SLA compliance endpoint."""
    with patch("backend.monitoring.status.api.get_sla_calculator") as mock_calc:
        mock_calc.return_value.check_sla_compliance.return_value = {
            "status": "compliant",
            "month": "2025-12-01",
            "project": "live2",
            "tier": "hobby",
            "metrics": {
                "uptime": 99.9,
                "uptime_target": 95.0,
                "response_time_p95": 1000,
                "response_time_p95_target": 10000,
                "error_rate": 0.1,
                "error_rate_target": 1.0
            },
            "violations": [],
            "warnings": [],
            "compliant": True
        }
        
        with TestClient(app) as client:
            response = client.get("/status/sla?project=live2&tier=hobby")
            assert response.status_code == 200
            data = response.json()
            assert data["status"] == "compliant"
            assert data["compliant"] == True

