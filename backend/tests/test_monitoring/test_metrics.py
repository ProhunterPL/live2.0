"""
Tests for monitoring metrics.
"""

import pytest
from datetime import datetime, date, timedelta
from unittest.mock import Mock, MagicMock

from backend.monitoring.metrics.response_time import ResponseTimeTracker
from backend.monitoring.metrics.uptime import UptimeTracker
from backend.monitoring.metrics.error_rate import ErrorRateTracker


@pytest.fixture
def mock_redis():
    """Create mock Redis client."""
    redis = MagicMock()
    redis.zadd = MagicMock()
    redis.zrange = MagicMock(return_value=[])
    redis.zrangebyscore = MagicMock(return_value=[])
    redis.zcard = MagicMock(return_value=0)
    redis.get = MagicMock(return_value=None)
    redis.expire = MagicMock()
    redis.ping = MagicMock()
    return redis


def test_response_time_tracker_record(mock_redis):
    """Test response time recording."""
    tracker = ResponseTimeTracker(redis_client=mock_redis)
    
    tracker.record_response_time("/api/test", "hobby", 150.5)
    
    mock_redis.zadd.assert_called_once()
    mock_redis.expire.assert_called_once()


def test_response_time_tracker_percentiles(mock_redis):
    """Test percentile calculation."""
    # Mock Redis to return sample data
    mock_redis.zrange.return_value = ["100", "150", "200", "250", "300"]
    
    tracker = ResponseTimeTracker(redis_client=mock_redis)
    percentiles = tracker.get_percentiles("/api/test", "hobby", [50, 95, 99])
    
    assert "p50" in percentiles
    assert "p95" in percentiles
    assert "p99" in percentiles


def test_uptime_tracker_record(mock_redis):
    """Test uptime tracking."""
    tracker = UptimeTracker(redis_client=mock_redis)
    
    tracker.record_health_check("hobby", True)
    
    mock_redis.zadd.assert_called_once()
    mock_redis.expire.assert_called_once()


def test_uptime_tracker_calculate(mock_redis):
    """Test uptime calculation."""
    # Mock Redis to return healthy checks
    mock_redis.zrangebyscore.return_value = [("1", 1000.0), ("1", 2000.0), ("1", 3000.0)]
    
    tracker = UptimeTracker(redis_client=mock_redis)
    uptime = tracker.calculate_uptime(
        "hobby",
        date.today() - timedelta(days=30),
        date.today()
    )
    
    assert uptime == 100.0  # All healthy


def test_error_rate_tracker_record(mock_redis):
    """Test error rate tracking."""
    tracker = ErrorRateTracker(redis_client=mock_redis)
    
    tracker.record_request("hobby", False)  # No error
    
    mock_redis.zadd.assert_called()
    mock_redis.expire.assert_called()


def test_error_rate_tracker_calculate(mock_redis):
    """Test error rate calculation."""
    # Mock Redis to return mixed results
    mock_redis.zrangebyscore.return_value = ["0", "0", "1", "0"]  # 1 error out of 4
    
    tracker = ErrorRateTracker(redis_client=mock_redis)
    error_rate = tracker.get_error_rate("hobby", date.today())
    
    assert error_rate == 25.0  # 25% error rate

