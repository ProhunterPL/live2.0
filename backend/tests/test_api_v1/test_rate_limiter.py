"""
Tests for API v1 rate limiter.
"""

import pytest
from unittest.mock import Mock, MagicMock
from datetime import datetime

from backend.api.v1.rate_limiter import RateLimiter
from backend.api.v1.auth import User


def test_rate_limiter_tier_quotas():
    """Test tier quotas are correctly defined."""
    limiter = RateLimiter(redis_client=Mock())
    
    assert "hobby" in limiter.TIER_QUOTAS
    assert "research" in limiter.TIER_QUOTAS
    assert "pro" in limiter.TIER_QUOTAS
    assert "enterprise" in limiter.TIER_QUOTAS
    
    assert limiter.TIER_QUOTAS["hobby"]["reactions_per_month"] == 10_000
    assert limiter.TIER_QUOTAS["enterprise"]["reactions_per_month"] == float("inf")


def test_check_quota_allowed():
    """Test check_quota returns allowed=True when quota available."""
    mock_redis = Mock()
    mock_redis.get.return_value = "5000"  # Used 5000
    
    limiter = RateLimiter(redis_client=mock_redis)
    user = User(
        id="user_123",
        api_key="test_key",
        tier="research",
        subscription_status="active"
    )
    
    result = limiter.check_quota(user, resource_type="reactions")
    
    assert result["allowed"] is True
    assert result["used"] == 5000
    assert result["quota"] == 100_000
    assert result["remaining"] == 95_000


def test_check_quota_exceeded():
    """Test check_quota raises 429 when quota exceeded."""
    from fastapi import HTTPException
    
    mock_redis = Mock()
    mock_redis.get.return_value = "100000"  # Used all quota
    
    limiter = RateLimiter(redis_client=mock_redis)
    user = User(
        id="user_123",
        api_key="test_key",
        tier="research",
        subscription_status="active"
    )
    
    with pytest.raises(HTTPException) as exc_info:
        limiter.check_quota(user, resource_type="reactions")
    
    assert exc_info.value.status_code == 429
    assert "Quota exceeded" in exc_info.value.detail


def test_check_quota_enterprise_unlimited():
    """Test enterprise tier has unlimited quota."""
    mock_redis = Mock()
    mock_redis.get.return_value = "1000000"  # Used 1M
    
    limiter = RateLimiter(redis_client=mock_redis)
    user = User(
        id="user_123",
        api_key="test_key",
        tier="enterprise",
        subscription_status="active"
    )
    
    result = limiter.check_quota(user, resource_type="reactions")
    
    assert result["allowed"] is True
    assert result["quota"] is None  # Unlimited
    assert result["remaining"] is None  # Unlimited


def test_track_usage():
    """Test track_usage increments usage in Redis."""
    mock_redis = Mock()
    mock_redis.get.return_value = "100"
    
    limiter = RateLimiter(redis_client=mock_redis)
    user = User(
        id="user_123",
        api_key="test_key",
        tier="research",
        subscription_status="active"
    )
    
    limiter.track_usage(user, resource_type="api_calls", amount=5)
    
    # Verify incrby was called
    mock_redis.incrby.assert_called_once()
    # Verify expire was called (for TTL)
    mock_redis.expire.assert_called_once()

