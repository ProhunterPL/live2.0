"""
Tests for usage tracking.
"""

import pytest
from sqlalchemy.orm import Session
from unittest.mock import Mock, MagicMock
from datetime import date, datetime, timedelta

from backend.billing.usage_tracker import UsageTracker
from backend.billing.models import User, Usage
from backend.billing.auth import create_user
# Use conftest db fixture
from backend.api.v1.rate_limiter import RateLimiter


# Use conftest db fixture


@pytest.fixture
def mock_redis():
    """Create mock Redis client."""
    redis_client = Mock()
    redis_client.get = Mock(return_value="0")
    redis_client.incrby = Mock(return_value=1)
    redis_client.expire = Mock(return_value=True)
    redis_client.delete = Mock(return_value=1)
    return redis_client


@pytest.fixture
def mock_rate_limiter():
    """Create mock rate limiter."""
    # Create a real RateLimiter instance but with mocked Redis
    from backend.api.v1.rate_limiter import RateLimiter
    from unittest.mock import Mock
    
    limiter = RateLimiter(redis_client=Mock())
    # Mock check_quota method
    limiter.check_quota = Mock(return_value={
        "allowed": True,
        "used": 0,
        "quota": 100000,
        "remaining": 100000
    })
    return limiter


@pytest.fixture
def test_user(db: Session):
    """Create test user."""
    import uuid
    unique_email = f"test_usage_{uuid.uuid4().hex[:8]}@example.com"
    return create_user(
        db=db,
        email=unique_email,
        password="password",
        tier="research"
    )


def test_track_reaction(db: Session, mock_redis, mock_rate_limiter, test_user):
    """Test reaction tracking."""
    tracker = UsageTracker(db, mock_redis, mock_rate_limiter)
    
    result = tracker.track_reaction(test_user.id)
    
    assert result is True
    mock_redis.incrby.assert_called_once()
    mock_redis.expire.assert_called_once()


def test_track_reaction_quota_exceeded(db: Session, mock_redis, mock_rate_limiter, test_user):
    """Test reaction tracking when quota exceeded."""
    # Create new tracker with mocked rate limiter that returns quota exceeded
    from backend.billing.usage_tracker import UsageTracker
    from backend.api.v1.rate_limiter import RateLimiter
    from unittest.mock import Mock
    
    # Create a mock rate limiter that returns quota exceeded
    mock_limiter = Mock(spec=RateLimiter)
    mock_limiter.check_quota = Mock(return_value={
        "allowed": False,
        "used": 100000,
        "quota": 100000,
        "remaining": 0
    })
    
    tracker = UsageTracker(db, mock_redis, mock_limiter)
    
    result = tracker.track_reaction(test_user.id)
    
    assert result is False
    mock_redis.incrby.assert_not_called()


def test_track_api_call(db: Session, mock_redis, mock_rate_limiter, test_user):
    """Test API call tracking."""
    tracker = UsageTracker(db, mock_redis, mock_rate_limiter)
    
    result = tracker.track_api_call(test_user.id, "molecules")
    
    assert result is True
    mock_redis.incrby.assert_called_once()


def test_get_usage(db: Session, mock_redis, mock_rate_limiter, test_user):
    """Test getting usage."""
    # Mock Redis to return usage
    mock_redis.get = Mock(side_effect=lambda key: {
        f"usage:{test_user.id}:2025-12:reactions": "5000",
        f"usage:{test_user.id}:2025-12:api_calls": "1200"
    }.get(key, "0"))
    
    tracker = UsageTracker(db, mock_redis, mock_rate_limiter)
    
    usage = tracker.get_usage(test_user.id)
    
    assert usage["reactions"] == 5000
    assert usage["api_calls"] == 1200
    assert usage["reactions_quota"] == 100000
    assert usage["percentage_used"] > 0


def test_check_quota(db: Session, mock_redis, mock_rate_limiter, test_user):
    """Test quota checking."""
    from datetime import datetime
    
    current_month = datetime.now().strftime("%Y-%m")
    mock_redis.get = Mock(side_effect=lambda key: {
        f"usage:{test_user.id}:{current_month}:reactions": "5000",
        f"usage:{test_user.id}:{current_month}:api_calls": "1200"
    }.get(key, "0"))
    
    tracker = UsageTracker(db, mock_redis, mock_rate_limiter)
    
    quota = tracker.check_quota(test_user.id)
    
    assert quota["reactions_remaining"] == 95000
    assert quota["api_calls_remaining"] is not None


def test_reset_monthly_usage(db: Session, mock_redis, mock_rate_limiter, test_user):
    """Test monthly usage reset."""
    from datetime import date, timedelta
    
    # Calculate previous month (always use last month for test)
    today = date.today()
    # Get first day of current month, then subtract 1 day to get last day of previous month
    first_of_current = today.replace(day=1)
    last_of_prev = first_of_current - timedelta(days=1)
    prev_month = last_of_prev.replace(day=1)  # First day of previous month
    
    prev_month_str = prev_month.strftime("%Y-%m")
    
    mock_redis.get = Mock(side_effect=lambda key: {
        f"usage:{test_user.id}:{prev_month_str}:reactions": "10000",
        f"usage:{test_user.id}:{prev_month_str}:api_calls": "5000"
    }.get(key, "0"))
    
    tracker = UsageTracker(db, mock_redis, mock_rate_limiter)
    
    tracker.reset_monthly_usage(test_user.id)
    
    # Check Usage record was created
    usage_record = db.query(Usage).filter(
        Usage.user_id == test_user.id,
        Usage.date == prev_month
    ).first()
    
    assert usage_record is not None
    assert usage_record.reactions_count == 10000
    assert usage_record.api_calls_count == 5000

