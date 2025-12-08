"""
Usage tracking for billing module.

Tracks usage per user per month using Redis (real-time) and PostgreSQL (historical).
"""

from typing import Dict, Optional
from datetime import datetime, date, timedelta
from sqlalchemy.orm import Session
import redis
import uuid
import logging

from backend.billing.models import User, Usage
from backend.api.v1.rate_limiter import RateLimiter

logger = logging.getLogger(__name__)


class UsageTracker:
    """
    Tracks usage per user per month.
    
    Uses Redis for real-time tracking and PostgreSQL for historical data.
    """
    
    def __init__(
        self, 
        db: Session,
        redis_client: redis.Redis,
        rate_limiter: RateLimiter
    ):
        """
        Initialize usage tracker.
        
        Args:
            db: SQLAlchemy database session
            redis_client: Redis client
            rate_limiter: RateLimiter instance (for quota info)
        """
        self.db = db
        self.redis = redis_client
        self.rate_limiter = rate_limiter
    
    def track_reaction(self, user_id: uuid.UUID) -> bool:
        """
        Track one reaction, return True if within quota.
        
        Args:
            user_id: User ID
        
        Returns:
            True if within quota, False if exceeded
        """
        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return False
        
        # Check quota using rate limiter
        try:
            quota_info = self.rate_limiter.check_quota(user, resource_type="reactions")
            if not quota_info["allowed"]:
                return False
        except Exception as e:
            logger.warning(f"Rate limiter check failed: {e}")
            # If rate limiter fails, allow (graceful degradation)
            pass
        
        # Track in Redis
        current_month = datetime.now().strftime("%Y-%m")
        usage_key = f"usage:{user_id}:{current_month}:reactions"
        
        try:
            self.redis.incrby(usage_key, 1)
            # Set TTL to end of month
            now = datetime.now()
            next_month = (now.replace(day=1) + timedelta(days=32)).replace(day=1)
            ttl_seconds = int((next_month - now).total_seconds())
            self.redis.expire(usage_key, ttl_seconds)
        except Exception as e:
            logger.warning(f"Failed to track usage in Redis: {e}")
            # Don't fail the request if tracking fails
        
        return True
    
    def track_api_call(self, user_id: uuid.UUID, endpoint: str) -> bool:
        """
        Track API call, return True if within quota.
        
        Args:
            user_id: User ID
            endpoint: Endpoint name
        
        Returns:
            True if within quota, False if exceeded
        """
        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return False
        
        # Check quota
        try:
            quota_info = self.rate_limiter.check_quota(user, resource_type="api_calls")
            if not quota_info["allowed"]:
                return False
        except Exception as e:
            logger.warning(f"Rate limiter check failed: {e}")
            pass
        
        # Track in Redis
        current_month = datetime.now().strftime("%Y-%m")
        usage_key = f"usage:{user_id}:{current_month}:api_calls"
        
        try:
            self.redis.incrby(usage_key, 1)
            # Set TTL
            now = datetime.now()
            next_month = (now.replace(day=1) + timedelta(days=32)).replace(day=1)
            ttl_seconds = int((next_month - now).total_seconds())
            self.redis.expire(usage_key, ttl_seconds)
        except Exception as e:
            logger.warning(f"Failed to track usage in Redis: {e}")
        
        return True
    
    def get_usage(
        self, 
        user_id: uuid.UUID, 
        month: Optional[date] = None
    ) -> Dict:
        """
        Get usage for month (default: current month).
        
        Args:
            user_id: User ID
            month: Month date (default: current month)
        
        Returns:
            Dict with usage data
        """
        if month is None:
            month = date.today().replace(day=1)
        
        # Get from Redis (real-time)
        current_month = month.strftime("%Y-%m")
        reactions_key = f"usage:{user_id}:{current_month}:reactions"
        api_calls_key = f"usage:{user_id}:{current_month}:api_calls"
        
        try:
            reactions_count = int(self.redis.get(reactions_key) or 0)
            api_calls_count = int(self.redis.get(api_calls_key) or 0)
        except Exception as e:
            logger.warning(f"Failed to get usage from Redis: {e}")
            reactions_count = 0
            api_calls_count = 0
        
        # Get from PostgreSQL (historical, if exists)
        usage_record = self.db.query(Usage).filter(
            Usage.user_id == user_id,
            Usage.date == month
        ).first()
        
        # Merge data (Redis takes precedence for current month)
        if usage_record and month < date.today().replace(day=1):
            # Historical month - use PostgreSQL
            reactions_count = usage_record.reactions_count
            api_calls_count = usage_record.api_calls_count
        
        # Get quota from user tier
        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return {
                "reactions": 0,
                "reactions_quota": None,
                "api_calls": 0,
                "api_calls_quota": None,
                "percentage_used": 0.0
            }
        
        tier_quota = RateLimiter.TIER_QUOTAS.get(user.tier, RateLimiter.TIER_QUOTAS["hobby"])
        reactions_quota = tier_quota.get("reactions_per_month")
        api_calls_quota = tier_quota.get("api_calls_per_month")
        
        # Calculate percentage
        if reactions_quota and reactions_quota != float("inf"):
            percentage_used = (reactions_count / reactions_quota) * 100
        else:
            percentage_used = 0.0
        
        return {
            "reactions": reactions_count,
            "reactions_quota": reactions_quota if reactions_quota != float("inf") else None,
            "api_calls": api_calls_count,
            "api_calls_quota": api_calls_quota if api_calls_quota != float("inf") else None,
            "percentage_used": percentage_used
        }
    
    def check_quota(self, user_id: uuid.UUID) -> Dict:
        """
        Check remaining quota.
        
        Args:
            user_id: User ID
        
        Returns:
            Dict with quota info
        """
        usage = self.get_usage(user_id)
        return {
            "reactions_remaining": (
                usage["reactions_quota"] - usage["reactions"]
                if usage["reactions_quota"] is not None else None
            ),
            "api_calls_remaining": (
                usage["api_calls_quota"] - usage["api_calls"]
                if usage["api_calls_quota"] is not None else None
            ),
            **usage
        }
    
    def reset_monthly_usage(self, user_id: uuid.UUID):
        """
        Reset usage for new month (called on 1st of month).
        
        Aggregates Redis data to PostgreSQL and clears Redis.
        
        Args:
            user_id: User ID
        """
        # Get previous month usage from Redis
        prev_month = (date.today().replace(day=1) - timedelta(days=1)).replace(day=1)
        prev_month_str = prev_month.strftime("%Y-%m")
        
        reactions_key = f"usage:{user_id}:{prev_month_str}:reactions"
        api_calls_key = f"usage:{user_id}:{prev_month_str}:api_calls"
        
        try:
            reactions_count = int(self.redis.get(reactions_key) or 0)
            api_calls_count = int(self.redis.get(api_calls_key) or 0)
        except Exception as e:
            logger.warning(f"Failed to get usage from Redis for reset: {e}")
            reactions_count = 0
            api_calls_count = 0
        
        # Get user tier
        user = self.db.query(User).filter(User.id == user_id).first()
        tier = user.tier if user else "hobby"
        
        # Save to PostgreSQL
        usage_record = self.db.query(Usage).filter(
            Usage.user_id == user_id,
            Usage.date == prev_month
        ).first()
        
        if usage_record:
            usage_record.reactions_count = reactions_count
            usage_record.api_calls_count = api_calls_count
            usage_record.tier = tier
        else:
            usage_record = Usage(
                user_id=user_id,
                date=prev_month,
                reactions_count=reactions_count,
                api_calls_count=api_calls_count,
                tier=tier
            )
            self.db.add(usage_record)
        
        self.db.commit()
        
        # Clear Redis keys (they will expire anyway, but clear explicitly)
        try:
            self.redis.delete(reactions_key)
            self.redis.delete(api_calls_key)
        except Exception as e:
            logger.warning(f"Failed to clear Redis keys: {e}")

