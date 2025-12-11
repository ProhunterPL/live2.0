"""
Rate limiter for API v1.

Tracks usage per user per tier and enforces quotas.
"""

import os
import logging
from typing import Dict, Optional, TYPE_CHECKING
from datetime import datetime, timedelta
from fastapi import HTTPException, status

if TYPE_CHECKING:
    import redis

try:
    import redis
except ImportError:
    redis = None  # type: ignore

from backend.api.v1.auth import User

logger = logging.getLogger(__name__)


class RateLimiter:
    """
    Rate limiter per tier with Redis backend.
    
    Tracks usage per user per month and enforces quotas.
    """
    
    TIER_QUOTAS = {
        "hobby": {
            "reactions_per_month": 10_000,
            "api_calls_per_month": 1_000
        },
        "research": {
            "reactions_per_month": 100_000,
            "api_calls_per_month": 10_000
        },
        "pro": {
            "reactions_per_month": 1_000_000,
            "api_calls_per_month": 100_000
        },
        "enterprise": {
            "reactions_per_month": float("inf"),
            "api_calls_per_month": float("inf")
        }
    }
    
    def __init__(self, redis_client: Optional["redis.Redis"] = None):
        """
        Initialize rate limiter.
        
        Args:
            redis_client: Redis client (if None, creates new connection)
        """
        if redis is None:
            logger.debug("Redis module not available (optional). Rate limiting will not work.")
            # Create a mock Redis client
            class MockRedis:
                def get(self, *args, **kwargs): return None
                def incrby(self, *args, **kwargs): pass
                def expire(self, *args, **kwargs): pass
            self.redis = MockRedis()  # type: ignore
        elif redis_client is None:
            from backend.api.v1.config import (
                REDIS_HOST, REDIS_PORT, REDIS_DB_USAGE,
                REDIS_USERNAME, REDIS_PASSWORD, REDIS_DECODE_RESPONSES
            )
            try:
                redis_kwargs = {
                    "host": REDIS_HOST,
                    "port": REDIS_PORT,
                    "db": REDIS_DB_USAGE,
                    "decode_responses": REDIS_DECODE_RESPONSES
                }
                # Add username/password only if provided (for Redis Labs)
                if REDIS_USERNAME:
                    redis_kwargs["username"] = REDIS_USERNAME
                if REDIS_PASSWORD:
                    redis_kwargs["password"] = REDIS_PASSWORD
                
                self.redis = redis.Redis(**redis_kwargs)
            except Exception as e:
                logger.warning(f"Failed to connect to Redis: {e}. Rate limiting will not work.")
                # Create a mock Redis client
                class MockRedis:
                    def get(self, *args, **kwargs): return None
                    def incrby(self, *args, **kwargs): pass
                    def expire(self, *args, **kwargs): pass
                self.redis = MockRedis()  # type: ignore
        else:
            self.redis = redis_client
    
    def check_quota(
        self, 
        user: User, 
        resource_type: str = "api_calls"
    ) -> Dict[str, any]:
        """
        Check if user has quota remaining.
        
        Args:
            user: User object
            resource_type: "api_calls" or "reactions"
        
        Returns:
            Dict with quota info:
            {
                "allowed": bool,
                "used": int,
                "quota": int,
                "remaining": int
            }
        
        Raises:
            HTTPException 429: Quota exceeded
        """
        tier_quota = self.TIER_QUOTAS.get(user.tier, self.TIER_QUOTAS["hobby"])
        quota = tier_quota.get(f"{resource_type}_per_month", 0)
        
        # Get current month usage
        current_month = datetime.now().strftime("%Y-%m")
        usage_key = f"usage:{user.id}:{current_month}:{resource_type}"
        
        try:
            used = int(self.redis.get(usage_key) or 0)
        except Exception as e:
            logger.warning(f"Failed to get usage from Redis: {e}")
            used = 0
        
        remaining = quota - used if quota != float("inf") else float("inf")
        
        result = {
            "allowed": remaining > 0 or quota == float("inf"),
            "used": used,
            "quota": quota if quota != float("inf") else None,
            "remaining": remaining if remaining != float("inf") else None
        }
        
        if not result["allowed"]:
            quota_str = "unlimited" if quota == float("inf") else str(quota)
            raise HTTPException(
                status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                detail=f"Quota exceeded. Used {used}/{quota_str} {resource_type} this month."
            )
        
        return result
    
    def track_usage(
        self, 
        user: User, 
        resource_type: str = "api_calls",
        amount: int = 1
    ):
        """
        Track usage for user.
        
        Args:
            user: User object
            resource_type: "api_calls" or "reactions"
            amount: Amount to increment (default: 1)
        """
        current_month = datetime.now().strftime("%Y-%m")
        usage_key = f"usage:{user.id}:{current_month}:{resource_type}"
        
        try:
            # Increment usage
            self.redis.incrby(usage_key, amount)
            
            # Set expiry to end of month (TTL)
            now = datetime.now()
            next_month = (now.replace(day=1) + timedelta(days=32)).replace(day=1)
            ttl_seconds = int((next_month - now).total_seconds())
            self.redis.expire(usage_key, ttl_seconds)
        except Exception as e:
            logger.error(f"Failed to track usage in Redis: {e}")
            # Don't fail the request if tracking fails

