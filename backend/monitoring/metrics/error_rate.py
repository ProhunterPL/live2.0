"""
Error rate tracking for API endpoints.
"""

import os
import logging
from typing import Dict, Optional, TYPE_CHECKING
from datetime import datetime, date, timedelta

if TYPE_CHECKING:
    import redis

try:
    import redis
except ImportError:
    redis = None  # type: ignore

logger = logging.getLogger(__name__)

from backend.monitoring.config import (
    REDIS_HOST, REDIS_PORT, REDIS_DB_METRICS,
    REDIS_USERNAME, REDIS_PASSWORD, METRICS_RETENTION_DAYS
)


class ErrorRateTracker:
    """
    Tracks error rates per tier.
    
    Calculates error rate as percentage of failed requests.
    """
    
    def __init__(self, redis_client: Optional["redis.Redis"] = None):
        """
        Initialize error rate tracker.
        
        Args:
            redis_client: Redis client for storage (if None, creates new connection)
        """
        if redis_client:
            self.redis = redis_client
        else:
            if redis is None:
                logger.debug("Redis module not available for error tracking (optional)")
                self.redis = None
            else:
                try:
                    redis_kwargs = {
                        "host": REDIS_HOST,
                        "port": REDIS_PORT,
                        "db": REDIS_DB_METRICS,
                        "decode_responses": True
                    }
                    if REDIS_USERNAME:
                        redis_kwargs["username"] = REDIS_USERNAME
                    if REDIS_PASSWORD:
                        redis_kwargs["password"] = REDIS_PASSWORD
                    
                    self.redis = redis.Redis(**redis_kwargs)
                    self.redis.ping()
                except Exception as e:
                    logger.warning(f"Redis not available for error tracking: {e}")
                    self.redis = None
    
    def record_request(
        self,
        tier: str,
        is_error: bool,
        timestamp: Optional[datetime] = None
    ):
        """
        Record request (success or error).
        
        Args:
            tier: User tier
            is_error: Whether request resulted in error
            timestamp: Optional timestamp (default: now)
        """
        if not self.redis:
            return
        
        if timestamp is None:
            timestamp = datetime.utcnow()
        
        try:
            # Store in Redis sorted set
            # Score = timestamp, value = "1" for error, "0" for success
            key = f"errors:{tier}"
            self.redis.zadd(
                key,
                {str(1 if is_error else 0): timestamp.timestamp()}
            )
            
            # Also track total requests
            total_key = f"requests:{tier}"
            self.redis.zadd(
                total_key,
                {str(1): timestamp.timestamp()}
            )
            
            # Set expiry
            self.redis.expire(key, 86400 * METRICS_RETENTION_DAYS)
            self.redis.expire(total_key, 86400 * METRICS_RETENTION_DAYS)
        except Exception as e:
            logger.warning(f"Failed to record request: {e}")
    
    def get_error_rate(
        self,
        tier: str,
        month: Optional[date] = None,
        window_minutes: Optional[int] = None
    ) -> float:
        """
        Get error rate for tier.
        
        Args:
            tier: User tier
            month: Month to check (if None and window_minutes is None, uses current month)
            window_minutes: Time window in minutes (if provided, overrides month)
        
        Returns:
            Error rate as percentage (0-100)
        """
        if not self.redis:
            return 0.0
        
        try:
            if window_minutes:
                # Use time window
                cutoff_time = (datetime.utcnow() - timedelta(minutes=window_minutes)).timestamp()
                
                errors_key = f"errors:{tier}"
                requests_key = f"requests:{tier}"
                
                errors = self.redis.zrangebyscore(
                    errors_key,
                    cutoff_time,
                    "+inf",
                    withscores=False
                )
                total_requests = self.redis.zrangebyscore(
                    requests_key,
                    cutoff_time,
                    "+inf",
                    withscores=False
                )
            else:
                # Use month
                if month is None:
                    month = date.today().replace(day=1)
                
                start_ts = datetime.combine(month, datetime.min.time()).timestamp()
                # End of month
                if month.month == 12:
                    end_date = date(month.year + 1, 1, 1) - timedelta(days=1)
                else:
                    end_date = date(month.year, month.month + 1, 1) - timedelta(days=1)
                end_ts = datetime.combine(end_date, datetime.max.time()).timestamp()
                
                errors_key = f"errors:{tier}"
                requests_key = f"requests:{tier}"
                
                errors = self.redis.zrangebyscore(
                    errors_key,
                    start_ts,
                    end_ts,
                    withscores=False
                )
                total_requests = self.redis.zrangebyscore(
                    requests_key,
                    start_ts,
                    end_ts,
                    withscores=False
                )
            
            # Count errors (value "1" = error)
            error_count = sum(1 for e in errors if e == "1")
            total_count = len(total_requests)
            
            if total_count == 0:
                return 0.0
            
            return (error_count / total_count) * 100.0
        except Exception as e:
            logger.warning(f"Failed to get error rate: {e}")
            return 0.0
    
    def get_current_error_rate(
        self,
        tier: str,
        window_minutes: int = 60
    ) -> float:
        """
        Get current error rate for last N minutes.
        
        Args:
            tier: User tier
            window_minutes: Number of minutes to look back (default: 60)
        
        Returns:
            Error rate as percentage
        """
        return self.get_error_rate(tier, window_minutes=window_minutes)
