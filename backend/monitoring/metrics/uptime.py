"""
Uptime tracking for system availability.
"""

import os
import logging
from typing import Dict, Optional
from datetime import datetime, timedelta, date
import redis

logger = logging.getLogger(__name__)

from backend.monitoring.config import (
    REDIS_HOST, REDIS_PORT, REDIS_DB_METRICS,
    REDIS_USERNAME, REDIS_PASSWORD, METRICS_RETENTION_DAYS
)


class UptimeTracker:
    """
    Tracks system uptime per tier.
    
    Calculates uptime percentage based on health check failures.
    """
    
    def __init__(self, redis_client: Optional[redis.Redis] = None):
        """
        Initialize uptime tracker.
        
        Args:
            redis_client: Redis client for storage (if None, creates new connection)
        """
        if redis_client:
            self.redis = redis_client
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
                logger.warning(f"Redis not available for uptime tracking: {e}")
                self.redis = None
    
    def record_health_check(
        self,
        tier: str,
        is_healthy: bool,
        timestamp: Optional[datetime] = None
    ):
        """
        Record health check result.
        
        Args:
            tier: User tier
            is_healthy: Whether system is healthy
            timestamp: Optional timestamp (default: now)
        """
        if not self.redis:
            return
        
        if timestamp is None:
            timestamp = datetime.utcnow()
        
        try:
            # Store in Redis sorted set (1 = healthy, 0 = unhealthy)
            key = f"health_check:{tier}"
            self.redis.zadd(
                key,
                {str(1 if is_healthy else 0): timestamp.timestamp()}
            )
            
            # Keep last 30 days
            cutoff = (timestamp - timedelta(days=30)).timestamp()
            self.redis.zremrangebyscore(key, "-inf", cutoff)
            
            # Set expiry
            self.redis.expire(key, 86400 * METRICS_RETENTION_DAYS)
        except Exception as e:
            logger.warning(f"Failed to record health check: {e}")
    
    def calculate_uptime(
        self,
        tier: str,
        start_date: date,
        end_date: date
    ) -> float:
        """
        Calculate uptime percentage for date range.
        
        Args:
            tier: User tier
            start_date: Start date
            end_date: End date
        
        Returns:
            Uptime percentage (0-100)
        """
        if not self.redis:
            return 100.0  # Assume healthy if no data
        
        try:
            key = f"health_check:{tier}"
            start_ts = datetime.combine(start_date, datetime.min.time()).timestamp()
            end_ts = datetime.combine(end_date, datetime.max.time()).timestamp()
            
            # Get all health checks in range
            checks = self.redis.zrangebyscore(
                key,
                start_ts,
                end_ts,
                withscores=True
            )
            
            if not checks:
                return 100.0  # Assume healthy if no data
            
            healthy_count = sum(1 for val, _ in checks if val == "1")
            total_count = len(checks)
            
            return (healthy_count / total_count) * 100.0 if total_count > 0 else 100.0
        except Exception as e:
            logger.warning(f"Failed to calculate uptime: {e}")
            return 100.0
    
    def get_current_uptime(
        self,
        tier: str,
        window_days: int = 30
    ) -> float:
        """
        Get current uptime for last N days.
        
        Args:
            tier: User tier
            window_days: Number of days to look back (default: 30)
        
        Returns:
            Uptime percentage
        """
        end_date = date.today()
        start_date = end_date - timedelta(days=window_days)
        return self.calculate_uptime(tier, start_date, end_date)
