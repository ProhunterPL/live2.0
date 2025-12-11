"""
Response time tracking for API endpoints.
"""

import os
import logging
from typing import Dict, List, Optional
from datetime import datetime, timedelta
import redis

logger = logging.getLogger(__name__)

from backend.monitoring.config import (
    REDIS_HOST, REDIS_PORT, REDIS_DB_METRICS,
    REDIS_USERNAME, REDIS_PASSWORD, METRICS_RETENTION_DAYS
)


class ResponseTimeTracker:
    """
    Tracks response times per endpoint per tier.
    
    Calculates p50, p95, p99 percentiles.
    """
    
    def __init__(
        self,
        redis_client: Optional[redis.Redis] = None,
        window_size: int = 1000  # Number of samples per window
    ):
        """
        Initialize response time tracker.
        
        Args:
            redis_client: Redis client for storage (if None, creates new connection)
            window_size: Number of samples per percentile calculation window
        """
        self.window_size = window_size
        
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
                logger.warning(f"Redis not available for metrics: {e}")
                self.redis = None
    
    def record_response_time(
        self,
        endpoint: str,
        tier: str,
        response_time_ms: float,
        timestamp: Optional[datetime] = None
    ):
        """
        Record response time for endpoint.
        
        Args:
            endpoint: API endpoint path (e.g., "/api/v1/molecules")
            tier: User tier (hobby, research, pro, enterprise, free, starter, professional, law_firm)
            response_time_ms: Response time in milliseconds
            timestamp: Optional timestamp (default: now)
        """
        if not self.redis:
            return  # Redis not available
        
        if timestamp is None:
            timestamp = datetime.utcnow()
        
        try:
            # Store in Redis sorted set (score = timestamp, value = response_time)
            key = f"response_time:{endpoint}:{tier}"
            self.redis.zadd(
                key,
                {str(response_time_ms): timestamp.timestamp()}
            )
            
            # Keep only last window_size samples
            count = self.redis.zcard(key)
            if count > self.window_size:
                # Remove oldest samples
                self.redis.zremrangebyrank(key, 0, count - self.window_size - 1)
            
            # Set expiry (retention days)
            self.redis.expire(key, 86400 * METRICS_RETENTION_DAYS)
        except Exception as e:
            logger.warning(f"Failed to record response time: {e}")
    
    def get_percentiles(
        self,
        endpoint: str,
        tier: str,
        percentiles: List[float] = [50, 95, 99]
    ) -> Dict[str, float]:
        """
        Get percentile response times.
        
        Args:
            endpoint: API endpoint path
            tier: User tier
            percentiles: List of percentiles to calculate (default: [50, 95, 99])
        
        Returns:
            Dict mapping percentile to response time in ms
        """
        if not self.redis:
            return {}
        
        try:
            key = f"response_time:{endpoint}:{tier}"
            
            # Get all samples
            samples = self.redis.zrange(key, 0, -1, withscores=False)
            if not samples:
                return {}
            
            # Convert to float and sort
            response_times = sorted([float(s) for s in samples])
            
            if not response_times:
                return {}
            
            # Calculate percentiles
            result = {}
            for p in percentiles:
                index = int(len(response_times) * (p / 100.0))
                if index >= len(response_times):
                    index = len(response_times) - 1
                result[f"p{p}"] = response_times[index]
            
            return result
        except Exception as e:
            logger.warning(f"Failed to get percentiles: {e}")
            return {}
    
    def get_average(
        self,
        endpoint: str,
        tier: str,
        window_minutes: int = 60
    ) -> float:
        """
        Get average response time for time window.
        
        Args:
            endpoint: API endpoint path
            tier: User tier
            window_minutes: Time window in minutes (default: 60)
        
        Returns:
            Average response time in ms
        """
        if not self.redis:
            return 0.0
        
        try:
            key = f"response_time:{endpoint}:{tier}"
            cutoff_time = (datetime.utcnow() - timedelta(minutes=window_minutes)).timestamp()
            
            # Get samples from time window
            samples = self.redis.zrangebyscore(
                key,
                cutoff_time,
                "+inf",
                withscores=False
            )
            
            if not samples:
                return 0.0
            
            response_times = [float(s) for s in samples]
            return sum(response_times) / len(response_times)
        except Exception as e:
            logger.warning(f"Failed to get average: {e}")
            return 0.0
