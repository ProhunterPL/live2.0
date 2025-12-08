"""
FastAPI dependencies for API v1.
"""

from functools import lru_cache
import redis
from typing import Optional

from backend.api.v1.config import (
    REDIS_HOST, REDIS_PORT, REDIS_DB_JOBS, REDIS_DB_USAGE,
    REDIS_USERNAME, REDIS_PASSWORD, REDIS_DECODE_RESPONSES
)
from backend.api.v1.rate_limiter import RateLimiter
from backend.api.v1.jobs import JobProcessor
from backend.api.v1.storage import StorageManager

# Global instances (singleton pattern)
_redis_jobs: Optional[redis.Redis] = None
_redis_usage: Optional[redis.Redis] = None
_rate_limiter: Optional[RateLimiter] = None
_job_processor: Optional[JobProcessor] = None
_storage_manager: Optional[StorageManager] = None


def get_redis_jobs() -> redis.Redis:
    """Get Redis client for jobs (singleton)."""
    global _redis_jobs
    if _redis_jobs is None:
        try:
            redis_kwargs = {
                "host": REDIS_HOST,
                "port": REDIS_PORT,
                "db": REDIS_DB_JOBS,
                "decode_responses": REDIS_DECODE_RESPONSES
            }
            # Add username/password only if provided (for Redis Labs)
            if REDIS_USERNAME:
                redis_kwargs["username"] = REDIS_USERNAME
            if REDIS_PASSWORD:
                redis_kwargs["password"] = REDIS_PASSWORD
            
            _redis_jobs = redis.Redis(**redis_kwargs)
            # Test connection
            _redis_jobs.ping()
        except Exception as e:
            import logging
            logger = logging.getLogger(__name__)
            logger.warning(f"Redis not available: {e}. API v1 will work but jobs and rate limiting may fail.")
            # Return a mock Redis client that will fail gracefully
            class MockRedis:
                def get(self, *args, **kwargs): return None
                def setex(self, *args, **kwargs): pass
                def lpush(self, *args, **kwargs): pass
                def brpop(self, *args, **kwargs): return None
                def incrby(self, *args, **kwargs): pass
                def expire(self, *args, **kwargs): pass
                def ping(self): raise Exception("Redis not available")
            _redis_jobs = MockRedis()
    return _redis_jobs


def get_redis_usage() -> redis.Redis:
    """Get Redis client for usage tracking (singleton)."""
    global _redis_usage
    if _redis_usage is None:
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
            
            _redis_usage = redis.Redis(**redis_kwargs)
            # Test connection
            _redis_usage.ping()
        except Exception as e:
            import logging
            logger = logging.getLogger(__name__)
            logger.warning(f"Redis not available: {e}. API v1 will work but rate limiting may fail.")
            # Return a mock Redis client that will fail gracefully
            class MockRedis:
                def get(self, *args, **kwargs): return None
                def setex(self, *args, **kwargs): pass
                def incrby(self, *args, **kwargs): pass
                def expire(self, *args, **kwargs): pass
                def ping(self): raise Exception("Redis not available")
            _redis_usage = MockRedis()
    return _redis_usage


def get_rate_limiter() -> RateLimiter:
    """Get rate limiter instance (singleton)."""
    global _rate_limiter
    if _rate_limiter is None:
        _rate_limiter = RateLimiter(redis_client=get_redis_usage())
    return _rate_limiter


def get_job_processor() -> JobProcessor:
    """Get job processor instance (singleton)."""
    global _job_processor
    if _job_processor is None:
        from backend.api.v1.config import BASE_RESULTS_DIR
        _job_processor = JobProcessor(
            redis_client=get_redis_jobs(),
            base_results_dir=BASE_RESULTS_DIR,
            storage_manager=get_storage_manager()
        )
    return _job_processor


def get_storage_manager() -> StorageManager:
    """Get storage manager instance (singleton)."""
    global _storage_manager
    if _storage_manager is None:
        from backend.api.v1.config import (
            STORAGE_TYPE, STORAGE_BASE_PATH, S3_BUCKET, S3_REGION,
            S3_ACCESS_KEY_ID, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
        )
        _storage_manager = StorageManager(
            storage_type=STORAGE_TYPE,
            base_path=str(STORAGE_BASE_PATH),
            s3_bucket=S3_BUCKET,
            s3_region=S3_REGION,
            s3_access_key_id=S3_ACCESS_KEY_ID,
            s3_secret_access_key=S3_SECRET_ACCESS_KEY,
            s3_endpoint_url=S3_ENDPOINT_URL
        )
    return _storage_manager

