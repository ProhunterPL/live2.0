"""
Configuration for monitoring module.
"""

import os
from typing import Optional

# Sentry Configuration
SENTRY_DSN: Optional[str] = os.getenv("SENTRY_DSN")
SENTRY_ENVIRONMENT: str = os.getenv("SENTRY_ENVIRONMENT", "production")

# UptimeRobot Configuration
UPTIMEROBOT_API_KEY: Optional[str] = os.getenv("UPTIMEROBOT_API_KEY")
UPTIMEROBOT_LEGALLY_MONITOR_ID: Optional[str] = os.getenv("UPTIMEROBOT_LEGALLY_MONITOR_ID")
UPTIMEROBOT_LIVE2_MONITOR_ID: Optional[str] = os.getenv("UPTIMEROBOT_LIVE2_MONITOR_ID")

# Redis Configuration (for metrics storage)
REDIS_HOST: str = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT: int = int(os.getenv("REDIS_PORT", 6379))
REDIS_DB_METRICS: int = int(os.getenv("REDIS_DB_METRICS", 2))  # Different DB for metrics
REDIS_USERNAME: Optional[str] = os.getenv("REDIS_USERNAME")
REDIS_PASSWORD: Optional[str] = os.getenv("REDIS_PASSWORD")

# Monitoring Settings
MONITORING_ENABLED: bool = os.getenv("MONITORING_ENABLED", "true").lower() == "true"
METRICS_RETENTION_DAYS: int = int(os.getenv("METRICS_RETENTION_DAYS", "30"))


def setup_sentry(dsn: Optional[str] = None, environment: Optional[str] = None):
    """
    Setup Sentry error tracking.
    
    Args:
        dsn: Sentry DSN (if None, uses SENTRY_DSN env var)
        environment: Environment name (if None, uses SENTRY_ENVIRONMENT env var)
    """
    if not dsn:
        dsn = SENTRY_DSN
    
    if not dsn:
        return  # Sentry not configured
    
    try:
        import sentry_sdk
        from sentry_sdk.integrations.fastapi import FastApiIntegration
        from sentry_sdk.integrations.sqlalchemy import SqlalchemyIntegration
        
        sentry_sdk.init(
            dsn=dsn,
            environment=environment or SENTRY_ENVIRONMENT,
            integrations=[
                FastApiIntegration(),
                SqlalchemyIntegration(),
            ],
            traces_sample_rate=0.1,  # 10% of transactions
            profiles_sample_rate=0.1,  # 10% of profiles
        )
    except ImportError:
        import logging
        logger = logging.getLogger(__name__)
        logger.warning("sentry-sdk not installed. Install with: pip install sentry-sdk[fastapi]")
