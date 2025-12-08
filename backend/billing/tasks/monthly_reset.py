"""
Monthly usage reset task.

Resets usage for all users on 1st of month.
Should be run as cron job or scheduled task.
"""

import logging
from datetime import date
from sqlalchemy.orm import Session

from backend.billing.database import get_db_context
from backend.billing.models import User
from backend.billing.usage_tracker import UsageTracker
from backend.api.v1.dependencies import get_redis_usage, get_rate_limiter

logger = logging.getLogger(__name__)


def reset_all_users_usage():
    """
    Reset usage for all users (called on 1st of month).
    
    Aggregates Redis data to PostgreSQL and clears Redis for new month.
    """
    logger.info("Starting monthly usage reset...")
    
    with get_db_context() as db:
        # Get all active users
        users = db.query(User).filter(
            User.subscription_status.in_(["active", "trial"])
        ).all()
        
        redis_client = get_redis_usage()
        rate_limiter = get_rate_limiter()
        usage_tracker = UsageTracker(db, redis_client, rate_limiter)
        
        reset_count = 0
        for user in users:
            try:
                usage_tracker.reset_monthly_usage(user.id)
                reset_count += 1
            except Exception as e:
                logger.error(f"Failed to reset usage for user {user.id}: {e}")
        
        logger.info(f"Monthly usage reset complete: {reset_count} users processed")

