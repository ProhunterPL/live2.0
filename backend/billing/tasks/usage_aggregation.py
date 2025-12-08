"""
Usage aggregation task.

Periodically aggregates Redis usage data to PostgreSQL.
Can be run more frequently than monthly reset (e.g., daily).
"""

import logging
from datetime import date, timedelta
from sqlalchemy.orm import Session

from backend.billing.database import get_db_context
from backend.billing.models import User, Usage
from backend.api.v1.dependencies import get_redis_usage

logger = logging.getLogger(__name__)


def aggregate_usage_to_db():
    """
    Aggregate current month usage from Redis to PostgreSQL.
    
    This can be run daily to ensure data persistence.
    """
    logger.info("Starting usage aggregation...")
    
    with get_db_context() as db:
        # Get all active users
        users = db.query(User).filter(
            User.subscription_status.in_(["active", "trial"])
        ).all()
        
        redis_client = get_redis_usage()
        current_month = date.today().replace(day=1)
        
        aggregated_count = 0
        for user in users:
            try:
                # Get from Redis
                month_str = current_month.strftime("%Y-%m")
                reactions_key = f"usage:{user.id}:{month_str}:reactions"
                api_calls_key = f"usage:{user.id}:{month_str}:api_calls"
                
                reactions_count = int(redis_client.get(reactions_key) or 0)
                api_calls_count = int(redis_client.get(api_calls_key) or 0)
                
                # Update or create Usage record
                usage_record = db.query(Usage).filter(
                    Usage.user_id == user.id,
                    Usage.date == current_month
                ).first()
                
                if usage_record:
                    usage_record.reactions_count = reactions_count
                    usage_record.api_calls_count = api_calls_count
                    usage_record.tier = user.tier
                else:
                    usage_record = Usage(
                        user_id=user.id,
                        date=current_month,
                        reactions_count=reactions_count,
                        api_calls_count=api_calls_count,
                        tier=user.tier
                    )
                    db.add(usage_record)
                
                aggregated_count += 1
            except Exception as e:
                logger.error(f"Failed to aggregate usage for user {user.id}: {e}")
        
        db.commit()
        logger.info(f"Usage aggregation complete: {aggregated_count} users processed")

