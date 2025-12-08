"""
API key authentication for API v1.
"""

from typing import Optional
from fastapi import Request, HTTPException, status
from pydantic import BaseModel

logger = __import__("logging").getLogger(__name__)


class User(BaseModel):
    """
    User model (placeholder - will be replaced with billing module in 1.3).
    
    For now, this is a simple model for MVP.
    """
    id: str
    api_key: str
    tier: str  # "hobby", "research", "pro", "enterprise"
    subscription_status: str  # "active", "cancelled", "expired"


async def verify_api_key(request: Request) -> User:
    """
    Verify API key from X-API-Key header.
    
    Args:
        request: FastAPI request object
    
    Returns:
        User object with tier info
    
    Raises:
        HTTPException 401: Invalid or missing API key
    """
    api_key = request.headers.get("X-API-Key")
    if not api_key:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing API key. Provide X-API-Key header."
        )
    
    user = get_user_by_api_key(api_key)
    if not user or user.subscription_status != "active":
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or expired API key"
        )
    
    return user


def get_user_by_api_key(api_key: str) -> Optional[User]:
    """
    Get user by API key.
    
    Uses billing module database lookup if available, otherwise falls back to placeholder.
    
    Args:
        api_key: API key string
    
    Returns:
        User object or None
    """
    # Try billing module first (if available)
    try:
        from backend.billing.models import get_user_by_api_key as billing_get_user
        from backend.billing.database import SessionLocal
        
        db = SessionLocal()
        try:
            db_user = billing_get_user(api_key, db_session=db)
            if db_user:
                # Convert SQLAlchemy User to Pydantic User
                return User(
                    id=str(db_user.id),
                    api_key=db_user.api_key,
                    tier=db_user.tier,
                    subscription_status=db_user.subscription_status
                )
        finally:
            db.close()
    except ImportError:
        # Billing module not available yet
        logger.debug("Billing module not available, using placeholder")
    except Exception as e:
        logger.warning(f"Failed to lookup user from billing module: {e}")
    
    # Fallback: Placeholder implementation (for development)
    import redis
    from backend.api.v1.config import (
        REDIS_HOST, REDIS_PORT, REDIS_USERNAME, REDIS_PASSWORD, REDIS_DECODE_RESPONSES
    )
    
    try:
        redis_kwargs = {
            "host": REDIS_HOST,
            "port": REDIS_PORT,
            "db": 2,  # Separate DB for user data (temporary)
            "decode_responses": REDIS_DECODE_RESPONSES
        }
        # Add username/password only if provided (for Redis Labs)
        if REDIS_USERNAME:
            redis_kwargs["username"] = REDIS_USERNAME
        if REDIS_PASSWORD:
            redis_kwargs["password"] = REDIS_PASSWORD
        
        redis_client = redis.Redis(**redis_kwargs)
        
        user_data = redis_client.get(f"user:api_key:{api_key}")
        if user_data:
            import json
            data = json.loads(user_data)
            return User(**data)
    except Exception as e:
        logger.warning(f"Failed to lookup user by API key: {e}")
    
    # Fallback: Return None (will trigger 401)
    return None

