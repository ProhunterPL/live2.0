"""
FastAPI middleware for monitoring.
"""

import time
import logging
from fastapi import Request, Response
from starlette.middleware.base import BaseHTTPMiddleware
from typing import Callable, Optional

from backend.monitoring.metrics.response_time import ResponseTimeTracker
from backend.monitoring.metrics.error_rate import ErrorRateTracker

logger = logging.getLogger(__name__)


def get_tier_from_request(request: Request) -> str:
    """
    Extract user tier from request.
    
    Tries multiple methods:
    1. Request state (if set by auth middleware)
    2. API key lookup (if X-API-Key header present)
    3. Endpoint path inference
    4. Default to "free"
    
    Args:
        request: FastAPI request
    
    Returns:
        Tier string (free, hobby, research, pro, enterprise, starter, professional, law_firm)
    """
    # Method 1: Try request state (set by auth middleware)
    try:
        if hasattr(request.state, "user"):
            user = request.state.user
            tier = getattr(user, "tier", None)
            if tier:
                return tier
    except:
        pass
    
    # Method 2: Try API key lookup
    try:
        api_key = request.headers.get("X-API-Key")
    except (KeyError, AttributeError):
        api_key = None
    
    if api_key:
        try:
            from backend.api.v1.auth import get_user_by_api_key
            user = get_user_by_api_key(api_key)
            if user and hasattr(user, "tier"):
                return user.tier
        except Exception as e:
            logger.debug(f"Failed to lookup user by API key: {e}")
    
    # Method 3: Try billing module direct lookup
    if api_key:
        try:
            from backend.billing.models import User
            from backend.billing.database import SessionLocal
            
            db = SessionLocal()
            try:
                db_user = db.query(User).filter(User.api_key == api_key).first()
                if db_user:
                    return db_user.tier
            finally:
                db.close()
        except Exception as e:
            logger.debug(f"Failed to lookup user from billing module: {e}")
    
    # Method 4: Infer from endpoint path (heuristic)
    path = request.url.path
    if "/api/v1/" in path:
        # Live 2.0 API - default to hobby for unauthenticated
        return "hobby"
    elif "/api/subscription" in path or "/api/auth" in path:
        # Legally endpoints - default to free for unauthenticated
        return "free"
    
    # Default
    return "free"


class ResponseTimeMiddleware(BaseHTTPMiddleware):
    """Middleware to track response times and errors."""
    
    def __init__(
        self,
        app,
        response_time_tracker: ResponseTimeTracker,
        error_rate_tracker: ErrorRateTracker
    ):
        """
        Initialize middleware.
        
        Args:
            app: FastAPI app
            response_time_tracker: Response time tracker instance
            error_rate_tracker: Error rate tracker instance
        """
        super().__init__(app)
        self.tracker = response_time_tracker
        self.error_tracker = error_rate_tracker
    
    async def dispatch(self, request: Request, call_next: Callable) -> Response:
        """
        Process request and track metrics.
        
        Args:
            request: FastAPI request
            call_next: Next middleware/handler
        
        Returns:
            Response
        """
        start_time = time.time()
        
        # Get user tier from request
        tier = get_tier_from_request(request)
        
        endpoint = request.url.path
        
        # Process request
        try:
            response = await call_next(request)
            is_error = response.status_code >= 400
        except Exception as e:
            logger.error(f"Request failed: {e}")
            is_error = True
            raise
        finally:
            # Calculate response time
            response_time_ms = (time.time() - start_time) * 1000
            
            # Record response time
            try:
                self.tracker.record_response_time(endpoint, tier, response_time_ms)
            except Exception as e:
                logger.warning(f"Failed to record response time: {e}")
            
            # Record error
            try:
                self.error_tracker.record_request(tier, is_error)
            except Exception as e:
                logger.warning(f"Failed to record error: {e}")
        
        return response
