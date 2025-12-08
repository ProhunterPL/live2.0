"""
Usage tracking routes.
"""

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List
from datetime import date

from backend.billing.schemas import UsageResponse, UsageMonthResponse
from backend.billing.database import get_db
from backend.billing.usage_tracker import UsageTracker
from backend.billing.models import Usage as UsageModel
from backend.api.v1.dependencies import get_redis_usage, get_rate_limiter

router = APIRouter(prefix="/billing", tags=["billing"])


@router.get("/usage", response_model=UsageResponse)
async def get_usage(
    user_id: str,  # TODO: Get from JWT token or API key
    db: Session = Depends(get_db)
):
    """
    Get usage for current month and history.
    """
    from uuid import UUID
    
    try:
        user_uuid = UUID(user_id)
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid user ID")
    
    # Get usage tracker
    redis_client = get_redis_usage()
    rate_limiter = get_rate_limiter()
    usage_tracker = UsageTracker(db, redis_client, rate_limiter)
    
    # Get current month usage
    current_usage = usage_tracker.get_usage(user_uuid)
    current_month = date.today().replace(day=1)
    
    current_month_response = UsageMonthResponse(
        month=current_month.strftime("%Y-%m"),
        reactions=current_usage["reactions"],
        reactions_quota=current_usage["reactions_quota"],
        api_calls=current_usage["api_calls"],
        api_calls_quota=current_usage["api_calls_quota"],
        percentage_used=current_usage["percentage_used"]
    )
    
    # Get history (last 6 months)
    history = []
    for i in range(1, 7):  # Last 6 months
        month_date = (current_month.replace(day=1) - __import__("datetime").timedelta(days=30*i)).replace(day=1)
        month_usage = usage_tracker.get_usage(user_uuid, month_date)
        
        history.append(UsageMonthResponse(
            month=month_date.strftime("%Y-%m"),
            reactions=month_usage["reactions"],
            reactions_quota=month_usage["reactions_quota"],
            api_calls=month_usage["api_calls"],
            api_calls_quota=month_usage["api_calls_quota"],
            percentage_used=month_usage["percentage_used"]
        ))
    
    return UsageResponse(
        current_month=current_month_response,
        history=history
    )

