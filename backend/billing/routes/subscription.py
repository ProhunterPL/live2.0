"""
Subscription management routes.
"""

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from backend.billing.schemas import UpgradeRequest, CancelRequest, SubscriptionResponse
from backend.billing.database import get_db
from backend.billing.subscriptions import SubscriptionManager
from backend.billing.payments import PaymentProcessor
from backend.billing.usage_tracker import UsageTracker
from backend.billing.config import STRIPE_SECRET_KEY
from backend.billing.dependencies import get_current_user
from backend.billing.models import User
from backend.api.v1.dependencies import get_redis_usage, get_rate_limiter

router = APIRouter(prefix="/billing", tags=["billing"])


def get_subscription_manager(db: Session = Depends(get_db)) -> SubscriptionManager:
    """Get subscription manager instance."""
    payment_processor = PaymentProcessor(STRIPE_SECRET_KEY) if STRIPE_SECRET_KEY else None
    return SubscriptionManager(db, payment_processor)


@router.get("/subscription", response_model=SubscriptionResponse)
async def get_subscription(
    db: Session = Depends(get_db),
    subscription_manager: SubscriptionManager = Depends(get_subscription_manager),
    user: User = Depends(get_current_user),
):
    """
    Get current subscription.
    
    Returns subscription info with usage.
    """
    subscription = (
        subscription_manager.get_active_subscription(user.id)
        or subscription_manager.get_current_subscription(user.id)
    )
    if not subscription:
        raise HTTPException(status_code=404, detail="No subscription found")
    
    # Get usage
    redis_client = get_redis_usage()
    rate_limiter = get_rate_limiter()
    usage_tracker = UsageTracker(db, redis_client, rate_limiter)
    usage = usage_tracker.get_usage(user.id)
    
    return SubscriptionResponse(
        tier=subscription.tier,
        status=subscription.status,
        current_period_start=subscription.current_period_start,
        current_period_end=subscription.current_period_end,
        usage=usage
    )


@router.post("/upgrade")
async def upgrade_subscription(
    request: UpgradeRequest,
    db: Session = Depends(get_db),
    subscription_manager: SubscriptionManager = Depends(get_subscription_manager),
    user: User = Depends(get_current_user),
):
    """
    Upgrade/downgrade subscription tier.
    """
    try:
        subscription = subscription_manager.update_tier(user.id, request.tier)
        return {
            "success": True,
            "subscription": {
                "tier": subscription.tier,
                "status": subscription.status
            }
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/cancel")
async def cancel_subscription(
    request: CancelRequest,
    db: Session = Depends(get_db),
    subscription_manager: SubscriptionManager = Depends(get_subscription_manager),
    user: User = Depends(get_current_user),
):
    """
    Cancel subscription.
    """
    try:
        subscription = subscription_manager.cancel_subscription(
            user.id,
            cancel_at_period_end=request.cancel_at_period_end
        )
        return {
            "success": True,
            "subscription": {
                "status": subscription.status,
                "cancel_at_period_end": subscription.cancel_at_period_end
            }
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

