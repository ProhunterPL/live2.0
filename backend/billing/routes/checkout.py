"""
Stripe Checkout / Customer Portal routes for billing module.

MVP flow:
- user authenticates (JWT or API key)
- create checkout session for a tier (subscription mode)
- Stripe webhooks update subscription status in our DB
"""

from __future__ import annotations

import logging
from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel
from sqlalchemy.orm import Session

from backend.billing.config import (
    STRIPE_SECRET_KEY,
    STRIPE_PUBLISHABLE_KEY,
    STRIPE_PRICE_IDS,
)
from backend.billing.database import get_db
from backend.billing.dependencies import get_current_user
from backend.billing.models import User

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/billing", tags=["billing"])


class CheckoutSessionRequest(BaseModel):
    tier: str
    success_url: str
    cancel_url: str


class CheckoutSessionResponse(BaseModel):
    session_id: str
    url: str
    publishable_key: str


class PortalRequest(BaseModel):
    return_url: str


class PortalResponse(BaseModel):
    url: str


def _require_stripe():
    if not STRIPE_SECRET_KEY:
        raise HTTPException(
            status_code=status.HTTP_501_NOT_IMPLEMENTED,
            detail="Stripe not configured (missing STRIPE_SECRET_KEY)",
        )
    try:
        import stripe  # type: ignore
    except Exception:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Stripe SDK not available on server",
        )
    stripe.api_key = STRIPE_SECRET_KEY
    return stripe


@router.post("/checkout/session", response_model=CheckoutSessionResponse)
async def create_checkout_session(
    payload: CheckoutSessionRequest,
    db: Session = Depends(get_db),
    user: User = Depends(get_current_user),
):
    try:
        stripe = _require_stripe()

        price_id = STRIPE_PRICE_IDS.get(payload.tier)
        if not price_id:
            logger.error(f"Invalid tier requested: {payload.tier}, available: {list(STRIPE_PRICE_IDS.keys())}")
            raise HTTPException(status_code=400, detail=f"Invalid tier: {payload.tier}")

        # Ensure Stripe customer exists
        if not user.stripe_customer_id:
            try:
                customer = stripe.Customer.create(
                    email=user.email,
                    metadata={"user_id": str(user.id)},
                )
                user.stripe_customer_id = customer.id
                db.commit()
                db.refresh(user)
                logger.info(f"Created Stripe customer {customer.id} for user {user.id}")
            except Exception as e:
                logger.error(f"Failed to create Stripe customer: {e}")
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to create Stripe customer: {str(e)}"
                )

        # Create Checkout Session in subscription mode
        try:
            session = stripe.checkout.Session.create(
                mode="subscription",
                customer=user.stripe_customer_id,
                line_items=[{"price": price_id, "quantity": 1}],
                success_url=payload.success_url,
                cancel_url=payload.cancel_url,
                client_reference_id=str(user.id),
                metadata={"user_id": str(user.id), "tier": payload.tier},
                subscription_data={
                    "metadata": {"user_id": str(user.id), "tier": payload.tier},
                },
                allow_promotion_codes=True,
            )
            logger.info(f"Created checkout session {session.id} for user {user.id}, tier {payload.tier}")
        except Exception as e:
            logger.error(f"Failed to create Stripe checkout session: {e}")
            raise HTTPException(
                status_code=500,
                detail=f"Failed to create checkout session: {str(e)}"
            )

        return CheckoutSessionResponse(
            session_id=session.id,
            url=session.url,
            publishable_key=STRIPE_PUBLISHABLE_KEY or "",
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Unexpected error in create_checkout_session: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error: {str(e)}"
        )


@router.post("/portal", response_model=PortalResponse)
async def create_customer_portal(
    payload: PortalRequest,
    db: Session = Depends(get_db),
    user: User = Depends(get_current_user),
):
    stripe = _require_stripe()

    if not user.stripe_customer_id:
        raise HTTPException(status_code=400, detail="User has no Stripe customer yet")

    portal_session = stripe.billing_portal.Session.create(
        customer=user.stripe_customer_id,
        return_url=payload.return_url,
    )

    return PortalResponse(url=portal_session.url)


