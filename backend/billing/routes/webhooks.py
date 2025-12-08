"""
Stripe webhook routes.
"""

from fastapi import APIRouter, Depends, Request, HTTPException, status
from sqlalchemy.orm import Session
import logging

from backend.billing.database import get_db
from backend.billing.webhooks import WebhookHandler
from backend.billing.config import STRIPE_WEBHOOK_SECRET

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/billing/webhooks", tags=["webhooks"])


@router.post("/stripe")
async def stripe_webhook(
    request: Request,
    db: Session = Depends(get_db)
):
    """
    Handle Stripe webhook events.
    
    Verifies webhook signature and processes events.
    """
    # Get webhook payload
    payload = await request.body()
    sig_header = request.headers.get("stripe-signature")
    
    if not sig_header:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Missing stripe-signature header"
        )
    
    # Verify webhook signature
    try:
        import stripe
        event = stripe.Webhook.construct_event(
            payload, sig_header, STRIPE_WEBHOOK_SECRET
        )
    except ValueError as e:
        logger.error(f"Invalid webhook payload: {e}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Invalid payload"
        )
    except Exception as e:
        # Handle both stripe.error.SignatureVerificationError and generic errors
        if "SignatureVerificationError" in str(type(e).__name__):
            logger.error(f"Invalid webhook signature: {e}")
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Invalid signature"
            )
        else:
            logger.error(f"Stripe webhook error: {e}")
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Failed to verify webhook"
            )
    
    # Handle event
    webhook_handler = WebhookHandler(db)
    success = webhook_handler.handle_event(event)
    
    if success:
        return {"status": "ok"}
    else:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to process webhook"
        )

