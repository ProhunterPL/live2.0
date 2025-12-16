"""
Stripe webhook handlers for billing module.
"""

from typing import Dict
import logging
from sqlalchemy.orm import Session
from datetime import datetime, date

from backend.billing.models import User, Subscription
from backend.billing.config import STRIPE_PRICE_IDS

logger = logging.getLogger(__name__)


class WebhookHandler:
    """
    Handles Stripe webhook events.
    """
    
    def __init__(self, db: Session):
        """
        Initialize webhook handler.
        
        Args:
            db: SQLAlchemy database session
        """
        self.db = db
    
    def handle_event(self, event: Dict) -> bool:
        """
        Handle Stripe webhook event.
        
        Args:
            event: Stripe event dict
        
        Returns:
            True if handled successfully
        """
        event_type = event.get("type")
        data = event.get("data", {}).get("object", {})
        
        try:
            if event_type == "customer.subscription.created":
                self._handle_subscription_created(data)
            elif event_type == "customer.subscription.updated":
                self._handle_subscription_updated(data)
            elif event_type == "customer.subscription.deleted":
                self._handle_subscription_deleted(data)
            elif event_type == "invoice.payment_succeeded":
                self._handle_payment_succeeded(data)
            elif event_type == "invoice.payment_failed":
                self._handle_payment_failed(data)
            else:
                logger.info(f"Unhandled webhook event: {event_type}")
            
            return True
        except Exception as e:
            logger.error(f"Failed to handle webhook event {event_type}: {e}", exc_info=True)
            return False
    
    def _handle_subscription_created(self, subscription_data: Dict):
        """Handle subscription.created event."""
        stripe_subscription_id = subscription_data.get("id")
        customer_id = subscription_data.get("customer")
        status = subscription_data.get("status") or "active"

        # Try to derive tier from price id
        tier = None
        try:
            reverse_price = {v: k for k, v in STRIPE_PRICE_IDS.items() if v}
            items = subscription_data.get("items", {}).get("data", [])
            if items:
                price_id = items[0].get("price", {}).get("id")
                tier = reverse_price.get(price_id)
        except Exception:
            tier = None
        if not tier:
            tier = (subscription_data.get("metadata", {}) or {}).get("tier") or "hobby"

        # Convert Stripe timestamps (unix seconds) to dates if present
        def _to_date(ts):
            try:
                return datetime.utcfromtimestamp(int(ts)).date()
            except Exception:
                return date.today()
        
        # Find user by Stripe customer ID
        user = self.db.query(User).filter(
            User.stripe_customer_id == customer_id
        ).first()
        
        if user:
            # Update subscription record
            subscription = self.db.query(Subscription).filter(
                Subscription.stripe_subscription_id == stripe_subscription_id
            ).first()
            
            if subscription:
                subscription.status = status
            else:
                subscription = Subscription(
                    user_id=user.id,
                    tier=tier,
                    status=status,
                    current_period_start=_to_date(subscription_data.get("current_period_start")),
                    current_period_end=_to_date(subscription_data.get("current_period_end")),
                    cancel_at_period_end=bool(subscription_data.get("cancel_at_period_end") or False),
                    stripe_subscription_id=stripe_subscription_id,
                )
                self.db.add(subscription)

            # Update user fields
            user.tier = tier
            user.stripe_subscription_id = stripe_subscription_id
            user.subscription_status = "active" if status == "active" else user.subscription_status
            self.db.commit()
            logger.info(f"Upserted subscription {stripe_subscription_id} (tier={tier}, status={status})")
    
    def _handle_subscription_updated(self, subscription_data: Dict):
        """Handle subscription.updated event."""
        stripe_subscription_id = subscription_data.get("id")
        status = subscription_data.get("status")

        # Derive tier (best effort)
        tier = None
        try:
            reverse_price = {v: k for k, v in STRIPE_PRICE_IDS.items() if v}
            items = subscription_data.get("items", {}).get("data", [])
            if items:
                price_id = items[0].get("price", {}).get("id")
                tier = reverse_price.get(price_id)
        except Exception:
            tier = None
        if not tier:
            tier = (subscription_data.get("metadata", {}) or {}).get("tier")

        def _to_date(ts):
            try:
                return datetime.utcfromtimestamp(int(ts)).date()
            except Exception:
                return date.today()
        
        subscription = self.db.query(Subscription).filter(
            Subscription.stripe_subscription_id == stripe_subscription_id
        ).first()
        
        if subscription:
            subscription.status = status
            if tier:
                subscription.tier = tier
                subscription.user.tier = tier
            # keep period dates in sync if provided
            if subscription_data.get("current_period_start"):
                subscription.current_period_start = _to_date(subscription_data.get("current_period_start"))
            if subscription_data.get("current_period_end"):
                subscription.current_period_end = _to_date(subscription_data.get("current_period_end"))
            if status == "canceled":
                user = subscription.user
                user.subscription_status = "cancelled"
            elif status == "active":
                user = subscription.user
                user.subscription_status = "active"
            self.db.commit()
            logger.info(f"Updated subscription {stripe_subscription_id} to status: {status}")
        else:
            # If we don't have local record yet (e.g., Checkout created it), treat as created
            self._handle_subscription_created(subscription_data)
    
    def _handle_subscription_deleted(self, subscription_data: Dict):
        """Handle subscription.deleted event."""
        stripe_subscription_id = subscription_data.get("id")
        
        subscription = self.db.query(Subscription).filter(
            Subscription.stripe_subscription_id == stripe_subscription_id
        ).first()
        
        if subscription:
            subscription.status = "cancelled"
            user = subscription.user
            user.subscription_status = "cancelled"
            self.db.commit()
            logger.info(f"Deleted subscription {stripe_subscription_id}")
    
    def _handle_payment_succeeded(self, invoice_data: Dict):
        """Handle payment.succeeded event."""
        subscription_id = invoice_data.get("subscription")
        if subscription_id:
            subscription = self.db.query(Subscription).filter(
                Subscription.stripe_subscription_id == subscription_id
            ).first()
            
            if subscription:
                # Subscription is active
                subscription.status = "active"
                user = subscription.user
                user.subscription_status = "active"
                self.db.commit()
                logger.info(f"Payment succeeded for subscription {subscription_id}")
    
    def _handle_payment_failed(self, invoice_data: Dict):
        """Handle payment.failed event."""
        subscription_id = invoice_data.get("subscription")
        logger.warning(f"Payment failed for subscription: {subscription_id}")
        # Best-effort mark local subscription as past_due if exists
        if subscription_id:
            subscription = self.db.query(Subscription).filter(
                Subscription.stripe_subscription_id == subscription_id
            ).first()
            if subscription:
                subscription.status = "past_due"
                subscription.user.subscription_status = "expired"
                self.db.commit()
        # TODO: Implement retry logic / notification / grace period policy

