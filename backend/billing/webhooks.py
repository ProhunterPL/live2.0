"""
Stripe webhook handlers for billing module.
"""

from typing import Dict
import logging
from sqlalchemy.orm import Session

from backend.billing.models import User, Subscription

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
        
        # Find user by Stripe customer ID
        user = self.db.query(User).filter(
            User.stripe_customer_id == customer_id
        ).first()
        
        if user:
            # Update subscription record
            subscription = self.db.query(Subscription).filter(
                Subscription.user_id == user.id,
                Subscription.stripe_subscription_id == stripe_subscription_id
            ).first()
            
            if subscription:
                subscription.status = "active"
                user.subscription_status = "active"
                self.db.commit()
                logger.info(f"Updated subscription {stripe_subscription_id} to active")
    
    def _handle_subscription_updated(self, subscription_data: Dict):
        """Handle subscription.updated event."""
        stripe_subscription_id = subscription_data.get("id")
        status = subscription_data.get("status")
        
        subscription = self.db.query(Subscription).filter(
            Subscription.stripe_subscription_id == stripe_subscription_id
        ).first()
        
        if subscription:
            subscription.status = status
            if status == "canceled":
                user = subscription.user
                user.subscription_status = "cancelled"
            elif status == "active":
                user = subscription.user
                user.subscription_status = "active"
            self.db.commit()
            logger.info(f"Updated subscription {stripe_subscription_id} to status: {status}")
    
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
        # TODO: Implement retry logic or notification

