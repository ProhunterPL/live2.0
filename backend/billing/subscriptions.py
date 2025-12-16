"""
Subscription management for billing module.
"""

from typing import Optional
from datetime import date, timedelta
from sqlalchemy.orm import Session
import uuid
import logging

from backend.billing.models import User, Subscription
from backend.billing.payments import PaymentProcessor

logger = logging.getLogger(__name__)


class SubscriptionManager:
    """
    Manages user subscriptions and tier assignments.
    """
    
    TIER_PRICES = {
        "hobby": 29.00,      # USD/month
        "research": 199.00,
        "pro": 999.00,
        "enterprise": None   # Custom pricing
    }
    
    def __init__(self, db: Session, payment_processor: PaymentProcessor):
        """
        Initialize subscription manager.
        
        Args:
            db: SQLAlchemy database session
            payment_processor: PaymentProcessor instance
        """
        self.db = db
        self.payment_processor = payment_processor
    
    def create_subscription(
        self, 
        user_id: uuid.UUID, 
        tier: str,
        payment_method_id: Optional[str] = None
    ) -> Subscription:
        """
        Create subscription and Stripe customer/subscription.
        
        Args:
            user_id: User ID
            tier: Tier name (hobby, research, pro, enterprise)
            payment_method_id: Stripe payment method ID (optional)
        
        Returns:
            Subscription object
        
        Raises:
            ValueError: Invalid tier or user not found
        """
        if tier not in self.TIER_PRICES:
            raise ValueError(f"Invalid tier: {tier}")
        
        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            raise ValueError(f"User not found: {user_id}")
        
        # Create Stripe customer if not exists (only if payment processor available)
        if not user.stripe_customer_id and self.payment_processor:
            try:
                customer_id = self.payment_processor.create_customer(
                    email=user.email,
                    metadata={"user_id": str(user_id)}
                )
                user.stripe_customer_id = customer_id
                self.db.commit()
            except Exception as e:
                logger.error(f"Failed to create Stripe customer: {e}")
                # Continue without Stripe for trial users
                if tier != "hobby":  # Only allow trial for hobby tier
                    raise
        
        # Create Stripe subscription (if payment method provided and processor available)
        subscription_id = None
        if user.stripe_customer_id and payment_method_id and self.payment_processor:
            try:
                subscription_id = self.payment_processor.create_subscription(
                    customer_id=user.stripe_customer_id,
                    tier=tier,
                    payment_method_id=payment_method_id
                )
            except Exception as e:
                logger.error(f"Failed to create Stripe subscription: {e}")
                # Continue without Stripe subscription (trial mode)
        
        # Create subscription record
        now = date.today()
        subscription = Subscription(
            user_id=user_id,
            tier=tier,
            status="active" if subscription_id else "trial",
            current_period_start=now,
            current_period_end=now + timedelta(days=30),
            cancel_at_period_end=False,
            stripe_subscription_id=subscription_id
        )
        
        self.db.add(subscription)
        user.tier = tier
        user.subscription_status = "active" if subscription_id else "trial"
        if subscription_id:
            user.stripe_subscription_id = subscription_id
        self.db.commit()
        self.db.refresh(subscription)
        
        logger.info(f"Created subscription for user {user_id}: tier={tier}, status={subscription.status}")
        return subscription
    
    def update_tier(
        self, 
        user_id: uuid.UUID, 
        new_tier: str
    ) -> Subscription:
        """
        Upgrade/downgrade tier.
        
        Args:
            user_id: User ID
            new_tier: New tier name
        
        Returns:
            Updated Subscription object
        """
        if new_tier not in self.TIER_PRICES:
            raise ValueError(f"Invalid tier: {new_tier}")
        
        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            raise ValueError(f"User not found: {user_id}")
        
        subscription = self.get_active_subscription(user_id)
        if not subscription:
            # Create new subscription if none exists
            return self.create_subscription(user_id, new_tier)
        
        # Update Stripe subscription (if processor available)
        if subscription.stripe_subscription_id and self.payment_processor:
            try:
                self.payment_processor.update_subscription(
                    subscription_id=subscription.stripe_subscription_id,
                    new_tier=new_tier
                )
            except Exception as e:
                logger.error(f"Failed to update Stripe subscription: {e}")
                # Continue with local update
        
        # Update subscription record
        subscription.tier = new_tier
        user.tier = new_tier
        self.db.commit()
        self.db.refresh(subscription)
        
        return subscription
    
    def cancel_subscription(
        self, 
        user_id: uuid.UUID,
        cancel_at_period_end: bool = True
    ) -> Subscription:
        """
        Cancel subscription (at period end or immediately).
        
        Args:
            user_id: User ID
            cancel_at_period_end: If True, cancel at period end; if False, cancel immediately
        
        Returns:
            Updated Subscription object
        """
        subscription = self.get_active_subscription(user_id)
        if not subscription:
            raise ValueError(f"No active subscription for user {user_id}")
        
        user = self.db.query(User).filter(User.id == user_id).first()
        
        if cancel_at_period_end:
            # Cancel at period end
            subscription.cancel_at_period_end = True
            if subscription.stripe_subscription_id and self.payment_processor:
                try:
                    self.payment_processor.cancel_subscription_at_period_end(
                        subscription_id=subscription.stripe_subscription_id
                    )
                except Exception as e:
                    logger.error(f"Failed to cancel Stripe subscription: {e}")
        else:
            # Cancel immediately
            subscription.status = "cancelled"
            user.subscription_status = "cancelled"
            if subscription.stripe_subscription_id and self.payment_processor:
                try:
                    self.payment_processor.cancel_subscription_immediately(
                        subscription_id=subscription.stripe_subscription_id
                    )
                except Exception as e:
                    logger.error(f"Failed to cancel Stripe subscription immediately: {e}")
        
        self.db.commit()
        self.db.refresh(subscription)
        
        return subscription
    
    def get_active_subscription(
        self, 
        user_id: uuid.UUID
    ) -> Optional[Subscription]:
        """
        Get active subscription for user.
        
        Args:
            user_id: User ID
        
        Returns:
            Subscription object or None
        """
        return self.db.query(Subscription).filter(
            Subscription.user_id == user_id,
            Subscription.status == "active"
        ).first()

    def get_current_subscription(self, user_id: uuid.UUID) -> Optional[Subscription]:
        """
        Get the most recent subscription record for user (active/trial/etc.).

        This is used for UI/status endpoints where we want to show something even
        when there is no active subscription yet.
        """
        return (
            self.db.query(Subscription)
            .filter(Subscription.user_id == user_id)
            .order_by(Subscription.created_at.desc())
            .first()
        )

