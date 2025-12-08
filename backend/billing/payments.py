"""
Stripe payment processor for billing module.
"""

from typing import Optional, Dict
import logging

logger = logging.getLogger(__name__)

# Try to import stripe
try:
    import stripe
    STRIPE_AVAILABLE = True
except ImportError:
    STRIPE_AVAILABLE = False
    logger.warning("stripe package not available. Install with: pip install stripe")


class PaymentProcessor:
    """
    Stripe payment processor.
    
    Handles customer creation, subscriptions, and payment processing.
    """
    
    def __init__(self, stripe_secret_key: str):
        """
        Initialize payment processor.
        
        Args:
            stripe_secret_key: Stripe secret key
        """
        if not STRIPE_AVAILABLE:
            raise ImportError("stripe package required. Install with: pip install stripe")
        
        stripe.api_key = stripe_secret_key
        self.stripe = stripe
    
    def create_customer(
        self, 
        email: str, 
        metadata: Optional[Dict] = None
    ) -> str:
        """
        Create Stripe customer.
        
        Args:
            email: Customer email
            metadata: Optional metadata dict
        
        Returns:
            Stripe customer ID
        """
        try:
            customer = self.stripe.Customer.create(
                email=email,
                metadata=metadata or {}
            )
            return customer.id
        except Exception as e:
            logger.error(f"Failed to create Stripe customer: {e}")
            raise
    
    def create_subscription(
        self,
        customer_id: str,
        tier: str,
        payment_method_id: Optional[str] = None
    ) -> str:
        """
        Create Stripe subscription.
        
        Args:
            customer_id: Stripe customer ID
            tier: Tier name
            payment_method_id: Payment method ID (optional)
        
        Returns:
            Stripe subscription ID
        """
        from backend.billing.config import STRIPE_PRICE_IDS
        
        price_id = STRIPE_PRICE_IDS.get(tier)
        if not price_id:
            raise ValueError(f"No price ID for tier: {tier}")
        
        try:
            subscription_data = {
                "customer": customer_id,
                "items": [{"price": price_id}],
            }
            
            if payment_method_id:
                subscription_data["default_payment_method"] = payment_method_id
            
            subscription = self.stripe.Subscription.create(**subscription_data)
            return subscription.id
        except Exception as e:
            logger.error(f"Failed to create Stripe subscription: {e}")
            raise
    
    def update_subscription(
        self,
        subscription_id: str,
        new_tier: str
    ) -> str:
        """
        Update subscription tier.
        
        Args:
            subscription_id: Stripe subscription ID
            new_tier: New tier name
        
        Returns:
            Updated subscription ID
        """
        from backend.billing.config import STRIPE_PRICE_IDS
        
        price_id = STRIPE_PRICE_IDS.get(new_tier)
        if not price_id:
            raise ValueError(f"No price ID for tier: {new_tier}")
        
        try:
            subscription = self.stripe.Subscription.retrieve(subscription_id)
            self.stripe.Subscription.modify(
                subscription_id,
                items=[{
                    "id": subscription["items"]["data"][0].id,
                    "price": price_id
                }]
            )
            return subscription_id
        except Exception as e:
            logger.error(f"Failed to update Stripe subscription: {e}")
            raise
    
    def cancel_subscription_at_period_end(
        self,
        subscription_id: str
    ) -> str:
        """
        Cancel subscription at period end.
        
        Args:
            subscription_id: Stripe subscription ID
        
        Returns:
            Subscription ID
        """
        try:
            self.stripe.Subscription.modify(
                subscription_id,
                cancel_at_period_end=True
            )
            return subscription_id
        except Exception as e:
            logger.error(f"Failed to cancel subscription: {e}")
            raise
    
    def cancel_subscription_immediately(
        self,
        subscription_id: str
    ) -> str:
        """
        Cancel subscription immediately.
        
        Args:
            subscription_id: Stripe subscription ID
        
        Returns:
            Subscription ID
        """
        try:
            self.stripe.Subscription.delete(subscription_id)
            return subscription_id
        except Exception as e:
            logger.error(f"Failed to cancel subscription immediately: {e}")
            raise

