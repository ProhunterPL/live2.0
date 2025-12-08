"""
Tests for subscription management.
"""

import pytest
from sqlalchemy.orm import Session
from unittest.mock import Mock, MagicMock
from datetime import date, timedelta

from backend.billing.subscriptions import SubscriptionManager
from backend.billing.models import User, Subscription
from backend.billing.auth import create_user
# Use conftest db fixture


# Use conftest db fixture


@pytest.fixture
def mock_payment_processor():
    """Create mock payment processor."""
    # Create a real PaymentProcessor instance but with mocked Stripe
    from backend.billing.payments import PaymentProcessor
    from unittest.mock import patch
    
    with patch('backend.billing.payments.stripe') as mock_stripe:
        processor = PaymentProcessor("sk_test_123")
        # Mock Stripe methods
        mock_customer = Mock()
        mock_customer.id = "cus_test123"
        mock_stripe.Customer.create = Mock(return_value=mock_customer)
        
        mock_subscription = Mock()
        mock_subscription.id = "sub_test123"
        mock_stripe.Subscription.create = Mock(return_value=mock_subscription)
        mock_stripe.Subscription.modify = Mock(return_value=mock_subscription)
        mock_stripe.Subscription.delete = Mock(return_value=mock_subscription)
        mock_stripe.Subscription.retrieve = Mock(return_value=mock_subscription)
        
        yield processor


@pytest.fixture
def test_user(db: Session):
    """Create test user."""
    import uuid
    unique_email = f"test_subscription_{uuid.uuid4().hex[:8]}@example.com"
    return create_user(
        db=db,
        email=unique_email,
        password="password",
        tier="hobby"
    )


def test_create_subscription(db: Session, mock_payment_processor, test_user):
    """Test subscription creation."""
    from unittest.mock import patch
    
    # Mock STRIPE_PRICE_IDS
    with patch('backend.billing.config.STRIPE_PRICE_IDS', {"research": "price_research"}):
        manager = SubscriptionManager(db, mock_payment_processor)
        
        subscription = manager.create_subscription(
            user_id=test_user.id,
            tier="research"
        )
        
        assert subscription.id is not None
        assert subscription.tier == "research"
        assert subscription.status in ["active", "trial"]
        assert subscription.user_id == test_user.id
        assert subscription.current_period_start == date.today()
        assert subscription.current_period_end == date.today() + timedelta(days=30)


def test_create_subscription_invalid_tier(db: Session, mock_payment_processor, test_user):
    """Test subscription creation with invalid tier."""
    manager = SubscriptionManager(db, mock_payment_processor)
    
    with pytest.raises(ValueError, match="Invalid tier"):
        manager.create_subscription(
            user_id=test_user.id,
            tier="invalid_tier"
        )


def test_update_tier(db: Session, mock_payment_processor, test_user):
    """Test tier update."""
    from unittest.mock import patch
    
    # Mock STRIPE_PRICE_IDS
    with patch('backend.billing.config.STRIPE_PRICE_IDS', {"hobby": "price_hobby", "pro": "price_pro"}):
        manager = SubscriptionManager(db, mock_payment_processor)
        
        # Create initial subscription
        subscription = manager.create_subscription(
            user_id=test_user.id,
            tier="hobby"
        )
        
        # Upgrade to pro
        updated = manager.update_tier(test_user.id, "pro")
        
        assert updated.tier == "pro"
        db.refresh(test_user)
        assert test_user.tier == "pro"


def test_cancel_subscription_at_period_end(db: Session, mock_payment_processor, test_user):
    """Test subscription cancellation at period end."""
    from unittest.mock import patch
    
    with patch('backend.billing.config.STRIPE_PRICE_IDS', {"research": "price_research"}):
        manager = SubscriptionManager(db, mock_payment_processor)
        
        # Create subscription first
        subscription = manager.create_subscription(
            user_id=test_user.id,
            tier="research"
        )
        # Ensure it's active
        subscription.status = "active"
        db.commit()
        db.refresh(subscription)
        
        cancelled = manager.cancel_subscription(
            user_id=test_user.id,
            cancel_at_period_end=True
        )
        
        assert cancelled.cancel_at_period_end is True


def test_cancel_subscription_immediately(db: Session, mock_payment_processor, test_user):
    """Test immediate subscription cancellation."""
    from unittest.mock import patch
    
    with patch('backend.billing.config.STRIPE_PRICE_IDS', {"research": "price_research"}):
        manager = SubscriptionManager(db, mock_payment_processor)
        
        # Create subscription first
        subscription = manager.create_subscription(
            user_id=test_user.id,
            tier="research"
        )
        # Ensure it's active
        subscription.status = "active"
        db.commit()
        db.refresh(subscription)
        
        cancelled = manager.cancel_subscription(
            user_id=test_user.id,
            cancel_at_period_end=False
        )
        
        assert cancelled.status == "cancelled"
        db.refresh(test_user)
        assert test_user.subscription_status == "cancelled"


def test_get_active_subscription(db: Session, mock_payment_processor, test_user):
    """Test getting active subscription."""
    from unittest.mock import patch
    
    with patch('backend.billing.config.STRIPE_PRICE_IDS', {"pro": "price_pro"}):
        manager = SubscriptionManager(db, mock_payment_processor)
        
        # No subscription yet
        subscription = manager.get_active_subscription(test_user.id)
        assert subscription is None
        
        # Create subscription
        created = manager.create_subscription(
            user_id=test_user.id,
            tier="pro"
        )
        # Ensure it's active
        created.status = "active"
        db.commit()
        db.refresh(created)
        
        # Get active subscription
        subscription = manager.get_active_subscription(test_user.id)
        assert subscription is not None
        assert subscription.id == created.id

