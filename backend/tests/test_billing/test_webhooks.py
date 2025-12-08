"""
Tests for Stripe webhook handling.
"""

import pytest
from sqlalchemy.orm import Session
from unittest.mock import Mock
from datetime import date, timedelta

from backend.billing.webhooks import WebhookHandler
from backend.billing.models import User, Subscription
from backend.billing.auth import create_user
# Use conftest db fixture


# Use conftest db fixture


@pytest.fixture
def test_user(db: Session):
    """Create test user with Stripe customer ID."""
    import uuid
    unique_email = f"test_webhook_{uuid.uuid4().hex[:8]}@example.com"
    user = create_user(
        db=db,
        email=unique_email,
        password="password",
        tier="research"
    )
    user.stripe_customer_id = "cus_test123"
    db.commit()
    return user


def test_handle_subscription_created(db: Session, test_user):
    """Test subscription.created event."""
    handler = WebhookHandler(db)
    
    # Ensure user has stripe_customer_id (handler looks for user by this)
    # Use unique customer ID to avoid conflicts with other tests
    import uuid
    unique_customer_id = f"cus_test_{uuid.uuid4().hex[:8]}"
    test_user.stripe_customer_id = unique_customer_id
    db.commit()
    db.refresh(test_user)
    
    # Create subscription first
    unique_sub_id = f"sub_test_{uuid.uuid4().hex[:8]}"
    subscription = Subscription(
        user_id=test_user.id,
        tier="research",
        status="trialing",
        current_period_start=date.today(),
        current_period_end=date.today() + timedelta(days=30),
        stripe_subscription_id=unique_sub_id
    )
    db.add(subscription)
    db.commit()
    
    event = {
        "type": "customer.subscription.created",
        "data": {
            "object": {
                "id": unique_sub_id,
                "customer": unique_customer_id
            }
        }
    }
    
    result = handler.handle_event(event)
    assert result is True
    
    # Verify subscription status updated
    db.refresh(subscription)
    db.refresh(test_user)
    # Handler should find user by stripe_customer_id, then subscription by stripe_subscription_id
    assert subscription.status == "active"
    assert test_user.subscription_status == "active"


def test_handle_subscription_updated(db: Session, test_user):
    """Test subscription.updated event."""
    handler = WebhookHandler(db)
    
    subscription = Subscription(
        user_id=test_user.id,
        tier="research",
        status="active",
        current_period_start=date.today(),
        current_period_end=date.today() + timedelta(days=30),
        stripe_subscription_id="sub_test123_updated"
    )
    db.add(subscription)
    db.commit()
    
    event = {
        "type": "customer.subscription.updated",
        "data": {
            "object": {
                "id": "sub_test123_updated",
                "status": "canceled"
            }
        }
    }
    
    result = handler.handle_event(event)
    assert result is True
    
    # Verify subscription and user status updated
    db.refresh(subscription)
    db.refresh(test_user)
    assert subscription.status == "canceled"
    assert test_user.subscription_status == "cancelled"


def test_handle_payment_succeeded(db: Session, test_user):
    """Test payment.succeeded event."""
    handler = WebhookHandler(db)
    
    subscription = Subscription(
        user_id=test_user.id,
        tier="research",
        status="trialing",
        current_period_start=date.today(),
        current_period_end=date.today() + timedelta(days=30),
        stripe_subscription_id="sub_test123_payment"
    )
    db.add(subscription)
    db.commit()
    
    event = {
        "type": "invoice.payment_succeeded",
        "data": {
            "object": {
                "subscription": "sub_test123_payment"
            }
        }
    }
    
    result = handler.handle_event(event)
    assert result is True
    
    # Verify subscription is active
    db.refresh(subscription)
    db.refresh(test_user)
    assert subscription.status == "active"
    assert test_user.subscription_status == "active"


def test_handle_unhandled_event(db: Session):
    """Test unhandled event type."""
    handler = WebhookHandler(db)
    
    event = {
        "type": "unknown.event.type",
        "data": {
            "object": {}
        }
    }
    
    result = handler.handle_event(event)
    assert result is True  # Should not raise exception

