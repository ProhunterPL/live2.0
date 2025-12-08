"""
Tests for payment processing (Stripe).
"""

import pytest
from unittest.mock import Mock, patch, MagicMock

from backend.billing.payments import PaymentProcessor


@pytest.fixture
def mock_stripe():
    """Mock Stripe module."""
    with patch('backend.billing.payments.stripe') as mock_stripe:
        yield mock_stripe


@pytest.fixture
def payment_processor(mock_stripe):
    """Create payment processor with mocked Stripe."""
    processor = PaymentProcessor("sk_test_123")
    return processor


def test_create_customer(payment_processor, mock_stripe):
    """Test customer creation."""
    mock_customer = Mock()
    mock_customer.id = "cus_test123"
    mock_stripe.Customer.create = Mock(return_value=mock_customer)
    
    customer_id = payment_processor.create_customer(
        email="test@example.com",
        metadata={"user_id": "123"}
    )
    
    assert customer_id == "cus_test123"
    mock_stripe.Customer.create.assert_called_once()


def test_create_subscription(payment_processor, mock_stripe):
    """Test subscription creation."""
    mock_subscription = Mock()
    mock_subscription.id = "sub_test123"
    mock_stripe.Subscription.create = Mock(return_value=mock_subscription)
    
    # Mock price IDs from config
    with patch('backend.billing.config.STRIPE_PRICE_IDS', {"research": "price_research"}):
        subscription_id = payment_processor.create_subscription(
            customer_id="cus_test123",
            tier="research",
            payment_method_id="pm_test123"
        )
    
    assert subscription_id == "sub_test123"
    mock_stripe.Subscription.create.assert_called_once()


def test_update_subscription(payment_processor, mock_stripe):
    """Test subscription update."""
    mock_subscription = Mock()
    mock_subscription.id = "sub_test123"
    mock_item = Mock()
    mock_item.id = "si_test123"
    mock_subscription.__getitem__ = Mock(return_value={"data": [mock_item]})
    mock_stripe.Subscription.retrieve = Mock(return_value=mock_subscription)
    mock_stripe.Subscription.modify = Mock(return_value=mock_subscription)
    
    # Mock price IDs from config
    with patch('backend.billing.config.STRIPE_PRICE_IDS', {"pro": "price_pro"}):
        subscription_id = payment_processor.update_subscription(
            subscription_id="sub_test123",
            new_tier="pro"
        )
    
    assert subscription_id == "sub_test123"
    mock_stripe.Subscription.modify.assert_called_once()


def test_cancel_subscription_at_period_end(payment_processor, mock_stripe):
    """Test subscription cancellation at period end."""
    mock_stripe.Subscription.modify = Mock(return_value={"id": "sub_test123"})
    
    subscription_id = payment_processor.cancel_subscription_at_period_end("sub_test123")
    
    assert subscription_id == "sub_test123"
    mock_stripe.Subscription.modify.assert_called_once_with(
        "sub_test123",
        cancel_at_period_end=True
    )


def test_cancel_subscription_immediately(payment_processor, mock_stripe):
    """Test immediate subscription cancellation."""
    mock_stripe.Subscription.delete = Mock(return_value={"id": "sub_test123"})
    
    subscription_id = payment_processor.cancel_subscription_immediately("sub_test123")
    
    assert subscription_id == "sub_test123"
    mock_stripe.Subscription.delete.assert_called_once_with("sub_test123")

