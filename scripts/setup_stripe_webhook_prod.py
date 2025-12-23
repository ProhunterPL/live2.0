#!/usr/bin/env python3
"""
Setup Stripe webhook endpoint for production using Stripe API.

This script creates a webhook endpoint in Stripe Dashboard programmatically.
Alternative to manual setup in Stripe Dashboard.
"""

import sys
import os
from dotenv import load_dotenv

load_dotenv()

def setup_webhook():
    """Create webhook endpoint in Stripe."""
    print("=" * 60)
    print("Stripe Webhook Setup (Production)")
    print("=" * 60)
    
    try:
        import stripe
    except ImportError:
        print("\nERROR: stripe module not installed")
        print("Install with: pip install stripe")
        return False
    
    # Get Stripe secret key
    stripe_secret_key = os.getenv("STRIPE_SECRET_KEY")
    if not stripe_secret_key:
        print("\nERROR: STRIPE_SECRET_KEY not set in .env")
        return False
    
    stripe.api_key = stripe_secret_key
    
    # Get webhook URL from user or env
    webhook_url = os.getenv("STRIPE_WEBHOOK_URL")
    if not webhook_url:
        print("\nEnter your webhook URL (e.g., https://your-domain.com/api/v1/billing/webhooks/stripe):")
        webhook_url = input("> ").strip()
        if not webhook_url:
            print("ERROR: Webhook URL is required")
            return False
    
    print(f"\n[1/2] Creating webhook endpoint...")
    print(f"URL: {webhook_url}")
    
    try:
        # Create webhook endpoint
        endpoint = stripe.WebhookEndpoint.create(
            url=webhook_url,
            enabled_events=[
                "customer.subscription.created",
                "customer.subscription.updated",
                "customer.subscription.deleted",
                "invoice.payment_succeeded",
                "invoice.payment_failed",
            ],
        )
        
        print(f"[2/2] Webhook endpoint created!")
        print(f"\nWebhook ID: {endpoint.id}")
        print(f"Webhook Secret: {endpoint.secret}")
        print(f"\nAdd to your .env file:")
        print(f"STRIPE_WEBHOOK_SECRET={endpoint.secret}")
        
        return True
        
    except stripe.error.StripeError as e:
        print(f"\nERROR: Stripe API error: {e}")
        return False
    except Exception as e:
        print(f"\nERROR: Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = setup_webhook()
    sys.exit(0 if success else 1)

