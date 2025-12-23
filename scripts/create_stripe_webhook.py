#!/usr/bin/env python3
"""
Create Stripe webhook endpoint using Stripe API.
"""

import sys
import os
from dotenv import load_dotenv

load_dotenv()

def create_webhook():
    """Create webhook endpoint in Stripe."""
    print("=" * 60)
    print("Creating Stripe Webhook Endpoint")
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
    
    # For local development, use Stripe CLI forwarding
    # For production, use actual URL
    print("\nChoose webhook type:")
    print("1. Local development (use Stripe CLI forwarding)")
    print("2. Production (create endpoint in Stripe Dashboard)")
    
    choice = input("Enter choice (1 or 2): ").strip()
    
    if choice == "1":
        print("\n" + "=" * 60)
        print("Local Development Setup")
        print("=" * 60)
        print("\nRun this command in a separate terminal:")
        print("  stripe listen --forward-to http://localhost:8001/api/v1/billing/webhooks/stripe")
        print("\nThe webhook secret (whsec_...) will be displayed.")
        print("Copy it to your .env file as STRIPE_WEBHOOK_SECRET")
        return True
    
    elif choice == "2":
        webhook_url = os.getenv("STRIPE_WEBHOOK_URL")
        if not webhook_url:
            print("\nEnter your production webhook URL:")
            print("(e.g., https://your-domain.com/api/v1/billing/webhooks/stripe)")
            webhook_url = input("> ").strip()
            if not webhook_url:
                print("ERROR: Webhook URL is required")
                return False
        
        print(f"\n[1/2] Creating webhook endpoint...")
        print(f"URL: {webhook_url}")
        
        try:
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
            print(f"\n" + "=" * 60)
            print("Add to your .env file:")
            print(f"STRIPE_WEBHOOK_SECRET={endpoint.secret}")
            print("=" * 60)
            
            return True
            
        except stripe.error.StripeError as e:
            print(f"\nERROR: Stripe API error: {e}")
            return False
    
    else:
        print("ERROR: Invalid choice")
        return False


if __name__ == "__main__":
    success = create_webhook()
    sys.exit(0 if success else 1)

