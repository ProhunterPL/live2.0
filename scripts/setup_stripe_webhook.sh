#!/bin/bash
# Setup Stripe webhook endpoint using Stripe CLI
# This script configures a local webhook endpoint for testing

echo "============================================================"
echo "Stripe Webhook Setup (Local Development)"
echo "============================================================"

# Check if Stripe CLI is installed
if ! command -v stripe &> /dev/null; then
    echo "ERROR: Stripe CLI not found. Install from: https://stripe.com/docs/stripe-cli"
    exit 1
fi

echo ""
echo "Step 1: Login to Stripe (if not already logged in)"
echo "Run: stripe login"
echo ""

# Webhook endpoint URL (local development)
LOCAL_WEBHOOK_URL="http://localhost:8001/api/v1/billing/webhooks/stripe"

echo "Step 2: Creating webhook endpoint..."
echo "URL: $LOCAL_WEBHOOK_URL"
echo ""

# Create webhook endpoint
stripe listen --forward-to "$LOCAL_WEBHOOK_URL" \
    --events customer.subscription.created,customer.subscription.updated,customer.subscription.deleted,invoice.payment_succeeded,invoice.payment_failed

echo ""
echo "============================================================"
echo "Webhook Secret will be displayed above (whsec_...)"
echo "Copy it to your .env file as STRIPE_WEBHOOK_SECRET"
echo "============================================================"

