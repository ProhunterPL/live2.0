# Setup Stripe webhook endpoint using Stripe CLI (PowerShell version)
# This script configures a local webhook endpoint for testing

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "Stripe Webhook Setup (Local Development)" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

# Check if Stripe CLI is installed
$stripeCmd = Get-Command stripe -ErrorAction SilentlyContinue
if (-not $stripeCmd) {
    Write-Host "ERROR: Stripe CLI not found. Install from: https://stripe.com/docs/stripe-cli" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "Step 1: Login to Stripe (if not already logged in)" -ForegroundColor Yellow
Write-Host "Run: stripe login" -ForegroundColor Yellow
Write-Host ""

# Webhook endpoint URL (local development)
$LOCAL_WEBHOOK_URL = "http://localhost:8001/api/v1/billing/webhooks/stripe"

Write-Host "Step 2: Creating webhook endpoint..." -ForegroundColor Yellow
Write-Host "URL: $LOCAL_WEBHOOK_URL" -ForegroundColor Yellow
Write-Host ""
Write-Host "Starting Stripe webhook listener..." -ForegroundColor Green
Write-Host "Press Ctrl+C to stop" -ForegroundColor Yellow
Write-Host ""

# Create webhook endpoint
stripe listen --forward-to "$LOCAL_WEBHOOK_URL" `
    --events customer.subscription.created,customer.subscription.updated,customer.subscription.deleted,invoice.payment_succeeded,invoice.payment_failed

Write-Host ""
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "Webhook Secret will be displayed above (whsec_...)" -ForegroundColor Yellow
Write-Host "Copy it to your .env file as STRIPE_WEBHOOK_SECRET" -ForegroundColor Yellow
Write-Host "============================================================" -ForegroundColor Cyan

