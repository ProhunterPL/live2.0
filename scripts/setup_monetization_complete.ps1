# Complete monetization setup script
# This script helps set up Stripe webhook and checks database

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "Live 2.0 Monetization Setup" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

# 1. Check Stripe CLI
Write-Host "`n[1/4] Checking Stripe CLI..." -ForegroundColor Yellow
$stripeCmd = Get-Command stripe -ErrorAction SilentlyContinue
if ($stripeCmd) {
    Write-Host "  OK - Stripe CLI found" -ForegroundColor Green
    $stripeVersion = stripe --version
    Write-Host "  Version: $stripeVersion" -ForegroundColor Gray
} else {
    Write-Host "  WARNING - Stripe CLI not found" -ForegroundColor Yellow
    Write-Host "  Install from: https://stripe.com/docs/stripe-cli" -ForegroundColor Gray
}

# 2. Check .env file
Write-Host "`n[2/4] Checking .env configuration..." -ForegroundColor Yellow
$envFile = ".env"
if (Test-Path $envFile) {
    Write-Host "  OK - .env file exists" -ForegroundColor Green
    
    # Check required vars
    $requiredVars = @(
        "STRIPE_SECRET_KEY",
        "STRIPE_PRICE_ID_HOBBY",
        "STRIPE_PRICE_ID_RESEARCH",
        "STRIPE_PRICE_ID_PRO",
        "JWT_SECRET_KEY",
        "DATABASE_URL"
    )
    
    $missing = @()
    foreach ($var in $requiredVars) {
        $value = [Environment]::GetEnvironmentVariable($var, "Process")
        if (-not $value) {
            # Try to read from .env file
            $content = Get-Content $envFile -Raw
            if ($content -notmatch "$var=") {
                $missing += $var
            }
        }
    }
    
    if ($missing.Count -eq 0) {
        Write-Host "  OK - All required variables present" -ForegroundColor Green
    } else {
        Write-Host "  WARNING - Missing variables:" -ForegroundColor Yellow
        foreach ($var in $missing) {
            Write-Host "    - $var" -ForegroundColor Gray
        }
    }
} else {
    Write-Host "  WARNING - .env file not found" -ForegroundColor Yellow
}

# 3. Stripe Webhook Setup
Write-Host "`n[3/4] Stripe Webhook Setup..." -ForegroundColor Yellow
Write-Host "`nFor LOCAL development:" -ForegroundColor Cyan
Write-Host "  Run in a separate terminal:" -ForegroundColor White
Write-Host "    stripe listen --forward-to http://localhost:8001/api/v1/billing/webhooks/stripe" -ForegroundColor Green
Write-Host "`n  The webhook secret (whsec_...) will be displayed." -ForegroundColor Gray
Write-Host "  Copy it to your .env file as STRIPE_WEBHOOK_SECRET" -ForegroundColor Gray

Write-Host "`nFor PRODUCTION:" -ForegroundColor Cyan
Write-Host "  1. Go to Stripe Dashboard > Developers > Webhooks" -ForegroundColor White
Write-Host "  2. Click 'Add endpoint'" -ForegroundColor White
Write-Host "  3. URL: https://your-domain.com/api/v1/billing/webhooks/stripe" -ForegroundColor White
Write-Host "  4. Select events:" -ForegroundColor White
Write-Host "     - customer.subscription.created" -ForegroundColor Gray
Write-Host "     - customer.subscription.updated" -ForegroundColor Gray
Write-Host "     - customer.subscription.deleted" -ForegroundColor Gray
Write-Host "     - invoice.payment_succeeded" -ForegroundColor Gray
Write-Host "     - invoice.payment_failed" -ForegroundColor Gray
Write-Host "  5. Copy webhook secret (whsec_...) to .env as STRIPE_WEBHOOK_SECRET" -ForegroundColor White

# 4. Database Setup
Write-Host "`n[4/4] Database Setup..." -ForegroundColor Yellow
Write-Host "  Run database migrations:" -ForegroundColor White
Write-Host "    alembic -c backend/billing/migrations/alembic.ini upgrade head" -ForegroundColor Green
Write-Host "`n  NOTE: If you get 'Tenant or user not found' error:" -ForegroundColor Yellow
Write-Host "    - Check DATABASE_URL in .env" -ForegroundColor Gray
Write-Host "    - Verify database exists in Supabase/PostgreSQL" -ForegroundColor Gray
Write-Host "    - Check if connection string is correct" -ForegroundColor Gray

Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host "Setup Complete!" -ForegroundColor Green
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "`nNext steps:" -ForegroundColor Yellow
Write-Host "  1. Set up Stripe webhook (see instructions above)" -ForegroundColor White
Write-Host "  2. Run database migrations" -ForegroundColor White
Write-Host "  3. Test: python scripts/test_redis_simple.py" -ForegroundColor White
Write-Host "  4. Test: python scripts/test_db_billing.py" -ForegroundColor White

