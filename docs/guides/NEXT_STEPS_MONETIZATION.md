# Następne Kroki — Monetyzacja Live 2.0

**Status:** ✅ Backend gotowy | ✅ Testy przechodzą | ⏳ Frontend integration | ⏳ E2E testing

---

## ✅ Co Jest Gotowe

1. **Backend:**
   - ✅ Authentication (register/login)
   - ✅ Billing routes (subscription, usage, checkout, portal)
   - ✅ Stripe integration (checkout, webhooks)
   - ✅ Database (PostgreSQL)
   - ✅ Redis (usage tracking)
   - ✅ All tests passing

2. **Infrastructure:**
   - ✅ Stripe Products/Prices configured
   - ✅ Webhook secret configured
   - ✅ JWT secret configured

---

## 🎯 Następne Kroki (Priorytet)

### 1. Smoke Test E2E (TERAZ)

Przetestuj pełny flow:
```powershell
python scripts/smoke_test_monetization.py
```

**Co testuje:**
- Register → Login → Checkout Session → API Access

**Po teście:**
- Otwórz checkout URL w przeglądarce
- Użyj Stripe test card: `4242 4242 4242 4242`
- Sprawdź czy webhook aktualizuje subscription status

### 2. Frontend Integration (Wysoki Priorytet)

**Dodaj do `frontend/src/components/APIv1Jobs.tsx`:**

#### A. Przycisk "Upgrade" (dla trial users)
```tsx
{user && user.subscription_status === 'trial' && (
  <button onClick={handleUpgrade}>
    Upgrade to {user.tier === 'hobby' ? 'Research' : 'Pro'}
  </button>
)}
```

#### B. Przycisk "Manage Billing" (dla active users)
```tsx
{user && user.subscription_status === 'active' && (
  <button onClick={handleManageBilling}>
    Manage Billing
  </button>
)}
```

#### C. Wyświetlanie Usage Stats
```tsx
// Fetch usage from /api/v1/billing/usage
// Display: reactions used/quota, api_calls used/quota
```

**Funkcje do zaimplementowania:**
- `handleUpgrade()` - call `/api/v1/billing/checkout/session` → redirect to `session.url`
- `handleManageBilling()` - call `/api/v1/billing/portal` → redirect to `portal.url`
- `loadUsage()` - call `/api/v1/billing/usage` → display stats

### 3. Paywall w API v1 (Średni Priorytet)

**Sprawdź czy wszystkie płatne endpointy wymagają aktywnej subskrypcji:**

- ✅ `/api/v1/datasets/generate` - już ma `verify_api_key` (sprawdza `subscription_status`)
- ✅ `/api/v1/jobs/*` - już ma `verify_api_key`
- ⚠️ Sprawdź czy `subscription_status != "active"` blokuje dostęp

**Test:**
```bash
# Z trial subscription
curl -X POST http://localhost:8001/api/v1/datasets/generate \
  -H "X-API-Key: <trial_user_api_key>" \
  -H "Content-Type: application/json" \
  -d '{"dataset_type": "test"}'
# Powinno zwrócić 401/403 jeśli trial nie ma dostępu
```

### 4. Dokumentacja dla Użytkowników (Niski Priorytet)

**Stwórz:**
- `docs/guides/USER_BILLING_GUIDE.md` - jak upgradeować, zarządzać subskrypcją
- `docs/api/BILLING_API.md` - API reference dla billing endpoints

### 5. Monitoring & Alerts (Średni Priorytet)

**Dodaj logowanie:**
- Webhook events (success/failure)
- Subscription status changes
- Payment failures

**Alerts:**
- Payment failed → email notification
- Subscription expired → user notification

---

## 🧪 Test Scenarios

### Scenario 1: New User Flow
1. Register → get trial subscription
2. Use API (should work with trial)
3. Upgrade to paid tier
4. Complete Stripe checkout
5. Webhook activates subscription
6. Continue using API

### Scenario 2: Existing User Upgrade
1. Login with existing account
2. Check current subscription (trial/active)
3. Upgrade tier
4. Complete checkout
5. Verify tier change in DB

### Scenario 3: Payment Failure
1. Simulate payment failure webhook
2. Verify subscription status → "expired" or "past_due"
3. Verify API access blocked
4. Test retry payment flow

---

## 📋 Checklist Przed Launch

### Backend
- [x] All tests passing
- [x] Stripe configured
- [x] Database migrations run
- [ ] E2E smoke test passed
- [ ] Paywall verified (trial vs active)

### Frontend
- [ ] Upgrade button implemented
- [ ] Manage billing button implemented
- [ ] Usage stats displayed
- [ ] Error handling for payment failures

### Documentation
- [ ] User guide (how to upgrade)
- [ ] API reference (billing endpoints)
- [ ] Troubleshooting guide

### Production
- [ ] Stripe webhook endpoint configured (production URL)
- [ ] Monitoring setup (webhook events, payment failures)
- [ ] Email notifications configured
- [ ] Backup/restore procedure for billing DB

---

## 🚀 Quick Start (Dla Implementer)

### 1. Smoke Test
```powershell
python scripts/smoke_test_monetization.py
```

### 2. Frontend Integration
Zobacz: `frontend/src/components/APIv1Jobs.tsx` - dodaj przyciski Upgrade/Manage

### 3. Test Paywall
```powershell
# Test z trial user
python -c "import requests; r = requests.post('http://localhost:8001/api/v1/datasets/generate', headers={'X-API-Key': 'trial_key'}, json={'dataset_type': 'test'}); print(r.status_code)"
```

---

**Ostatnia aktualizacja:** 2025-12-23

