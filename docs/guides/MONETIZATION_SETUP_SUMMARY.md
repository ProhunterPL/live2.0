# Podsumowanie Konfiguracji Monetyzacji

**Data:** 2025-12-23  
**Status:** ✅ Redis działa | ⚠️ DB connection issue | ✅ Stripe CLI skonfigurowany

---

## ✅ Co Jest Gotowe

1. **Stripe Products/Prices** — dodane w Stripe Dashboard
2. **Stripe CLI** — zainstalowane i zalogowane (v1.34.0)
3. **JWT Secret Key** — wygenerowany i dodany do `.env`
4. **Redis** — działa (Redis Labs)
5. **Backend Code** — wszystkie endpointy zaimplementowane

---

## ⚠️ Problemy do Rozwiązania

### 1. Database Connection Error

**Błąd:**
```
FATAL: Tenant or user not found
```

**Przyczyna:** `DATABASE_URL` w `.env` wskazuje na nieistniejącą bazę lub nieprawidłowy connection string.

**Rozwiązanie:**
1. Sprawdź `DATABASE_URL` w `.env`
2. Upewnij się, że baza istnieje w Supabase/PostgreSQL
3. Jeśli używasz Supabase, użyj **Connection Pooling URL** (port 6543) zamiast direct connection (port 5432)
4. Format dla Supabase Pooling:
   ```
   DATABASE_URL=postgresql://postgres.[project-ref]:[password]@aws-0-[region].pooler.supabase.com:6543/postgres
   ```

**Test:**
```bash
python scripts/test_db_billing.py
```

### 2. Stripe Webhook Setup

**Dla lokalnego developmentu:**
```bash
# W osobnym terminalu
stripe listen --forward-to http://localhost:8001/api/v1/billing/webhooks/stripe
```

Skopiuj wyświetlony `whsec_...` do `.env` jako `STRIPE_WEBHOOK_SECRET`.

**Dla produkcji:**
1. Stripe Dashboard → Developers → Webhooks → Add endpoint
2. URL: `https://your-domain.com/api/v1/billing/webhooks/stripe`
3. Events:
   - `customer.subscription.created`
   - `customer.subscription.updated`
   - `customer.subscription.deleted`
   - `invoice.payment_succeeded`
   - `invoice.payment_failed`
4. Skopiuj webhook secret do `.env`

---

## 📋 Checklist

- [x] Stripe Products/Prices utworzone
- [x] `STRIPE_PRICE_ID_*` w `.env`
- [x] `STRIPE_SECRET_KEY` w `.env`
- [x] `JWT_SECRET_KEY` w `.env` (zmieniony z default)
- [x] Stripe CLI zainstalowane i zalogowane
- [ ] `STRIPE_WEBHOOK_SECRET` w `.env` (po skonfigurowaniu webhook)
- [ ] `DATABASE_URL` poprawiony (rozwiązać błąd connection)
- [ ] Migracje uruchomione (`alembic upgrade head`)

---

## 🧪 Testy

### Redis
```bash
python scripts/test_redis_simple.py
```
**Status:** ✅ PASSED

### Database
```bash
python scripts/test_db_billing.py
```
**Status:** ❌ FAILED (connection error)

### Stripe Webhook (po setup)
```bash
# Test webhook endpoint
curl -X POST http://localhost:8001/api/v1/billing/webhooks/stripe \
  -H "stripe-signature: ..." \
  -d '{"type": "test"}'
```

---

## 🔧 Następne Kroki

1. **Napraw DATABASE_URL** — sprawdź connection string w Supabase Dashboard
2. **Skonfiguruj Stripe webhook** — użyj Stripe CLI (local) lub Dashboard (prod)
3. **Uruchom migracje** — po naprawieniu DB connection
4. **Smoke test E2E** — register → checkout → webhook → subscription active

---

## 📞 Pomoc

- **Stripe CLI docs:** https://stripe.com/docs/stripe-cli
- **Supabase Connection Strings:** https://supabase.com/docs/guides/database/connecting-to-postgres
- **Billing Guide:** `docs/guides/BILLING_MONETIZATION_MVP_GUIDE.md`

---

**Ostatnia aktualizacja:** 2025-12-23

