---
date: 2025-12-15
label: guide
---
# Billing & Monetization MVP (Stripe + API v1) — Deployment Guide

## Cel

Uruchomić **płatne subskrypcje** dla API v1 (SDaaS) z minimalnym, bezpiecznym flow:

- Rejestracja/logowanie → dostajesz JWT + API key
- Checkout (Stripe) → webhook → aktualizacja `subscription_status` w DB
- API v1 działa tylko dla userów z aktywną subskrypcją (paywall na `X-API-Key`)

## Endpointy (MVP)

Wszystkie poniższe endpointy są pod `/api/v1` (bo API v1 jest mountowane na `/api/v1`).

- Auth:
  - `POST /api/v1/auth/register`
  - `POST /api/v1/auth/login`
- Billing:
  - `GET /api/v1/billing/subscription` (JWT lub API key)
  - `GET /api/v1/billing/usage` (JWT lub API key)
  - `POST /api/v1/billing/checkout/session` (JWT lub API key)
  - `POST /api/v1/billing/portal` (JWT lub API key)
  - `POST /api/v1/billing/webhooks/stripe` (Stripe → webhook)

## Wymagane ENV (prod/staging)

Minimalny zestaw dla monetyzacji:

- `DATABASE_URL` (Postgres dla billing)
- `REDIS_HOST`, `REDIS_PORT` (+ opcjonalnie `REDIS_USERNAME`, `REDIS_PASSWORD`)
- `JWT_SECRET_KEY`
- `STRIPE_SECRET_KEY`
- `STRIPE_PUBLISHABLE_KEY`
- `STRIPE_WEBHOOK_SECRET`
- `STRIPE_PRICE_ID_HOBBY`
- `STRIPE_PRICE_ID_RESEARCH`
- `STRIPE_PRICE_ID_PRO`

## DB — migracje (zalecane)

Repo ma Alembica dla billing:

```bash
alembic -c backend/billing/migrations/alembic.ini upgrade head
```

Alternatywa (lokalnie / dev): `init_db()` z `backend/billing/database.py` (bez wersjonowania).

## Stripe — konfiguracja

1. Utwórz Products/Prices w Stripe dla tierów.
2. Wstaw price IDs do ENV (`STRIPE_PRICE_ID_*`).
3. Skonfiguruj webhook endpoint:
   - URL: `/api/v1/billing/webhooks/stripe`
   - Events:
     - `customer.subscription.created`
     - `customer.subscription.updated`
     - `customer.subscription.deleted`
     - `invoice.payment_succeeded`
     - `invoice.payment_failed`
4. Wstaw `STRIPE_WEBHOOK_SECRET` do ENV.

## Smoke test (ręczny)

1. `POST /api/v1/auth/register` → odbierz `access_token` i `api_key`.
2. `POST /api/v1/billing/checkout/session` (Bearer/JWT) → redirect usera na `url`.
3. Po sukcesie płatności Stripe wyśle webhook → DB powinna mieć `subscription.status=active`, a user `subscription_status=active`.
4. Wywołaj płatny endpoint API v1 z `X-API-Key`:
   - np. `POST /api/v1/datasets/generate`


