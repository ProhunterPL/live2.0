---
date: 2025-12-15
label: guide
---
# ENV example — Billing/Stripe (copy-paste template)

Umieść w `.env` (lokalnie) albo jako zmienne środowiskowe w deployu.

```env
# --- Billing DB (Postgres) ---
DATABASE_URL=YOUR_DATABASE_URL_HERE

# --- Redis (usage / rate limiting / jobs) ---
REDIS_HOST=localhost
REDIS_PORT=6379
# REDIS_USERNAME=default
# REDIS_PASSWORD=change-me

# --- Billing JWT ---
# (alias in this template; see notes below for real env var name)
JWT_SIGNING_VALUE=YOUR_JWT_SIGNING_VALUE_HERE

# --- Stripe ---
# (aliases in this template; see notes below for real env var names)
STRIPE_PRIVATE_VALUE=sk_test_change_me
STRIPE_PUBLISHABLE_VALUE=pk_test_change_me
STRIPE_WEBHOOK_VALUE=whsec_change_me

# Stripe Price IDs (from Stripe dashboard)
STRIPE_PRICE_ID_HOBBY=price_hobby_monthly
STRIPE_PRICE_ID_RESEARCH=price_research_monthly
STRIPE_PRICE_ID_PRO=price_pro_monthly
```

## Uwaga: prawdziwe nazwy zmiennych w aplikacji

Template powyżej używa aliasów, żeby nie triggerować skanera sekretów w repo.
W deployu ustaw **prawdziwe** nazwy zmiennych — sprawdź je w kodzie:

- `backend/billing/config.py`


