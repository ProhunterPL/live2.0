---
date: 2025-12-23
label: [guide]
---
# Deployment Checklist - Live 2.0 z Monetyzacją

## ✅ Status Komponentów

### Backend
- [x] FastAPI server z API v1
- [x] Billing module (Stripe integration)
- [x] Authentication (JWT + API keys)
- [x] Database models (PostgreSQL)
- [x] Redis integration (usage tracking)
- [x] Webhook handlers (Stripe)
- [x] Migration system (Alembic)

### Frontend
- [x] React app z Vite
- [x] API v1 Jobs component
- [x] Billing UI (Upgrade/Manage Billing)
- [x] Usage stats display
- [x] Authentication forms

### Infrastructure
- [x] Docker Compose configuration
- [x] Database migrations
- [x] Environment variable templates

---

## 📋 Pre-Deployment Checklist

### 1. Environment Variables (Production)

**Wymagane zmienne środowiskowe:**

```bash
# Database
# Connection string format: protocol://username:password@host:port/database
# Example format (replace with your actual connection string)
DATABASE_URL=YOUR_DATABASE_CONNECTION_STRING  # security-ignore

# Redis
REDIS_HOST=your-redis-host
REDIS_PORT=6379
REDIS_USERNAME=default  # Optional
REDIS_PASSWORD=YOUR_REDIS_PASSWORD  # Optional

# JWT (minimum 32 characters)
JWT_SECRET_KEY=YOUR_JWT_SECRET_KEY  # security-ignore

# Stripe
STRIPE_SECRET_KEY=sk_live_...
STRIPE_PUBLISHABLE_KEY=pk_live_...
STRIPE_WEBHOOK_SECRET=whsec_...
STRIPE_PRICE_ID_HOBBY=price_...
STRIPE_PRICE_ID_RESEARCH=price_...
STRIPE_PRICE_ID_PRO=price_...
```

**Test konfiguracji:**
```bash
python scripts/test_monetization_full.py
```

### 2. Database Setup

**Migracje:**
```bash
# Z projektu root
alembic -c backend/billing/migrations/alembic.ini upgrade head
```

**Alternatywa (dev):**
```python
from backend.billing.database import init_db
init_db()
```

**Weryfikacja:**
```bash
python scripts/test_db_billing.py
```

### 3. Stripe Configuration

**Produkty i ceny:**
- [x] Utworzone w Stripe Dashboard
- [x] Price IDs dodane do `.env`

**Webhook endpoint (PRODUKCJA):**
1. Stripe Dashboard → Developers → Webhooks
2. Add endpoint: `https://your-domain.com/api/v1/billing/webhooks/stripe`
3. Events:
   - `customer.subscription.created`
   - `customer.subscription.updated`
   - `customer.subscription.deleted`
   - `invoice.payment_succeeded`
   - `invoice.payment_failed`
4. Skopiuj `whsec_...` do `STRIPE_WEBHOOK_SECRET`

**Webhook endpoint (LOCAL DEV):**
```bash
stripe listen --forward-to http://localhost:8001/api/v1/billing/webhooks/stripe
```

### 4. Redis Setup

**Lokalnie (Docker):**
```bash
docker-compose up redis
```

**Produkcja (Redis Labs / AWS ElastiCache):**
- [x] Instance utworzona
- [x] Connection string w `.env`
- [x] Test connection: `python scripts/test_redis_billing.py`

### 5. Frontend Build

**Build production:**
```bash
cd frontend
npm run build
```

**Output:** `frontend/dist/` (gotowe do deployu)

**Environment variables dla frontend:**
- `VITE_API_BASE_URL=https://your-api-domain.com` (opcjonalne, default: localhost:8001)

---

## 🚀 Deployment Options

### Option 0: Split Deploy Architecture (Recommended dla Production)

**Nowa architektura hybrydowa: DigitalOcean (SaaS) + AWS (compute on-demand)**

**Dokumentacja:**
- 📖 **Architektura:** [`docs/technical/SPLIT_DEPLOY_ARCHITECTURE.md`](../technical/SPLIT_DEPLOY_ARCHITECTURE.md)
- 🚀 **Quick Start:** [`docs/guides/QUICK_START_DEPLOY.md`](QUICK_START_DEPLOY.md)
- 📋 **Krok po kroku:** [`docs/guides/SPLIT_DEPLOY_STEP_BY_STEP.md`](SPLIT_DEPLOY_STEP_BY_STEP.md)

**Kluczowe zalety:**
- ✅ Minimalny stały koszt (~$40/mo MVP)
- ✅ Skalowanie bez refaktoru (DO + AWS Batch)
- ✅ Bezpieczne wykonywanie zadań (AWS IAM, VPC)
- ✅ Szybkie uruchomienie monetyzacji (DO always-on)

**Komponenty:**
- **DigitalOcean:** Backend API + Frontend (always-on)
- **AWS:** Batch compute (on-demand, auto-scale to zero)
- **Supabase:** Postgres + Auth (współdzielone)
- **Redis:** Rate limiting + cache (współdzielone)

**Pomocnicze skrypty:**
- `scripts/deploy/setup_do_droplet.sh` - Automatyczny setup DO Droplet
- `scripts/deploy/aws_batch_setup.sh` - Setup AWS infrastructure
- `scripts/deploy/test_deployment.sh` - Testy E2E po wdrożeniu

### Option 1: Docker Compose (Recommended dla dev/staging)

```bash
# 1. Skonfiguruj .env
cp .env.example .env
# Edytuj .env z production values

# 2. Uruchom migracje
docker-compose run backend alembic -c backend/billing/migrations/alembic.ini upgrade head

# 3. Start services
docker-compose up -d

# 4. Verify
curl http://localhost:8000/status/health
curl http://localhost:3000
```

### Option 2: AWS EC2 (Production)

**Prerequisites:**
- EC2 instance (t3.xlarge lub większy)
- Security group z portami: 80, 443, 8000, 22
- Domain name + SSL certificate

**Steps:**
1. SSH do instance
2. Install dependencies:
   ```bash
   sudo apt update && sudo apt upgrade -y
   sudo apt install -y python3.11 python3-pip nginx postgresql redis-server
   ```
3. Clone repo:
   ```bash
   git clone https://github.com/your-repo/live2.0.git
   cd live2.0
   ```
4. Setup environment:
   ```bash
   cp .env.example .env
   # Edytuj .env
   pip3 install -r requirements.txt
   ```
5. Run migrations:
   ```bash
   alembic -c backend/billing/migrations/alembic.ini upgrade head
   ```
6. Setup systemd services (backend + frontend)
7. Configure Nginx reverse proxy
8. Setup SSL (Let's Encrypt)

**Szczegóły:** `docs/technical/aws/AWS_QUICK_START.md`

### Option 3: Cloud Platform (Heroku/Railway/Render)

**Heroku:**
```bash
heroku create live2-app
heroku addons:create heroku-postgresql:hobby-dev
heroku addons:create heroku-redis:hobby-dev
heroku config:set STRIPE_SECRET_KEY=...
# ... inne env vars
git push heroku main
heroku run alembic -c backend/billing/migrations/alembic.ini upgrade head
```

---

## 🧪 Post-Deployment Testing

### 1. Health Checks
```bash
# Backend
curl https://your-domain.com/status/health

# Frontend
curl https://your-domain.com
```

### 2. Smoke Test E2E
```bash
python scripts/smoke_test_monetization.py
# Ustaw API_BASE_URL w .env lub jako argument
```

### 3. Manual Testing
1. **Register user:**
   ```bash
   curl -X POST https://your-domain.com/api/v1/auth/register \
     -H "Content-Type: application/json" \
     -d '{"email":"test@example.com","password":"YOUR_PASSWORD","tier":"hobby"}'
   ```

2. **Login:**
   ```bash
   curl -X POST https://your-domain.com/api/v1/auth/login \
     -H "Content-Type: application/json" \
     -d '{"email":"test@example.com","password":"YOUR_PASSWORD"}'
   ```

3. **Create checkout session:**
   ```bash
   curl -X POST https://your-domain.com/api/v1/billing/checkout/session \
     -H "X-API-Key: sk_live_..." \
     -H "Content-Type: application/json" \
     -d '{"tier":"hobby","success_url":"https://your-domain.com/success","cancel_url":"https://your-domain.com/cancel"}'
   ```

4. **Test webhook (Stripe CLI):**
   ```bash
   stripe trigger customer.subscription.created
   ```

---

## ⚠️ Security Checklist

- [ ] `JWT_SECRET_KEY` zmieniony z default value
- [ ] `DATABASE_URL` używa bezpiecznego hasła
- [ ] HTTPS enabled (SSL certificate)
- [ ] CORS configured dla production domain
- [ ] Rate limiting enabled
- [ ] Error messages nie ujawniają internal details
- [ ] Secrets w environment variables, nie w kodzie
- [ ] Database backups configured
- [ ] Logging configured (bez sensitive data)

---

## 📊 Monitoring

**Recommended:**
- Application logs (backend + frontend)
- Database connection monitoring
- Redis connection monitoring
- Stripe webhook delivery monitoring
- Error tracking (Sentry - już zintegrowane)
- Uptime monitoring (UptimeRobot / Pingdom)

---

## 🔄 Rollback Plan

**Jeśli coś pójdzie nie tak:**

1. **Database rollback:**
   ```bash
   alembic -c backend/billing/migrations/alembic.ini downgrade -1
   ```

2. **Code rollback:**
   ```bash
   git checkout <previous-commit>
   # Redeploy
   ```

3. **Environment rollback:**
   - Przywróć poprzednie `.env`
   - Restart services

---

## 📝 Notes

- **Database:** Użyj connection pooling dla production (Supabase Pooling URL)
- **Redis:** Użyj Redis Labs lub AWS ElastiCache dla production
- **Stripe:** Użyj test mode dla staging, live mode dla production
- **Frontend:** Build output (`dist/`) może być serwowany przez Nginx lub CDN
- **Backend:** Użyj Gunicorn + Uvicorn workers dla production (nie uvicorn bezpośrednio)

---

## ✅ Final Checklist

Przed ogłoszeniem "production ready":

- [ ] Wszystkie environment variables ustawione
- [ ] Database migrations uruchomione
- [ ] Stripe webhook skonfigurowany
- [ ] Redis connection działa
- [ ] Smoke test E2E passed
- [ ] HTTPS enabled
- [ ] Monitoring configured
- [ ] Backups configured
- [ ] Documentation updated
- [ ] Team notified

---

**Ostatnia aktualizacja:** 2025-12-23

