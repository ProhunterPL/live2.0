---
date: 2025-12-23
label: [guide, deployment]
---

# Konfiguracja DigitalOcean Droplet

**Przewodnik konfiguracji istniejącego Droplet na DO**

---

## ✅ Co Masz Już Gotowe

- ✅ Droplet na DigitalOcean (działający)
- ✅ Domena: live2.world
- ✅ AWS Access Keys
- ✅ Redis
- ✅ Supabase (Postgres + Storage)
- ✅ S3 Storage

---

## 🔧 KROK 1: Połącz się z Droplet

```bash
# SSH do Droplet
ssh root@YOUR_DROPLET_IP
# lub jeśli masz domenę skonfigurowaną:
ssh root@live2.world
```

---

## 🔧 KROK 2: Setup Backend na Droplet

### 2.1 Sklonuj Repozytorium

```bash
# Na Droplet
cd /opt
git clone https://github.com/YOUR_REPO/live2.0.git
cd live2.0
```

### 2.2 Uruchom Setup Script

```bash
# Uruchom automatyczny setup
bash scripts/deploy/setup_do_droplet.sh
```

Lub ręcznie:

```bash
# Install dependencies
apt update && apt upgrade -y
apt install -y python3.11 python3.11-venv python3-pip nginx certbot python3-certbot-nginx

# Setup Python venv
python3.11 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

---

## 🔧 KROK 3: Konfiguracja .env

```bash
# Na Droplet
cd /opt/live2.0
nano .env
```

**Wypełnij wszystkie wartości:**

```bash
# ============================================
# DATABASE (Supabase)
# ============================================
DATABASE_URL=postgresql://postgres.YOUR_PROJECT:YOUR_PASSWORD@aws-0-us-east-1.pooler.supabase.com:6543/postgres

# ============================================
# REDIS
# ============================================
REDIS_HOST=your-redis-host
REDIS_PORT=6379
REDIS_USERNAME=default
REDIS_PASSWORD=YOUR_REDIS_PASSWORD

# ============================================
# JWT
# ============================================
JWT_SECRET_KEY=YOUR_SECURE_RANDOM_32_CHAR_KEY_HERE

# ============================================
# STRIPE
# ============================================
STRIPE_SECRET_KEY=sk_live_YOUR_KEY
STRIPE_PUBLISHABLE_KEY=pk_live_YOUR_KEY
STRIPE_WEBHOOK_SECRET=whsec_YOUR_SECRET
STRIPE_PRICE_ID_HOBBY=price_YOUR_HOBBY_ID
STRIPE_PRICE_ID_RESEARCH=price_YOUR_RESEARCH_ID
STRIPE_PRICE_ID_PRO=price_YOUR_PRO_ID

# ============================================
# AWS (dla Batch - użyj credentials z live2-do-orchestrator)
# ============================================
AWS_ACCESS_KEY_ID=YOUR_AWS_ACCESS_KEY
AWS_SECRET_ACCESS_KEY=YOUR_AWS_SECRET_KEY
AWS_REGION=us-east-1
AWS_BATCH_JOB_QUEUE=live2-job-queue
AWS_BATCH_JOB_DEFINITION=live2-simulation

# ============================================
# SUPABASE (dla job updates)
# ============================================
SUPABASE_URL=https://YOUR_PROJECT.supabase.co
SUPABASE_SERVICE_KEY=YOUR_SERVICE_ROLE_KEY

# ============================================
# APP CONFIG
# ============================================
ENV=prod
API_BASE_URL=https://live2.world
```

**Zabezpiecz plik:**
```bash
chmod 600 .env
```

---

## 🔧 KROK 4: Database Migrations

```bash
# Na Droplet
cd /opt/live2.0
source venv/bin/activate

# Run migrations
alembic -c backend/billing/migrations/alembic.ini upgrade head
```

---

## 🔧 KROK 5: Utwórz Tabele w Supabase

Otwórz Supabase Dashboard → SQL Editor i wykonaj SQL z:
`docs/technical/SPLIT_DEPLOY_ARCHITECTURE.md` (sekcja "Model Danych")

**Kluczowe tabele:**
- `jobs` - status jobów
- `job_artifacts` - artefakty z S3
- RLS policies - bezpieczeństwo

---

## 🔧 KROK 6: Systemd Service

```bash
# Na Droplet
sudo nano /etc/systemd/system/live2-backend.service
```

**Zawartość:**
```ini
[Unit]
Description=Live 2.0 Backend API
After=network.target

[Service]
Type=simple
User=root
WorkingDirectory=/opt/live2.0
Environment="PATH=/opt/live2.0/venv/bin"
ExecStart=/opt/live2.0/venv/bin/gunicorn \
  -w 4 \
  -k uvicorn.workers.UvicornWorker \
  -b 127.0.0.1:8000 \
  backend.api.server:app
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

**Aktywuj:**
```bash
sudo systemctl daemon-reload
sudo systemctl enable live2-backend
sudo systemctl start live2-backend

# Sprawdź status
sudo systemctl status live2-backend
```

---

## 🔧 KROK 7: Nginx Configuration

```bash
# Na Droplet
sudo nano /etc/nginx/sites-available/live2
```

**Zawartość:**
```nginx
server {
    listen 80;
    server_name live2.world www.live2.world;

    # Frontend (static files)
    location / {
        root /opt/live2.0/frontend/dist;
        try_files $uri $uri/ /index.html;
    }

    # Backend API
    location /api {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Health check
    location /status {
        proxy_pass http://127.0.0.1:8000;
    }
}
```

**Aktywuj:**
```bash
sudo ln -s /etc/nginx/sites-available/live2 /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default
sudo nginx -t
sudo systemctl reload nginx
```

---

## 🔧 KROK 8: SSL (Let's Encrypt)

```bash
# Na Droplet
sudo certbot --nginx -d live2.world -d www.live2.world
```

**Auto-renewal jest już skonfigurowany przez certbot.**

---

## 🔧 KROK 9: DNS Configuration

Upewnij się, że DNS wskazuje na Twój Droplet:

```
A     @           -> YOUR_DROPLET_IP
A     www         -> YOUR_DROPLET_IP
```

---

## ✅ Weryfikacja

```bash
# Health check
curl https://live2.world/status/health

# Powinno zwrócić: {"status":"healthy"}
```

---

## 🛡️ Guardrails: Tylko Płacący Użytkownicy

**Ważne:** Backend musi sprawdzać subskrypcję przed uruchomieniem joba!

Implementacja w `backend/jobs/aws_batch.py`:

```python
async def validate_job_request(user: User, request: JobRequest):
    """Validate job request - ONLY paying users can run jobs"""
    
    # Check subscription status
    if user.subscription_status != 'active':
        raise HTTPException(
            403, 
            "Active subscription required. Please upgrade your plan."
        )
    
    # Check tier limits
    tier_limits = TIER_LIMITS[user.tier]
    
    # ... rest of validation
```

**To zapewnia:**
- ✅ Tylko płacący użytkownicy mogą uruchamiać joby
- ✅ AWS Batch uruchamia się tylko dla płacących
- ✅ Po zakończeniu → automatycznie scale to zero
- ✅ Zero kosztów gdy brak płacących użytkowników

---

## 📋 Checklist

- [ ] Backend deployed na DO
- [ ] .env skonfigurowany (wszystkie wartości)
- [ ] Database migrations uruchomione
- [ ] Supabase tables utworzone
- [ ] Systemd service działa
- [ ] Nginx skonfigurowany
- [ ] SSL działa (Let's Encrypt)
- [ ] DNS wskazuje na Droplet
- [ ] Health check działa
- [ ] Guardrails włączone (tylko płacący użytkownicy)

---

**Ostatnia aktualizacja:** 2025-12-23

