---
date: 2025-12-23
label: [guide, deployment]
---

# Split Deploy - Krok po Kroku Wdrożenie

**Praktyczny przewodnik wdrożenia architektury Split Deploy dla Live 2.0**

---

## 📋 Prerequisites

Przed rozpoczęciem upewnij się, że masz:

- [ ] Konto DigitalOcean (z dostępem do API lub dashboard)
- [ ] Konto AWS (z dostępem do Batch, S3, ECR, IAM)
- [ ] Projekt Supabase (Postgres + Auth)
- [ ] Redis (Redis Labs lub własny)
- [ ] Konto Stripe (produkty i ceny utworzone)
- [ ] Domain name (opcjonalnie, ale zalecane)
- [ ] SSH key pair (dla DO)
- [ ] AWS CLI skonfigurowane (`aws configure`)

---

## 🚀 KROK 1: DigitalOcean - Utwórz Droplet

### 1.1 Utwórz Droplet

**Opcja A: Przez Dashboard**
1. Zaloguj się do DigitalOcean
2. Create → Droplets
3. Wybierz:
   - **Image:** Ubuntu 22.04 LTS
   - **Plan:** Basic ($24/mo: 2 vCPU, 4GB RAM, 80GB SSD)
   - **Region:** Najbliższy użytkownikom (np. NYC3, SFO3)
   - **Authentication:** SSH keys (dodaj swój klucz)
4. Create Droplet

**Opcja B: Przez API/CLI**
```bash
# Zainstaluj doctl (DO CLI)
# Windows: choco install doctl
# Mac: brew install doctl
# Linux: https://docs.digitalocean.com/reference/doctl/how-to/install/

# Login
doctl auth init

# Utwórz Droplet
doctl compute droplet create live2-api \
  --image ubuntu-22-04-x64 \
  --size s-2vcpu-4gb \
  --region nyc3 \
  --ssh-keys YOUR_SSH_KEY_ID \
  --wait
```

### 1.2 Zapisuj IP Droplet

```bash
# Pobierz IP
doctl compute droplet list

# Zapisuj IP (będzie potrzebne później)
export DO_IP="YOUR_DROPLET_IP"
```

### 1.3 Połącz się przez SSH

```bash
ssh root@$DO_IP
# lub
ssh root@YOUR_DROPLET_IP
```

---

## 🔧 KROK 2: DigitalOcean - Setup Systemu

### 2.1 Update Systemu

```bash
# Na Droplet (SSH)
apt update && apt upgrade -y
```

### 2.2 Zainstaluj Dependencies

```bash
# Python 3.11
apt install -y python3.11 python3.11-venv python3-pip

# Nginx (reverse proxy)
apt install -y nginx

# Git
apt install -y git

# Certbot (dla SSL)
apt install -y certbot python3-certbot-nginx

# Inne narzędzia
apt install -y curl wget build-essential
```

### 2.3 Utwórz Użytkownika (opcjonalnie, ale zalecane)

```bash
# Zamiast root, użyj dedykowanego użytkownika
adduser live2
usermod -aG sudo live2

# Przełącz się na użytkownika
su - live2
```

### 2.4 Sklonuj Repozytorium

```bash
# Na Droplet
cd /opt
git clone https://github.com/YOUR_REPO/live2.0.git
cd live2.0

# Lub jeśli repo jest prywatne:
# git clone https://YOUR_TOKEN@github.com/YOUR_REPO/live2.0.git
```

---

## 🔐 KROK 3: DigitalOcean - Konfiguracja Environment Variables

### 3.1 Utwórz Plik .env

```bash
# Na Droplet
cd /opt/live2.0
cp .env.example .env  # Jeśli istnieje
# Lub utwórz nowy:
nano .env
```

### 3.2 Wypełnij .env (Production Values)

```bash
# .env na Droplet
# ============================================
# DATABASE (Supabase)
# ============================================
DATABASE_URL=postgresql://postgres.YOUR_PROJECT:YOUR_PASSWORD@aws-0-us-east-1.pooler.supabase.com:6543/postgres

# ============================================
# REDIS
# ============================================
REDIS_HOST=redis-13645.c73.us-east-1-2.ec2.cloud.redislabs.com
REDIS_PORT=13645
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
# AWS (dla Batch)
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
SUPABASE_SERVICE_KEY=YOUR_SERVICE_ROLE_KEY  # Service key (bypass RLS)

# ============================================
# APP CONFIG
# ============================================
ENV=prod
API_BASE_URL=https://api.live2.com  # Twój domain
```

### 3.3 Zabezpiecz .env

```bash
# Tylko root/owner może czytać
chmod 600 .env
```

---

## 🗄️ KROK 4: Supabase - Setup Database Schema

### 4.1 Połącz się z Supabase

1. Zaloguj się do Supabase Dashboard
2. Project Settings → Database → Connection string
3. Skopiuj connection string (dla migrations)

### 4.2 Utwórz Tabele (SQL Editor)

Otwórz SQL Editor w Supabase i wykonaj:

```sql
-- Tabela: jobs
CREATE TABLE IF NOT EXISTS jobs (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  user_id UUID NOT NULL REFERENCES auth.users(id),
  
  -- Status
  status TEXT NOT NULL DEFAULT 'queued',
  -- queued|running|succeeded|failed|canceled
  
  -- Configuration
  params JSONB NOT NULL,
  idempotency_key TEXT UNIQUE,
  
  -- Timing
  created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
  started_at TIMESTAMPTZ,
  finished_at TIMESTAMPTZ,
  
  -- AWS
  aws_batch_job_id TEXT,
  aws_batch_job_arn TEXT,
  
  -- Cost tracking
  cost_estimate DECIMAL(10, 4),
  cost_actual DECIMAL(10, 4),
  
  -- Results
  error TEXT,
  error_details JSONB,
  
  -- Metadata
  progress INTEGER DEFAULT 0,
  metadata JSONB,
  
  CONSTRAINT status_check CHECK (status IN ('queued', 'running', 'succeeded', 'failed', 'canceled'))
);

-- Indeksy
CREATE INDEX IF NOT EXISTS idx_jobs_user_id ON jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);
CREATE INDEX IF NOT EXISTS idx_jobs_created_at ON jobs(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_jobs_aws_batch_job_id ON jobs(aws_batch_job_id);

-- Tabela: job_artifacts
CREATE TABLE IF NOT EXISTS job_artifacts (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  job_id UUID NOT NULL REFERENCES jobs(id) ON DELETE CASCADE,
  
  artifact_type TEXT NOT NULL,
  s3_key TEXT NOT NULL,
  s3_bucket TEXT NOT NULL DEFAULT 'live2-artifacts',
  
  size_bytes BIGINT,
  checksum TEXT,
  
  created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
  
  CONSTRAINT artifact_type_check CHECK (artifact_type IN ('snapshot', 'result', 'log', 'metadata'))
);

-- Indeksy
CREATE INDEX IF NOT EXISTS idx_job_artifacts_job_id ON job_artifacts(job_id);
CREATE INDEX IF NOT EXISTS idx_job_artifacts_type ON job_artifacts(artifact_type);

-- Row Level Security (RLS)
ALTER TABLE jobs ENABLE ROW LEVEL SECURITY;
ALTER TABLE job_artifacts ENABLE ROW LEVEL SECURITY;

-- Policies: Users can view own jobs
CREATE POLICY "Users can view own jobs"
  ON jobs FOR SELECT
  USING (auth.uid() = user_id);

CREATE POLICY "Users can create own jobs"
  ON jobs FOR INSERT
  WITH CHECK (auth.uid() = user_id);

-- Policies: Users can view artifacts of own jobs
CREATE POLICY "Users can view artifacts of own jobs"
  ON job_artifacts FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM jobs
      WHERE jobs.id = job_artifacts.job_id
      AND jobs.user_id = auth.uid()
    )
  );
```

### 4.3 Weryfikacja

```sql
-- Sprawdź czy tabele istnieją
SELECT table_name FROM information_schema.tables 
WHERE table_schema = 'public' 
AND table_name IN ('jobs', 'job_artifacts');

-- Sprawdź RLS
SELECT tablename, rowsecurity FROM pg_tables 
WHERE schemaname = 'public' 
AND tablename IN ('jobs', 'job_artifacts');
```

---

## ☁️ KROK 5: AWS - Setup S3 + ECR + IAM

### 5.1 S3 Bucket

```bash
# Lokalnie (z AWS CLI)
aws s3 mb s3://live2-artifacts --region us-east-1

# Lifecycle policy (delete after 90 days)
cat > lifecycle.json <<EOF
{
  "Rules": [
    {
      "Id": "DeleteOldArtifacts",
      "Status": "Enabled",
      "Expiration": {
        "Days": 90
      }
    }
  ]
}
EOF

aws s3api put-bucket-lifecycle-configuration \
  --bucket live2-artifacts \
  --lifecycle-configuration file://lifecycle.json
```

### 5.2 ECR Repository

```bash
# Utwórz ECR repository
aws ecr create-repository \
  --repository-name live2-simulation \
  --region us-east-1

# Zapisuj URI (będzie potrzebne później)
export ECR_URI=$(aws ecr describe-repositories \
  --repository-names live2-simulation \
  --region us-east-1 \
  --query 'repositories[0].repositoryUri' \
  --output text)

echo "ECR URI: $ECR_URI"
```

### 5.3 IAM Role dla Job Execution

```bash
# Utwórz trust policy
cat > trust-policy.json <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "ecs-tasks.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

# Utwórz role
aws iam create-role \
  --role-name Live2JobRole \
  --assume-role-policy-document file://trust-policy.json

# Attach policies
aws iam attach-role-policy \
  --role-name Live2JobRole \
  --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess

# Lub bardziej restrykcyjna policy (tylko do live2-artifacts):
cat > s3-policy.json <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:PutObject",
        "s3:GetObject"
      ],
      "Resource": "arn:aws:s3:::live2-artifacts/*"
    }
  ]
}
EOF

aws iam put-role-policy \
  --role-name Live2JobRole \
  --policy-name S3ArtifactsAccess \
  --policy-document file://s3-policy.json
```

### 5.4 IAM User dla DO Orchestrator

```bash
# Utwórz user
aws iam create-user --user-name live2-do-orchestrator

# Utwórz policy dla Batch access
cat > batch-policy.json <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "batch:SubmitJob",
        "batch:DescribeJobs",
        "batch:CancelJob"
      ],
      "Resource": "*"
    }
  ]
}
EOF

aws iam put-user-policy \
  --user-name live2-do-orchestrator \
  --policy-name BatchAccess \
  --policy-document file://batch-policy.json

# Utwórz access key
aws iam create-access-key --user-name live2-do-orchestrator

# ZAPISZ Access Key ID i Secret Access Key!
# Będą potrzebne w .env na DO
```

---

## ⚙️ KROK 6: AWS - Setup Batch

### 6.1 Compute Environment

```bash
# Utwórz compute environment (Managed, Spot)
cat > compute-env.json <<EOF
{
  "computeEnvironmentName": "live2-compute-env",
  "type": "MANAGED",
  "state": "ENABLED",
  "computeResources": {
    "type": "SPOT",
    "minvCpus": 0,
    "maxvCpus": 32,
    "desiredvCpus": 0,
    "instanceTypes": ["c5.2xlarge", "c5.4xlarge", "c5.9xlarge"],
    "subnets": [],
    "securityGroupIds": [],
    "bidPercentage": 100,
    "spotIamFleetRole": "arn:aws:iam::YOUR_ACCOUNT:role/aws-ec2-spot-fleet-role"
  },
  "serviceRole": "arn:aws:iam::YOUR_ACCOUNT:role/AWSBatchServiceRole"
}
EOF

# UWAGA: Najpierw utwórz AWSBatchServiceRole i aws-ec2-spot-fleet-role
# (można przez AWS Console → IAM → Roles)

aws batch create-compute-environment --cli-input-json file://compute-env.json
```

**Alternatywa: Przez AWS Console**
1. AWS Console → Batch → Compute environments
2. Create → Managed
3. Name: `live2-compute-env`
4. Instance types: `c5.2xlarge`, `c5.4xlarge`, `c5.9xlarge`
5. Min vCPUs: 0, Max vCPUs: 32
6. Spot: Yes, Bid: 100%

### 6.2 Job Queue

```bash
# Utwórz job queue
aws batch create-job-queue \
  --job-queue-name live2-job-queue \
  --priority 1 \
  --compute-environment-order order=1,computeEnvironment=live2-compute-env \
  --state ENABLED
```

**Alternatywa: Przez AWS Console**
1. Batch → Job queues → Create
2. Name: `live2-job-queue`
3. Priority: 1
4. Compute environment: `live2-compute-env`

### 6.3 Job Definition (Placeholder - zaktualizujemy później)

```bash
# Najpierw potrzebujemy ECR image (Krok 7)
# Na razie utwórz placeholder

aws batch register-job-definition \
  --job-definition-name live2-simulation \
  --type container \
  --container-properties '{
    "image": "YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/live2-simulation:latest",
    "vcpus": 8,
    "memory": 16384,
    "jobRoleArn": "arn:aws:iam::YOUR_ACCOUNT:role/Live2JobRole",
    "environment": [
      {"name": "JOB_ID", "value": "test"},
      {"name": "USER_ID", "value": "test"}
    ]
  }' \
  --retry-strategy attempts=2 \
  --timeout attemptDurationSeconds=3600
```

---

## 🐳 KROK 7: Build & Push Docker Image

### 7.1 Utwórz Dockerfile dla Symulacji

```bash
# Lokalnie (na Twoim komputerze)
cd /path/to/live2.0

# Utwórz docker/simulation.Dockerfile
cat > docker/simulation.Dockerfile <<EOF
FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy backend
COPY backend/ ./backend/

# Set Python path
ENV PYTHONPATH=/app

# Default command (może być nadpisane przez Batch)
CMD ["python", "-m", "backend.sim.run_simulation"]
EOF
```

### 7.2 Build Image

```bash
# Lokalnie
docker build -t live2-simulation:latest -f docker/simulation.Dockerfile .
```

### 7.3 Push do ECR

```bash
# Login do ECR
aws ecr get-login-password --region us-east-1 | \
  docker login --username AWS --password-stdin $ECR_URI

# Tag
docker tag live2-simulation:latest $ECR_URI:latest
docker tag live2-simulation:latest $ECR_URI:v1.0.0

# Push
docker push $ECR_URI:latest
docker push $ECR_URI:v1.0.0
```

### 7.4 Zaktualizuj Job Definition

```bash
# Zaktualizuj Job Definition z prawdziwym image URI
aws batch register-job-definition \
  --job-definition-name live2-simulation \
  --type container \
  --container-properties "{
    \"image\": \"$ECR_URI:latest\",
    \"vcpus\": 8,
    \"memory\": 16384,
    \"jobRoleArn\": \"arn:aws:iam::YOUR_ACCOUNT:role/Live2JobRole\",
    \"environment\": [
      {\"name\": \"SUPABASE_URL\", \"value\": \"YOUR_SUPABASE_URL\"},
      {\"name\": \"SUPABASE_SERVICE_KEY\", \"value\": \"YOUR_SERVICE_KEY\"},
      {\"name\": \"S3_BUCKET\", \"value\": \"live2-artifacts\"},
      {\"name\": \"ENV\", \"value\": \"prod\"}
    ]
  }" \
  --retry-strategy attempts=2 \
  --timeout attemptDurationSeconds=3600
```

---

## 🚀 KROK 8: DigitalOcean - Deploy Backend

### 8.1 Setup Python Environment

```bash
# Na Droplet (SSH)
cd /opt/live2.0

# Utwórz venv
python3.11 -m venv venv
source venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

### 8.2 Run Migrations

```bash
# Na Droplet
cd /opt/live2.0
source venv/bin/activate

# Run migrations
alembic -c backend/billing/migrations/alembic.ini upgrade head
```

### 8.3 Utwórz Systemd Service

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
User=live2
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

**Aktywuj service:**
```bash
sudo systemctl daemon-reload
sudo systemctl enable live2-backend
sudo systemctl start live2-backend

# Sprawdź status
sudo systemctl status live2-backend
```

### 8.4 Test Backend

```bash
# Na Droplet
curl http://localhost:8000/status/health

# Powinno zwrócić: {"status":"healthy"}
```

---

## 🌐 KROK 9: DigitalOcean - Setup Nginx + SSL

### 9.1 Konfiguracja Nginx

```bash
# Na Droplet
sudo nano /etc/nginx/sites-available/live2
```

**Zawartość:**
```nginx
server {
    listen 80;
    server_name api.live2.com;  # Twój domain

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
sudo nginx -t  # Test config
sudo systemctl reload nginx
```

### 9.2 SSL (Let's Encrypt)

```bash
# Na Droplet
sudo certbot --nginx -d api.live2.com

# Auto-renewal (już skonfigurowane przez certbot)
sudo certbot renew --dry-run
```

---

## 🧪 KROK 10: Testy E2E

### 10.1 Health Check

```bash
# Z Twojego komputera
curl https://api.live2.com/status/health
```

### 10.2 Test Auth Flow

```bash
# Register
curl -X POST https://api.live2.com/api/v1/auth/register \
  -H "Content-Type: application/json" \
  -d '{
    "email": "test@example.com",
    "password": "SecurePassword123!",
    "tier": "hobby"
  }'

# Login
curl -X POST https://api.live2.com/api/v1/auth/login \
  -H "Content-Type: application/json" \
  -d '{
    "email": "test@example.com",
    "password": "SecurePassword123!"
  }'

# Zapisuj token z response
```

### 10.3 Test Job Submission

```bash
# Submit job (wymaga aktywnej subskrypcji)
curl -X POST https://api.live2.com/api/v1/jobs \
  -H "Authorization: Bearer YOUR_TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "simulation_config": {
      "mode": "open_chemistry",
      "config": {
        "num_particles": 100,
        "steps": 1000
      }
    }
  }'

# Sprawdź status
curl https://api.live2.com/api/v1/jobs/{job_id} \
  -H "Authorization: Bearer YOUR_TOKEN"
```

---

## ✅ Checklist Wdrożenia

### Pre-Deployment
- [ ] DO Droplet utworzony i dostępny
- [ ] Domain skonfigurowany (DNS → DO IP)
- [ ] Supabase project + tables utworzone
- [ ] Redis dostępny (connection string)
- [ ] Stripe produkty + webhook skonfigurowane
- [ ] AWS account + S3 + ECR + Batch setup

### Deployment
- [ ] Backend API deployed na DO
- [ ] Frontend build deployed (opcjonalnie)
- [ ] Database migrations uruchomione
- [ ] AWS Batch job definition utworzona
- [ ] Docker image w ECR
- [ ] IAM roles/permissions skonfigurowane
- [ ] Nginx + SSL skonfigurowane

### Testing
- [ ] Health checks działają
- [ ] Auth flow (register/login)
- [ ] Stripe checkout flow
- [ ] Job submission (DO → AWS)
- [ ] Job execution (AWS → Supabase update)
- [ ] Artifact download (presigned URLs)

---

## 🆘 Troubleshooting

### Backend nie startuje
```bash
# Sprawdź logi
sudo journalctl -u live2-backend -f

# Sprawdź .env
cat /opt/live2.0/.env

# Test connection do Supabase
python3 -c "from backend.billing.database import get_db; next(get_db())"
```

### AWS Batch job nie startuje
```bash
# Sprawdź job status
aws batch describe-jobs --jobs JOB_ID

# Sprawdź CloudWatch logs
aws logs tail /aws/batch/job --follow
```

### Nginx 502 Bad Gateway
```bash
# Sprawdź czy backend działa
curl http://localhost:8000/status/health

# Sprawdź nginx error log
sudo tail -f /var/log/nginx/error.log
```

---

## 📚 Następne Kroki

Po udanym wdrożeniu:

1. **Monitoring:** Setup CloudWatch alarms, DO monitoring
2. **Backups:** Configure Supabase backups
3. **CI/CD:** Setup GitHub Actions dla auto-deploy
4. **Documentation:** Update API docs, user guides

---

**Ostatnia aktualizacja:** 2025-12-23

