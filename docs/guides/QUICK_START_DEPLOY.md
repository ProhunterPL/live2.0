---
date: 2025-12-23
label: [guide, quick-start]
---

# Quick Start - Split Deploy Wdrożenie

**Szybki start dla wdrożenia architektury Split Deploy**

---

## 🎯 Co Musisz Mieć

- [ ] DigitalOcean account
- [ ] AWS account
- [ ] Supabase project
- [ ] Redis (Redis Labs lub własny)
- [ ] Stripe account (produkty utworzone)
- [ ] Domain name (opcjonalnie)

---

## ⚡ Szybki Start (5 Minut)

### 1. DigitalOcean - Utwórz Droplet

```bash
# Przez CLI (lub dashboard)
doctl compute droplet create live2-api \
  --image ubuntu-22-04-x64 \
  --size s-2vcpu-4gb \
  --region nyc3 \
  --ssh-keys YOUR_SSH_KEY_ID
```

**Zapisuj IP Droplet!**

### 2. Setup Droplet (Automatyczny)

```bash
# SSH do Droplet
ssh root@YOUR_DROPLET_IP

# Pobierz i uruchom setup script
cd /opt
git clone https://github.com/YOUR_REPO/live2.0.git
cd live2.0
bash scripts/deploy/setup_do_droplet.sh
```

### 3. Konfiguracja .env

```bash
# Na Droplet
nano /opt/live2.0/.env

# Wypełnij wszystkie wartości (patrz: docs/guides/SPLIT_DEPLOY_STEP_BY_STEP.md)
```

### 4. Supabase - Utwórz Tabele

1. Otwórz Supabase Dashboard → SQL Editor
2. Wykonaj SQL z `docs/technical/SPLIT_DEPLOY_ARCHITECTURE.md` (sekcja "Model Danych")

### 5. AWS - Setup Infrastructure

```bash
# Lokalnie (z AWS CLI)
bash scripts/deploy/aws_batch_setup.sh

# Utwórz access key dla DO
aws iam create-access-key --user-name live2-do-orchestrator
# Zapisuj Access Key ID i Secret!
```

### 6. Build & Push Docker Image

```bash
# Lokalnie
docker build -t live2-simulation:latest -f docker/simulation.Dockerfile .

# Login do ECR
aws ecr get-login-password --region us-east-1 | \
  docker login --username AWS --password-stdin \
  YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com

# Push
docker tag live2-simulation:latest YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/live2-simulation:latest
docker push YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/live2-simulation:latest
```

### 7. AWS Batch - Utwórz Compute Environment + Queue

**Przez AWS Console:**
1. Batch → Compute environments → Create
2. Name: `live2-compute-env`
3. Type: Managed, Spot
4. Instance types: `c5.2xlarge`, `c5.4xlarge`
5. Min: 0, Max: 32 vCPUs

**Job Queue:**
1. Batch → Job queues → Create
2. Name: `live2-job-queue`
3. Compute environment: `live2-compute-env`

**Job Definition:**
```bash
# Zaktualizuj scripts/deploy/aws_batch_job_definition.json z Twoimi wartościami
aws batch register-job-definition \
  --cli-input-json file://scripts/deploy/aws_batch_job_definition.json
```

### 8. Start Backend na DO

```bash
# Na Droplet
cd /opt/live2.0
source venv/bin/activate

# Migracje
alembic -c backend/billing/migrations/alembic.ini upgrade head

# Start service
sudo systemctl enable live2-backend
sudo systemctl start live2-backend

# Sprawdź status
sudo systemctl status live2-backend
```

### 9. Nginx + SSL

```bash
# Na Droplet
sudo certbot --nginx -d your-domain.com
```

### 10. Test

```bash
# Lokalnie
export API_BASE_URL=https://your-domain.com
bash scripts/deploy/test_deployment.sh
```

---

## 📚 Pełna Dokumentacja

- **Szczegółowy przewodnik:** [`SPLIT_DEPLOY_STEP_BY_STEP.md`](SPLIT_DEPLOY_STEP_BY_STEP.md)
- **Architektura:** [`../technical/SPLIT_DEPLOY_ARCHITECTURE.md`](../technical/SPLIT_DEPLOY_ARCHITECTURE.md)

---

## 🆘 Problemy?

### Backend nie startuje
```bash
sudo journalctl -u live2-backend -f
```

### AWS Batch job nie działa
```bash
aws batch describe-jobs --jobs JOB_ID
aws logs tail /aws/batch/job --follow
```

---

**Czas wdrożenia:** ~2-3 godziny (pierwszy raz)

