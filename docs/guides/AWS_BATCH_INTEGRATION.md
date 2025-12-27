---
date: 2025-12-27
label: [guide, aws, batch, integration]
---

# AWS Batch Integration - Kompletny Przewodnik

**Status:** ✅ Zintegrowane z backendem

---

## 📋 Co Zostało Zaimplementowane

### 1. AWS Batch Client (`backend/api/v1/aws_batch.py`)

Moduł do zarządzania jobami AWS Batch:

- ✅ `submit_job()` - Submituje job do AWS Batch
- ✅ `get_job_status()` - Sprawdza status joba w Batch
- ✅ `cancel_job()` - Anuluje job w Batch
- ✅ `list_job_artifacts()` - Listuje artefakty z S3
- ✅ `generate_presigned_url()` - Generuje presigned URLs dla S3

### 2. Zaktualizowany JobProcessor (`backend/api/v1/jobs.py`)

- ✅ Integracja z AWS Batch dla `run_simulation` jobs
- ✅ Automatyczne aktualizowanie statusu z Batch
- ✅ Metoda `cancel_job()` do anulowania jobów
- ✅ Obsługa `aws_batch_job_id` w job data

### 3. Nowe Endpointy API (`backend/api/v1/routes/jobs.py`)

#### POST `/api/v1/jobs`
Uruchamia nowy job (simulation lub dataset generation).

**Request:**
```json
{
  "job_type": "run_simulation",
  "params": {
    "simulation_config": {
      "mode": "open_chemistry",
      "config": {
        "num_particles": 1000,
        "temperature": 300,
        "steps": 50000
      }
    }
  },
  "idempotency_key": "optional-key",
  "webhook_url": "https://example.com/webhook"
}
```

**Response:**
```json
{
  "job_id": "job_abc123",
  "status": "queued",
  "estimated_time": 3600
}
```

#### POST `/api/v1/jobs/{job_id}/cancel`
Anuluje uruchomiony job.

**Response:**
```json
{
  "job_id": "job_abc123",
  "status": "cancelled",
  "progress": 0,
  "error": null,
  "created_at": "2025-12-27T10:00:00Z",
  "completed_at": "2025-12-27T10:05:00Z"
}
```

#### GET `/api/v1/jobs/{job_id}/artifacts`
Listuje artefakty (S3 objects) dla zakończonego joba.

**Response:**
```json
{
  "job_id": "job_abc123",
  "artifacts": [
    {
      "key": "prod/user123/job_abc123/results.h5",
      "size": 1048576,
      "url": "https://s3.amazonaws.com/...?presigned=...",
      "last_modified": "2025-12-27T10:30:00Z"
    }
  ]
}
```

---

## 🔧 Konfiguracja

### Environment Variables (DO Droplet `.env`)

```bash
# AWS Batch
AWS_ACCESS_KEY_ID=YOUR_ACCESS_KEY
AWS_SECRET_ACCESS_KEY=YOUR_SECRET_KEY
AWS_REGION=eu-central-1
AWS_BATCH_JOB_QUEUE=live2-job-queue
AWS_BATCH_JOB_DEFINITION=live2-simulation

# S3
AWS_S3_BUCKET=live2-artifacts

# Supabase (opcjonalnie, dla job updates)
SUPABASE_URL=https://YOUR_PROJECT.supabase.co
SUPABASE_SERVICE_KEY=YOUR_SERVICE_KEY

# Redis
REDIS_HOST=YOUR_REDIS_HOST
REDIS_PORT=6379
REDIS_PASSWORD=YOUR_REDIS_PASSWORD
```

---

## 🧪 Testowanie

### 1. Test Integracji AWS Batch

```bash
# Na DO Droplet (gdzie są zainstalowane zależności)
cd /opt/live2.0
source venv/bin/activate
python scripts/test_aws_batch_integration.py
```

### 2. Test przez API

```bash
# Start job
curl -X POST http://localhost:8000/api/v1/jobs \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "job_type": "run_simulation",
    "params": {
      "simulation_config": {
        "mode": "open_chemistry",
        "config": {
          "num_particles": 100,
          "temperature": 300,
          "steps": 1000
        }
      }
    }
  }'

# Check status
curl -X GET http://localhost:8000/api/v1/jobs/{job_id} \
  -H "Authorization: Bearer YOUR_API_KEY"

# Cancel job
curl -X POST http://localhost:8000/api/v1/jobs/{job_id}/cancel \
  -H "Authorization: Bearer YOUR_API_KEY"

# List artifacts
curl -X GET http://localhost:8000/api/v1/jobs/{job_id}/artifacts \
  -H "Authorization: Bearer YOUR_API_KEY"
```

---

## 🔄 Flow Joba

1. **User submits job** → `POST /api/v1/jobs`
2. **Backend creates job** → Redis: `job:{job_id}`
3. **Backend submits to AWS Batch** → `batch_client.submit_job()`
4. **Backend stores batch_job_id** → Redis: `job:{job_id}` → `aws_batch_job_id`
5. **AWS Batch starts job** → EC2 instance spins up (minvCpus=0 → scale up)
6. **Job runs** → Docker container executes simulation
7. **Job uploads artifacts** → S3: `s3://live2-artifacts/prod/{user_id}/{job_id}/...`
8. **Job completes** → Status updated in Batch
9. **Backend polls status** → `GET /api/v1/jobs/{job_id}` → Updates from Batch
10. **User downloads artifacts** → `GET /api/v1/jobs/{job_id}/artifacts` → Presigned URLs

---

## 📊 Status Mapping

| AWS Batch Status | Our Status | Description |
|-----------------|------------|-------------|
| SUBMITTED | queued | Job submitted, waiting |
| PENDING | queued | Job pending execution |
| RUNNABLE | queued | Job ready to run |
| RUNNING | running | Job executing |
| SUCCEEDED | completed | Job completed successfully |
| FAILED | failed | Job failed |
| CANCELLED | cancelled | Job cancelled |

---

## 🐛 Troubleshooting

### Problem: Job nie startuje w Batch

**Sprawdź:**
1. Compute Environment status: `aws batch describe-compute-environments`
2. Job Queue status: `aws batch describe-job-queues`
3. Job Definition: `aws batch describe-job-definitions --job-definition-name live2-simulation`
4. CloudWatch Logs: Sprawdź logi dla joba

### Problem: Job fails immediately

**Sprawdź:**
1. Docker image w ECR: `aws ecr describe-images --repository-name live2-simulation`
2. IAM Role permissions: `aws iam get-role --role-name Live2JobRole`
3. Environment variables w Job Definition
4. CloudWatch Logs dla szczegółów błędu

### Problem: Artifacts nie są dostępne

**Sprawdź:**
1. S3 bucket permissions: `aws s3 ls s3://live2-artifacts/prod/`
2. IAM Role ma uprawnienia do S3: `aws iam get-role-policy --role-name Live2JobRole --policy-name S3ArtifactsAccess`
3. Presigned URL expiration (domyślnie 1h)

---

## ✅ Checklist Deployment

- [x] AWS Batch infrastructure setup
- [x] Docker image pushed to ECR
- [x] Job Definition created
- [x] Backend integration code
- [x] API endpoints implemented
- [ ] Environment variables configured on DO Droplet
- [ ] Test job submission
- [ ] Test job cancellation
- [ ] Test artifact download
- [ ] Monitor CloudWatch logs
- [ ] Set up alerts for failed jobs

---

**Ostatnia aktualizacja:** 2025-12-27

