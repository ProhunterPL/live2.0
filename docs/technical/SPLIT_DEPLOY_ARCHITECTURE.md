---
date: 2025-12-23
label: [architecture, deployment]
---

# Split Deploy Architecture - Live 2.0

**Architektura hybrydowa: DigitalOcean (SaaS/monetyzacja) + AWS (compute on-demand)**

---

## 📐 Architektura Systemu

### Diagram Przepływu

```
┌─────────────────────────────────────────────────────────────────┐
│                        USERS / CLIENTS                          │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                    DIGITALOCEAN (Always-On)                     │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │  Frontend (React/Vite)                                   │  │
│  │  - Landing page                                          │  │
│  │  - Dashboard                                             │  │
│  │  - Billing UI                                            │  │
│  └──────────────────────────────────────────────────────────┘  │
│                             │                                    │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │  Backend API (FastAPI)                                   │  │
│  │  ├─ /api/v1/auth/*          (JWT + API keys)            │  │
│  │  ├─ /api/v1/billing/*       (Stripe integration)        │  │
│  │  ├─ /api/v1/jobs/*          (Job management)            │  │
│  │  └─ /api/v1/simulations/*   (Status/artifacts)          │  │
│  │                                                           │  │
│  │  Components:                                             │  │
│  │  - Job Orchestrator (enqueue → AWS trigger)             │  │
│  │  - Rate Limiter (Redis)                                  │  │
│  │  - Usage Tracker (Redis + Supabase)                      │  │
│  │  - Stripe Webhook Handler                                │  │
│  └──────────────────────────────────────────────────────────┘  │
│                             │                                    │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │  Reverse Proxy (Nginx/Caddy)                             │  │
│  │  - SSL termination                                       │  │
│  │  - Static files (frontend/dist)                          │  │
│  │  - API routing                                           │  │
│  └──────────────────────────────────────────────────────────┘  │
└────────────────────────────┬────────────────────────────────────┘
                             │
        ┌────────────────────┼────────────────────┐
        │                    │                    │
        ▼                    ▼                    ▼
┌──────────────┐    ┌──────────────┐    ┌──────────────┐
│  SUPABASE    │    │    REDIS     │    │    STRIPE    │
│              │    │              │    │              │
│  Postgres:   │    │  - Rate      │    │  - Products  │
│  - users     │    │    limiting  │    │  - Subscr.   │
│  - jobs      │    │  - Cache     │    │  - Webhooks  │
│  - usage     │    │  - Queue     │    │              │
│  - artifacts │    │  - Idemp.    │    │              │
│              │    │    keys      │    │              │
│  Auth:       │    │              │    │              │
│  - JWT       │    │              │    │              │
│  - RLS       │    │              │    │              │
└──────────────┘    └──────────────┘    └──────────────┘
        │                    │
        └────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                    AWS (Compute On-Demand)                      │
│                                                                 │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │  AWS Batch / ECS                                         │  │
│  │  ┌────────────────────────────────────────────────────┐  │  │
│  │  │  Job Container (Docker)                            │  │  │
│  │  │  - Simulation runner                               │  │  │
│  │  │  - Taichi GPU kernels                              │  │  │
│  │  │  - Chemistry analysis                              │  │  │
│  │  │  - Snapshot generation                             │  │  │
│  │  │  - S3 upload (artifacts)                           │  │  │
│  │  │  - Status update (Supabase/DO callback)            │  │  │
│  │  └────────────────────────────────────────────────────┘  │  │
│  │                                                           │  │
│  │  Compute Options:                                        │  │
│  │  - Spot Instances (default, ~70% savings)               │  │
│  │  - On-Demand (Pro tier, guaranteed)                     │  │
│  │  - Auto-terminate after completion                       │  │
│  └──────────────────────────────────────────────────────────┘  │
│                             │                                    │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │  S3 Bucket: live2-artifacts                              │  │
│  │  Structure:                                              │  │
│  │  s3://live2-artifacts/                                   │  │
│  │    {env}/                                                │  │
│  │      {user_id}/                                          │  │
│  │        {job_id}/                                         │  │
│  │          snapshots/                                      │  │
│  │          results/                                        │  │
│  │          logs/                                           │  │
│  │          metadata.json                                   │  │
│  └──────────────────────────────────────────────────────────┘  │
│                             │                                    │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │  CloudWatch Logs                                         │  │
│  │  - Job execution logs                                    │  │
│  │  - Error tracking                                        │  │
│  │  - Metrics (CPU, memory, duration)                       │  │
│  └──────────────────────────────────────────────────────────┘  │
│                             │                                    │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │  ECR (Elastic Container Registry)                        │  │
│  │  - live2-simulation:latest                               │  │
│  │  - Versioned tags (v1.0.0, v1.1.0, ...)                  │  │
│  │  - Rollback capability                                    │  │
│  └──────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🏗️ Komponenty i Konfiguracja

### 1. DigitalOcean (Always-On)

#### Minimalna Instancja
- **Type:** Basic Droplet
- **Size:** 2 vCPU, 4GB RAM, 80GB SSD
- **OS:** Ubuntu 22.04 LTS
- **Estimated Cost:** ~$24/mo

#### Upgrade Path
- **Current:** 2 vCPU / 4GB RAM → $24/mo
- **Scale-up:** 4 vCPU / 8GB RAM → $48/mo (gdy >1000 req/min)
- **Scale-out:** Load Balancer + 2x droplets → $48-72/mo (HA)

#### Stack Technologiczny
```
┌─────────────────────────────────────┐
│  Application Layer                  │
│  - FastAPI (Python 3.11)           │
│  - Gunicorn + Uvicorn workers (4)  │
│  - React Frontend (static build)    │
└─────────────────────────────────────┘
             │
┌─────────────────────────────────────┐
│  Reverse Proxy                       │
│  - Nginx (lub Caddy dla auto-SSL)   │
│  - SSL (Let's Encrypt)               │
│  - Static file serving               │
└─────────────────────────────────────┘
             │
┌─────────────────────────────────────┐
│  System                              │
│  - systemd services                  │
│  - Log rotation                      │
│  - Health checks                     │
└─────────────────────────────────────┘
```

#### Integracje DO

**Supabase (Postgres):**
- Connection pooling (Supabase Pooling URL)
- RLS (Row Level Security) dla multi-tenant
- Migracje: Alembic

**Redis:**
- Redis Labs (managed) lub własny na DO
- Użycie:
  - Rate limiting: `rate_limit:{user_id}:{endpoint}`
  - Job queue: `job_queue` (list)
  - Idempotency: `idempotency:{key}` (TTL 1h)
  - Cache: `cache:{key}` (TTL 5min)

**Stripe:**
- Webhook endpoint: `https://api.live2.com/api/v1/billing/webhooks/stripe`
- Events: subscription.*, invoice.payment_*

#### CI/CD (GitHub Actions)
```yaml
# .github/workflows/deploy-do.yml
- Deploy na push do main
- Build Docker image (opcjonalnie)
- SSH deploy (rsync + systemd restart)
- Health check verification
```

---

### 2. AWS (Compute On-Demand)

#### Rekomendacja: AWS Batch (zamiast ECS/EC2)

**Uzasadnienie:**
- ✅ Automatyczne zarządzanie infrastrukturą
- ✅ Wsparcie dla Spot Instances (70% oszczędności)
- ✅ Auto-scaling (0 → N → 0)
- ✅ Integracja z CloudWatch Logs
- ✅ Job queue management (wbudowane)
- ✅ Retry logic (wbudowane)
- ✅ Cost-effective dla batch workloads

**Alternatywy (odrzucone):**
- ❌ ECS Fargate: droższe, brak Spot
- ❌ EC2 per job: wymaga własnego orchestracji
- ❌ Lambda: timeout 15min, nie dla długich symulacji

#### Konfiguracja AWS Batch

**Compute Environment:**
- **Type:** Managed
- **Instance Types:** 
  - Spot: `c5.2xlarge`, `c5.4xlarge`, `c5.9xlarge` (flexible)
  - On-demand: `c5.4xlarge` (dla Pro tier)
- **Min vCPUs:** 0 (scale to zero)
- **Max vCPUs:** 32 (można zwiększyć)
- **Spot Bid:** 100% of On-Demand (max savings)

**Job Queue:**
- **Priority:** 1 (można dodać więcej dla różnych tierów)
- **State:** ENABLED

**Job Definition:**
```json
{
  "jobDefinitionName": "live2-simulation",
  "type": "container",
  "containerProperties": {
    "image": "YOUR_ACCOUNT.dkr.ecr.REGION.amazonaws.com/live2-simulation:latest",
    "vcpus": 8,
    "memory": 16384,
    "jobRoleArn": "arn:aws:iam::ACCOUNT:role/Live2JobRole",
    "environment": [
      {"name": "SUPABASE_URL", "value": "..."},
      {"name": "SUPABASE_KEY", "value": "..."},
      {"name": "REDIS_HOST", "value": "..."},
      {"name": "S3_BUCKET", "value": "live2-artifacts"},
      {"name": "JOB_ID", "value": "{{job_id}}"},
      {"name": "USER_ID", "value": "{{user_id}}"}
    ],
    "command": ["python", "-m", "backend.sim.run_simulation", "--job-id", "{{job_id}}"]
  },
  "retryStrategy": {
    "attempts": 2
  },
  "timeout": {
    "attemptDurationSeconds": 3600  # 1h max per attempt
  }
}
```

#### S3 Bucket Structure

```
s3://live2-artifacts/
├── prod/
│   └── {user_id}/
│       └── {job_id}/
│           ├── snapshots/
│           │   ├── snapshot_0000.json
│           │   └── snapshot_1000.json
│           ├── results/
│           │   ├── molecules.json
│           │   └── reactions.json
│           ├── logs/
│           │   └── simulation.log
│           └── metadata.json
└── staging/
    └── ...
```

**S3 Lifecycle Policy:**
- Delete after 90 days (staging)
- Delete after 365 days (prod)
- Transition to Glacier after 30 days (opcjonalnie)

#### IAM Roles

**Job Execution Role (`Live2JobRole`):**
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:PutObject",
        "s3:GetObject"
      ],
      "Resource": "arn:aws:s3:::live2-artifacts/{env}/{user_id}/{job_id}/*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "logs:CreateLogStream",
        "logs:PutLogEvents"
      ],
      "Resource": "arn:aws:logs:*:*:log-group:/aws/batch/job"
    }
  ]
}
```

**DO → AWS Access:**
- **Option 1 (Rekomendowane):** IAM User z Access Key
  - User: `live2-do-orchestrator`
  - Policy: `Batch:SubmitJob`, `Batch:DescribeJobs`, `S3:GetObject` (read-only)
  - Access Key w DO env vars (encrypted)

- **Option 2:** Cross-account role (bardziej bezpieczne, ale bardziej złożone)

#### VPC / Networking

**Minimalna konfiguracja:**
- Default VPC (wystarczy dla startu)
- Security Group dla Batch:
  - Outbound: All traffic (dla S3, CloudWatch, Supabase)
  - Inbound: None (nie potrzebne)

**Upgrade (opcjonalnie):**
- Dedicated VPC z private subnets
- NAT Gateway (jeśli potrzebny dostęp do private resources)

#### ECR (Container Registry)

**Repository:**
- Name: `live2-simulation`
- Image scanning: Enabled
- Lifecycle policy: Keep last 10 images

**Build & Push:**
```bash
# Build
docker build -t live2-simulation:latest -f docker/simulation.Dockerfile .

# Tag
docker tag live2-simulation:latest \
  YOUR_ACCOUNT.dkr.ecr.REGION.amazonaws.com/live2-simulation:latest
docker tag live2-simulation:latest \
  YOUR_ACCOUNT.dkr.ecr.REGION.amazonaws.com/live2-simulation:v1.0.0

# Push
aws ecr get-login-password --region REGION | \
  docker login --username AWS --password-stdin \
  YOUR_ACCOUNT.dkr.ecr.REGION.amazonaws.com
docker push YOUR_ACCOUNT.dkr.ecr.REGION.amazonaws.com/live2-simulation:latest
```

---

## 🔌 Interfejs DO ↔ AWS

### Kontrakt API (DO Backend)

#### POST /api/v1/jobs
**Start joba**

**Request:**
```json
{
  "simulation_config": {
    "mode": "open_chemistry",
    "config": {
      "num_particles": 1000,
      "temperature": 300,
      "steps": 50000
    }
  },
  "idempotency_key": "optional-unique-key"  // dla retry safety
}
```

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "queued",
  "estimated_cost": 0.15,
  "estimated_duration_minutes": 30,
  "created_at": "2025-12-23T10:00:00Z"
}
```

**Flow:**
1. Validate user subscription (active?)
2. Check rate limits (Redis)
3. Check idempotency (Redis: `idempotency:{key}`)
4. Create job record (Supabase: `jobs` table)
5. Submit AWS Batch job
6. Store `aws_batch_job_id` in Supabase
7. Return job_id

#### GET /api/v1/jobs/{job_id}
**Status joba**

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "running",  // queued|running|succeeded|failed|canceled
  "progress": 45,  // 0-100
  "created_at": "2025-12-23T10:00:00Z",
  "started_at": "2025-12-23T10:01:00Z",
  "finished_at": null,
  "cost_estimate": 0.15,
  "cost_actual": null,
  "error": null,
  "aws_batch_job_id": "abc123-def456-ghi789"
}
```

**Source:** Supabase `jobs` table (real-time)

#### POST /api/v1/jobs/{job_id}/cancel
**Anuluj job**

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "canceled",
  "canceled_at": "2025-12-23T10:15:00Z"
}
```

**Flow:**
1. Update status w Supabase → `canceled`
2. Cancel AWS Batch job (`aws batch cancel-job`)
3. Cleanup (S3 artifacts opcjonalnie)

#### GET /api/v1/jobs/{job_id}/artifacts
**Lista artefaktów (presigned URLs)**

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "artifacts": [
    {
      "type": "snapshot",
      "s3_key": "prod/user123/job456/snapshots/snapshot_0000.json",
      "presigned_url": "https://s3.amazonaws.com/...?signature=...",
      "size_bytes": 1024000,
      "created_at": "2025-12-23T10:05:00Z"
    },
    {
      "type": "result",
      "s3_key": "prod/user123/job456/results/molecules.json",
      "presigned_url": "https://s3.amazonaws.com/...?signature=...",
      "size_bytes": 512000,
      "created_at": "2025-12-23T10:30:00Z"
    }
  ]
}
```

**Flow:**
1. Query Supabase `job_artifacts` table
2. Generate presigned URLs (AWS SDK, 1h expiry)
3. Return list

---

### Model Danych (Supabase)

#### Tabela: `jobs`

```sql
CREATE TABLE jobs (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  user_id UUID NOT NULL REFERENCES auth.users(id),
  
  -- Status
  status TEXT NOT NULL DEFAULT 'queued',
  -- queued|running|succeeded|failed|canceled
  
  -- Configuration
  params JSONB NOT NULL,  -- simulation_config
  idempotency_key TEXT UNIQUE,  -- dla deduplikacji
  
  -- Timing
  created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
  started_at TIMESTAMPTZ,
  finished_at TIMESTAMPTZ,
  
  -- AWS
  aws_batch_job_id TEXT,  -- AWS Batch job ID
  aws_batch_job_arn TEXT,
  
  -- Cost tracking
  cost_estimate DECIMAL(10, 4),  -- USD
  cost_actual DECIMAL(10, 4),    -- USD (po zakończeniu)
  
  -- Results
  error TEXT,  -- error message jeśli failed
  error_details JSONB,  -- stack trace, etc.
  
  -- Metadata
  progress INTEGER DEFAULT 0,  -- 0-100
  metadata JSONB,  -- dodatkowe info
  
  CONSTRAINT status_check CHECK (status IN ('queued', 'running', 'succeeded', 'failed', 'canceled'))
);

CREATE INDEX idx_jobs_user_id ON jobs(user_id);
CREATE INDEX idx_jobs_status ON jobs(status);
CREATE INDEX idx_jobs_created_at ON jobs(created_at DESC);
CREATE INDEX idx_jobs_aws_batch_job_id ON jobs(aws_batch_job_id);
```

#### Tabela: `job_artifacts`

```sql
CREATE TABLE job_artifacts (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  job_id UUID NOT NULL REFERENCES jobs(id) ON DELETE CASCADE,
  
  artifact_type TEXT NOT NULL,
  -- snapshot|result|log|metadata
  
  s3_key TEXT NOT NULL,  -- pełna ścieżka S3
  s3_bucket TEXT NOT NULL DEFAULT 'live2-artifacts',
  
  size_bytes BIGINT,
  checksum TEXT,  -- MD5/SHA256
  
  created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
  
  CONSTRAINT artifact_type_check CHECK (artifact_type IN ('snapshot', 'result', 'log', 'metadata'))
);

CREATE INDEX idx_job_artifacts_job_id ON job_artifacts(job_id);
CREATE INDEX idx_job_artifacts_type ON job_artifacts(artifact_type);
```

#### Tabela: `usage` (rozszerzenie istniejącej)

```sql
-- Rozszerz istniejącą tabelę usage o job tracking
ALTER TABLE usage ADD COLUMN IF NOT EXISTS jobs_count INTEGER DEFAULT 0;
ALTER TABLE usage ADD COLUMN IF NOT EXISTS jobs_running_count INTEGER DEFAULT 0;
ALTER TABLE usage ADD COLUMN IF NOT EXISTS total_job_cost DECIMAL(10, 4) DEFAULT 0;
```

#### Tabela: `billing_state` (rozszerzenie istniejącej)

```sql
-- Upewnij się, że istnieje:
-- subscription_status, stripe_customer_id, tier, quota, limits
```

#### Row Level Security (RLS)

```sql
-- Jobs: użytkownik widzi tylko swoje joby
ALTER TABLE jobs ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Users can view own jobs"
  ON jobs FOR SELECT
  USING (auth.uid() = user_id);

CREATE POLICY "Users can create own jobs"
  ON jobs FOR INSERT
  WITH CHECK (auth.uid() = user_id);

-- Job artifacts: tylko właściciel joba
ALTER TABLE job_artifacts ENABLE ROW LEVEL SECURITY;

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

---

### Trigger / Transport: DO → AWS

#### Mechanizm: Bezpośredni AWS SDK Call (Rekomendowane)

**Uzasadnienie:**
- ✅ Prosty w implementacji
- ✅ Niska latencja
- ✅ Pełna kontrola nad retry logic
- ⚠️ Wymaga IAM credentials w DO (mitigacja: encrypted env vars)

**Implementacja (DO Backend):**

```python
# backend/jobs/aws_batch.py
import boto3
from botocore.exceptions import ClientError

class AWSBatchOrchestrator:
    def __init__(self):
        self.batch = boto3.client(
            'batch',
            aws_access_key_id=os.getenv('AWS_ACCESS_KEY_ID'),
            aws_secret_access_key=os.getenv('AWS_SECRET_ACCESS_KEY'),
            region_name=os.getenv('AWS_REGION', 'us-east-1')
        )
        self.job_queue = os.getenv('AWS_BATCH_JOB_QUEUE')
        self.job_definition = os.getenv('AWS_BATCH_JOB_DEFINITION')
    
    def submit_job(
        self,
        job_id: str,
        user_id: str,
        simulation_config: dict,
        tier: str  # dla wyboru Spot vs On-Demand
    ) -> str:
        """Submit job to AWS Batch, return AWS job ID"""
        
        # Determine compute type
        use_spot = tier != 'pro'  # Pro tier = On-Demand
        
        try:
            response = self.batch.submit_job(
                jobName=f"live2-{job_id[:8]}",
                jobQueue=self.job_queue,
                jobDefinition=self.job_definition,
                parameters={
                    'job_id': job_id,
                    'user_id': user_id
                },
                containerOverrides={
                    'environment': [
                        {'name': 'JOB_ID', 'value': job_id},
                        {'name': 'USER_ID', 'value': user_id},
                        {'name': 'SIMULATION_CONFIG', 'value': json.dumps(simulation_config)},
                        {'name': 'SUPABASE_URL', 'value': os.getenv('SUPABASE_URL')},
                        {'name': 'SUPABASE_SERVICE_KEY', 'value': os.getenv('SUPABASE_SERVICE_KEY')},
                        {'name': 'S3_BUCKET', 'value': 'live2-artifacts'},
                        {'name': 'ENV', 'value': os.getenv('ENV', 'prod')}
                    ]
                },
                retryStrategy={
                    'attempts': 2 if use_spot else 1  # Spot może być przerwany
                },
                timeout={
                    'attemptDurationSeconds': 3600  # 1h max
                }
            )
            
            return response['jobId']
            
        except ClientError as e:
            logger.error(f"Failed to submit AWS Batch job: {e}")
            raise
    
    def cancel_job(self, aws_batch_job_id: str):
        """Cancel running AWS Batch job"""
        try:
            self.batch.cancel_job(
                jobId=aws_batch_job_id,
                reason="User requested cancellation"
            )
        except ClientError as e:
            logger.error(f"Failed to cancel AWS Batch job: {e}")
            raise
```

**Alternatywy (odrzucone):**
- ❌ SQS: dodaje opóźnienie, wymaga worker w AWS
- ❌ EventBridge: overkill dla prostego triggera
- ❌ Lambda trigger: dodatkowa warstwa, koszt

---

### Callback: AWS → DO

#### Opcja 1: Supabase Update (Rekomendowane)

**Uzasadnienie:**
- ✅ Prosty (job container zapisuje status do Supabase)
- ✅ Real-time (Supabase Realtime subscriptions)
- ✅ Nie wymaga webhook endpoint w DO
- ✅ Idempotent (upsert)

**Implementacja (AWS Job Container):**

```python
# W job containerze (backend/sim/run_simulation.py)
from supabase import create_client

def update_job_status(
    job_id: str,
    status: str,
    progress: int = None,
    error: str = None,
    artifacts: list = None
):
    """Update job status in Supabase"""
    
    supabase = create_client(
        os.getenv('SUPABASE_URL'),
        os.getenv('SUPABASE_SERVICE_KEY')  # Service key (bypass RLS)
    )
    
    update_data = {
        'status': status,
        'updated_at': datetime.utcnow().isoformat()
    }
    
    if progress is not None:
        update_data['progress'] = progress
    
    if error:
        update_data['error'] = error
        update_data['error_details'] = {'timestamp': datetime.utcnow().isoformat()}
    
    if status == 'succeeded':
        update_data['finished_at'] = datetime.utcnow().isoformat()
    
    # Upsert job
    supabase.table('jobs').update(update_data).eq('id', job_id).execute()
    
    # Insert artifacts
    if artifacts:
        for artifact in artifacts:
            supabase.table('job_artifacts').insert({
                'job_id': job_id,
                'artifact_type': artifact['type'],
                's3_key': artifact['s3_key'],
                'size_bytes': artifact['size_bytes'],
                'checksum': artifact.get('checksum')
            }).execute()
```

**Frontend (Real-time updates):**
```typescript
// frontend/src/hooks/useJobStatus.ts
import { RealtimeChannel } from '@supabase/supabase-js'

const channel = supabase
  .channel(`job:${jobId}`)
  .on('postgres_changes', {
    event: 'UPDATE',
    schema: 'public',
    table: 'jobs',
    filter: `id=eq.${jobId}`
  }, (payload) => {
    setJobStatus(payload.new)
  })
  .subscribe()
```

#### Opcja 2: Webhook Callback (Alternatywa)

**Uzasadnienie:**
- ✅ Szybsze powiadomienia (bez polling)
- ⚠️ Wymaga webhook endpoint w DO
- ⚠️ Wymaga HMAC signature verification

**Implementacja:**

```python
# W job containerze
import requests
import hmac
import hashlib

def notify_do_webhook(job_id: str, status: str, secret: str):
    """Send webhook to DO backend"""
    
    payload = {
        'job_id': job_id,
        'status': status,
        'timestamp': datetime.utcnow().isoformat()
    }
    
    # Generate HMAC signature
    signature = hmac.new(
        secret.encode(),
        json.dumps(payload).encode(),
        hashlib.sha256
    ).hexdigest()
    
    requests.post(
        f"{os.getenv('DO_API_BASE_URL')}/api/v1/jobs/{job_id}/complete",
        json=payload,
        headers={
            'X-Webhook-Signature': signature,
            'X-Webhook-Source': 'aws-batch'
        },
        timeout=5
    )
```

**DO Backend (webhook handler):**
```python
@router.post("/jobs/{job_id}/complete")
async def job_complete_webhook(
    job_id: str,
    payload: dict,
    request: Request
):
    """Handle AWS → DO webhook"""
    
    # Verify HMAC signature
    signature = request.headers.get('X-Webhook-Signature')
    expected_sig = hmac.new(
        WEBHOOK_SECRET.encode(),
        json.dumps(payload).encode(),
        hashlib.sha256
    ).hexdigest()
    
    if not hmac.compare_digest(signature, expected_sig):
        raise HTTPException(401, "Invalid signature")
    
    # Update Supabase (idempotent)
    # ...
```

**Rekomendacja:** Użyj Opcji 1 (Supabase) dla MVP, Opcja 2 (webhook) dla production (szybsze powiadomienia).

---

### Idempotencja / Retry / Deduplikacja

#### Idempotency Keys

**Mechanizm:**
- Client może przesłać `idempotency_key` w `POST /jobs`
- DO backend sprawdza Redis: `idempotency:{key}`
- Jeśli istnieje → zwraca istniejący `job_id`
- Jeśli nie → tworzy job, zapisuje key w Redis (TTL 1h)

**Implementacja:**
```python
async def create_job(request: JobRequest, user: User):
    # Check idempotency
    if request.idempotency_key:
        cache_key = f"idempotency:{request.idempotency_key}"
        existing_job_id = redis.get(cache_key)
        if existing_job_id:
            return {"job_id": existing_job_id, "status": "queued"}
        
        # Store for 1h
        redis.setex(cache_key, 3600, job_id)
    
    # Create job...
```

#### Retry Logic

**AWS Batch:**
- Wbudowane retry (2 attempts dla Spot, 1 dla On-Demand)
- Job container powinien być idempotent (checkpoint/resume)

**DO Backend:**
- Retry dla AWS SDK calls (exponential backoff)
- Max 3 retries

#### Deduplikacja

**Rate Limiting (Redis):**
```python
# Per user, per endpoint
rate_limit_key = f"rate_limit:jobs:{user_id}"
current = redis.incr(rate_limit_key)
if current == 1:
    redis.expire(rate_limit_key, 3600)  # 1h window

if current > MAX_JOBS_PER_HOUR:
    raise HTTPException(429, "Rate limit exceeded")
```

---

## 🛡️ Guardrails Kosztowe

### Tabela Guardrails (Wartości Startowe)

| Guardrail | Hobby | Research | Pro | Implementacja |
|-----------|-------|----------|-----|---------------|
| **Max równoległe joby** | 1 | 2 | 5 | Redis counter + check przed submit |
| **Max joby/dzień** | 5 | 20 | 100 | Redis daily counter (TTL 24h) |
| **Max runtime joba** | 30 min | 2h | 6h | AWS Batch timeout |
| **Auto-terminate** | ✅ | ✅ | ✅ | AWS Batch timeout |
| **Spot vs On-Demand** | Spot only | Spot (fallback On-Demand) | On-Demand (opcjonalnie Spot) | Job definition override |
| **Max cost/job** | $0.50 | $5.00 | $50.00 | Estimate przed submit, alert jeśli przekroczone |
| **Rate limit (API)** | 10 req/min | 50 req/min | 200 req/min | Redis rate limiter |
| **Blokada (brak subskrypcji)** | ✅ | ✅ | ✅ | Check `subscription_status` przed submit |

### Implementacja Guardrails

#### 1. Pre-Submit Validation (DO Backend)

```python
async def validate_job_request(user: User, request: JobRequest):
    """Validate job request against guardrails"""
    
    # Check subscription
    if user.subscription_status != 'active':
        raise HTTPException(403, "Active subscription required")
    
    # Get tier limits
    tier_limits = TIER_LIMITS[user.tier]
    
    # Check parallel jobs
    running_count = await get_running_jobs_count(user.id)
    if running_count >= tier_limits['max_parallel']:
        raise HTTPException(429, f"Max {tier_limits['max_parallel']} parallel jobs allowed")
    
    # Check daily limit
    daily_count = await get_daily_jobs_count(user.id)
    if daily_count >= tier_limits['max_per_day']:
        raise HTTPException(429, f"Max {tier_limits['max_per_day']} jobs per day")
    
    # Estimate cost
    estimated_cost = estimate_job_cost(request.simulation_config)
    if estimated_cost > tier_limits['max_cost_per_job']:
        raise HTTPException(400, f"Job cost exceeds limit: ${tier_limits['max_cost_per_job']}")
    
    # Check runtime
    estimated_duration = estimate_job_duration(request.simulation_config)
    max_duration = tier_limits['max_runtime_minutes']
    if estimated_duration > max_duration:
        raise HTTPException(400, f"Job duration exceeds limit: {max_duration} minutes")
    
    return True
```

#### 2. AWS Budgets & Alerts

**Budżet AWS:**
- **Monthly Budget:** $100 (start)
- **Alert Thresholds:**
  - 50% → Email alert
  - 80% → Email + Slack
  - 100% → Email + Slack + Auto-stop (opcjonalnie)

**CloudWatch Alarms:**
- High cost per job (>$10)
- Unusual job duration (>6h)
- Failed job rate (>10%)

**Implementacja:**
```bash
# AWS CLI
aws budgets create-budget \
  --account-id YOUR_ACCOUNT \
  --budget file://budget.json \
  --notifications-with-subscribers file://notifications.json
```

#### 3. Rate Limiting (Redis)

```python
# backend/middleware/rate_limiter.py
from fastapi import Request, HTTPException
import redis

redis_client = redis.Redis(...)

async def rate_limit_middleware(request: Request, user: User):
    """Rate limit per user, per endpoint"""
    
    endpoint = request.url.path
    tier_limits = TIER_LIMITS[user.tier]
    limit = tier_limits['rate_limit_per_minute']
    
    key = f"rate_limit:{user.id}:{endpoint}"
    current = redis_client.incr(key)
    
    if current == 1:
        redis_client.expire(key, 60)  # 1 minute window
    
    if current > limit:
        raise HTTPException(
            429,
            f"Rate limit exceeded: {limit} requests per minute"
        )
    
    request.state.rate_limit_remaining = limit - current
```

---

## 🚀 Plan Wdrożenia (10 Kroków)

### Krok 1: DO - Deploy API + Stripe + Webhooki
**Cel:** Uruchomić podstawowy backend na DO

**Akcje:**
1. Utwórz Droplet (2 vCPU / 4GB RAM)
2. Setup Nginx + SSL (Let's Encrypt)
3. Deploy FastAPI backend (systemd service)
4. Skonfiguruj Stripe webhook endpoint
5. Test: Register → Login → Checkout

**Weryfikacja:**
```bash
curl https://api.live2.com/status/health
curl -X POST https://api.live2.com/api/v1/auth/register ...
```

**Czas:** 2-3h

---

### Krok 2: Supabase - Schemat Tabel + RLS
**Cel:** Przygotować model danych dla jobów

**Akcje:**
1. Utwórz tabele: `jobs`, `job_artifacts`
2. Dodaj indeksy
3. Włącz RLS + utwórz policies
4. Test: Insert job (service key)

**SQL:**
```sql
-- Wykonaj migracje z sekcji "Model Danych"
```

**Weryfikacja:**
```sql
SELECT * FROM jobs LIMIT 1;
SELECT * FROM job_artifacts LIMIT 1;
```

**Czas:** 1h

---

### Krok 3: Redis - Rate Limiting + Idempotency
**Cel:** Skonfigurować Redis dla guardrails

**Akcje:**
1. Skonfiguruj Redis connection (Redis Labs lub własny)
2. Test rate limiting (per user)
3. Test idempotency keys
4. Monitor Redis memory usage

**Weryfikacja:**
```python
# Test script
redis.set("test", "value")
assert redis.get("test") == "value"
```

**Czas:** 30min

---

### Krok 4: AWS - S3 + ECR + IAM Role
**Cel:** Przygotować AWS infrastructure

**Akcje:**
1. Utwórz S3 bucket: `live2-artifacts`
2. Setup lifecycle policy (delete after 90d)
3. Utwórz ECR repository: `live2-simulation`
4. Utwórz IAM role: `Live2JobRole` (S3 + CloudWatch)
5. Utwórz IAM user: `live2-do-orchestrator` (Batch access)
6. Test: Upload test file do S3

**Weryfikacja:**
```bash
aws s3 ls s3://live2-artifacts/
aws ecr describe-repositories --repository-names live2-simulation
```

**Czas:** 1-2h

---

### Krok 5: AWS - Batch Setup (Compute Environment + Queue)
**Cel:** Skonfigurować AWS Batch

**Akcje:**
1. Utwórz Compute Environment (Managed, Spot + On-Demand)
2. Utwórz Job Queue
3. Utwórz Job Definition (placeholder image)
4. Test: Submit test job (echo "hello")

**Weryfikacja:**
```bash
aws batch describe-compute-environments
aws batch describe-job-queues
aws batch submit-job --job-name test --job-queue ...
```

**Czas:** 1-2h

---

### Krok 6: AWS - Build & Push Docker Image
**Cel:** Przygotować container image dla symulacji

**Akcje:**
1. Utwórz `docker/simulation.Dockerfile`
2. Build image lokalnie
3. Push do ECR
4. Update Job Definition (nowy image)

**Dockerfile:**
```dockerfile
FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY backend/ ./backend/
CMD ["python", "-m", "backend.sim.run_simulation"]
```

**Weryfikacja:**
```bash
docker pull YOUR_ACCOUNT.dkr.ecr.REGION.amazonaws.com/live2-simulation:latest
```

**Czas:** 1h

---

### Krok 7: DO - Start Job (Trigger AWS) + Zapis batch_job_id
**Cel:** Zintegrować DO → AWS flow

**Akcje:**
1. Implementuj `AWSBatchOrchestrator` (backend/jobs/aws_batch.py)
2. Dodaj endpoint `POST /api/v1/jobs`
3. Implementuj guardrails validation
4. Test: Submit job → sprawdź AWS Batch console

**Weryfikacja:**
```bash
curl -X POST https://api.live2.com/api/v1/jobs \
  -H "X-API-Key: ..." \
  -d '{"simulation_config": {...}}'

# Sprawdź AWS Batch
aws batch describe-jobs --jobs <job_id>
```

**Czas:** 2-3h

---

### Krok 8: AWS → DO Callback (Update Status w Supabase)
**Cel:** Job container aktualizuje status w Supabase

**Akcje:**
1. W job containerze: implementuj `update_job_status()`
2. Test: Job zapisuje status → sprawdź Supabase
3. Frontend: Real-time subscription (opcjonalnie)

**Weryfikacja:**
```sql
-- W Supabase
SELECT * FROM jobs WHERE id = '<job_id>';
-- Powinien być status 'running' → 'succeeded'
```

**Czas:** 1-2h

---

### Krok 9: Presigned URLs do Wyników
**Cel:** Umożliwić pobieranie artefaktów z S3

**Akcje:**
1. Implementuj `GET /api/v1/jobs/{job_id}/artifacts`
2. Generate presigned URLs (AWS SDK, 1h expiry)
3. Test: Pobierz snapshot/result

**Weryfikacja:**
```bash
curl https://api.live2.com/api/v1/jobs/{job_id}/artifacts \
  -H "X-API-Key: ..."

# Powinien zwrócić listę presigned URLs
```

**Czas:** 1h

---

### Krok 10: Alerty + Budżety + Testy E2E
**Cel:** Monitoring i finalne testy

**Akcje:**
1. Setup AWS Budgets ($100/mo, alerty)
2. Setup CloudWatch Alarms (cost, duration, failures)
3. Test E2E: Register → Subscribe → Create Job → Wait → Download
4. Test scenariusze błędów:
   - Retry (Spot interruption)
   - Cancel job
   - Rate limit exceeded
   - Subscription expired

**Weryfikacja:**
```bash
# E2E test script
python scripts/test_e2e_split_deploy.py
```

**Czas:** 2-3h

---

## 📊 Szacunkowe Koszty (MVP)

### DigitalOcean
- **Droplet (2 vCPU / 4GB):** $24/mo
- **Total:** ~$24/mo

### AWS (przy założeniu: 100 jobów/miesiąc, średnio 30min/job)
- **Batch (Spot, c5.4xlarge):** ~$0.10/job × 100 = $10/mo
- **S3 Storage (100GB):** ~$2.30/mo
- **CloudWatch Logs:** ~$1/mo
- **ECR:** ~$0.10/mo
- **Total:** ~$13-15/mo

### Supabase
- **Free tier:** $0/mo (do 500MB DB, 2GB bandwidth)
- **Pro tier (jeśli potrzebne):** $25/mo

### Redis
- **Redis Labs (30MB):** $0/mo (free tier)
- **Lub własny na DO:** wliczone w Droplet

### Stripe
- **Transaction fees:** 2.9% + $0.30 (per payment)

### **Total MVP:** ~$37-64/mo (stały koszt)

---

## 🔄 Plan Skalowania

### Faza 1: MVP (0-100 użytkowników)
- DO: 1 Droplet (2 vCPU / 4GB)
- AWS: Batch (max 32 vCPUs)
- Koszt: ~$40/mo

### Faza 2: Growth (100-1000 użytkowników)
- DO: Upgrade do 4 vCPU / 8GB ($48/mo) lub Load Balancer + 2x droplets
- AWS: Zwiększ max vCPUs do 128
- Koszt: ~$60-100/mo

### Faza 3: Scale (1000+ użytkowników)
- DO: Load Balancer + Auto-scaling droplets
- AWS: Dedicated VPC, Reserved Instances (opcjonalnie)
- Supabase: Pro tier
- Koszt: ~$200-500/mo

---

## ✅ Checklist Wdrożenia

### Pre-Deployment
- [ ] DO Droplet utworzony
- [ ] Domain + SSL skonfigurowane
- [ ] Supabase project + tables
- [ ] Redis dostępny
- [ ] Stripe produkty + webhook
- [ ] AWS account + S3 + ECR + Batch

### Deployment
- [ ] Backend API deployed na DO
- [ ] Frontend build deployed
- [ ] Database migrations uruchomione
- [ ] AWS Batch job definition utworzona
- [ ] Docker image w ECR
- [ ] IAM roles/permissions skonfigurowane

### Testing
- [ ] Health checks działają
- [ ] Auth flow (register/login)
- [ ] Stripe checkout flow
- [ ] Job submission (DO → AWS)
- [ ] Job execution (AWS → Supabase update)
- [ ] Artifact download (presigned URLs)
- [ ] Guardrails (rate limits, quotas)
- [ ] Error scenarios (retry, cancel, failure)

### Monitoring
- [ ] AWS Budgets skonfigurowane
- [ ] CloudWatch Alarms skonfigurowane
- [ ] Log aggregation (DO + AWS)
- [ ] Uptime monitoring (UptimeRobot)

---

**Ostatnia aktualizacja:** 2025-12-23

