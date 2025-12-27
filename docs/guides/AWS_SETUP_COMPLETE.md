---
date: 2025-12-27
label: [guide, aws, summary]
---

# AWS Infrastructure Setup - Status

**Podsumowanie utworzonej infrastruktury AWS**

---

## ✅ Co Zostało Utworzone

### 1. S3 Bucket
- **Name:** `live2-artifacts`
- **Region:** `eu-central-1`
- **Status:** ✅ Utworzony
- **Lifecycle:** ⚠️ Wymaga ręcznej konfiguracji (problem z JSON)

### 2. ECR Repository
- **Name:** `live2-simulation`
- **URI:** `559089787622.dkr.ecr.eu-central-1.amazonaws.com/live2-simulation`
- **Region:** `eu-central-1`
- **Status:** ✅ Utworzony

### 3. IAM User
- **Name:** `live2-do-orchestrator`
- **Status:** ✅ Utworzony
- **Access Keys:** ✅ Utworzone

**ZAPISZ TE CREDENTIALS:**
```
AWS_ACCESS_KEY_ID=<YOUR_ACCESS_KEY_ID>
AWS_SECRET_ACCESS_KEY=<YOUR_SECRET_ACCESS_KEY>
```

**⚠️ WAŻNE:** 
- Credentials są dostępne w AWS Console → IAM → Users → `live2-do-orchestrator` → Security credentials
- Dodaj te credentials do `.env` na DO Droplet!

### 4. IAM Role
- **Name:** `Live2JobRole`
- **Status:** ⚠️ Wymaga ręcznej konfiguracji (problem z JSON)

---

## ⚠️ Co Wymaga Ręcznej Konfiguracji

### 1. S3 Lifecycle Policy

**Przez AWS Console:**
1. S3 → `live2-artifacts` → Management → Lifecycle rules
2. Create lifecycle rule
3. Name: `DeleteOldArtifacts`
4. Expiration: 90 days
5. Apply to all objects

**Lub przez AWS CLI:**
```bash
# Utwórz plik lifecycle.json lokalnie
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

### 2. IAM Role: Live2JobRole

**Przez AWS Console:**
1. IAM → Roles → Create role
2. Trusted entity: AWS service → ECS → ECS Task
3. Permissions: Create custom policy

**Policy JSON:**
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
      "Resource": "arn:aws:s3:::live2-artifacts/*"
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

4. Role name: `Live2JobRole`
5. Create role

### 3. IAM User Policy: live2-do-orchestrator

**Przez AWS Console:**
1. IAM → Users → `live2-do-orchestrator` → Permissions
2. Add permissions → Create inline policy

**Policy JSON:**
```json
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
```

3. Policy name: `BatchAccess`
4. Create policy

---

## 📋 Następne Kroki

### 1. Dodaj AWS Credentials do DO Droplet

Na Twoim DO Droplet, edytuj `/opt/live2.0/.env`:

```bash
# AWS (dla Batch)
# Pobierz credentials z AWS Console → IAM → Users → live2-do-orchestrator → Security credentials
AWS_ACCESS_KEY_ID=<YOUR_ACCESS_KEY_ID>
AWS_SECRET_ACCESS_KEY=<YOUR_SECRET_ACCESS_KEY>
AWS_REGION=eu-central-1
AWS_BATCH_JOB_QUEUE=live2-job-queue
AWS_BATCH_JOB_DEFINITION=live2-simulation
```

### 2. Build & Push Docker Image

```bash
# Lokalnie
docker build -t live2-simulation:latest -f docker/simulation.Dockerfile .

# Login do ECR
aws ecr get-login-password --region eu-central-1 | \
  docker login --username AWS --password-stdin \
  559089787622.dkr.ecr.eu-central-1.amazonaws.com

# Tag i push
docker tag live2-simulation:latest \
  559089787622.dkr.ecr.eu-central-1.amazonaws.com/live2-simulation:latest

docker push \
  559089787622.dkr.ecr.eu-central-1.amazonaws.com/live2-simulation:latest
```

### 3. Utwórz AWS Batch

**Compute Environment:**
1. AWS Console → Batch → Compute environments → Create
2. Name: `live2-compute-env`
3. Type: Managed, Spot
4. Instance types: `c5.2xlarge`, `c5.4xlarge`
5. **Min vCPUs: 0** (scale to zero!)
6. Max vCPUs: 32
7. Create

**Job Queue:**
1. Batch → Job queues → Create
2. Name: `live2-job-queue`
3. Compute environment: `live2-compute-env`
4. Create

**Job Definition:**
1. Batch → Job definitions → Create
2. Name: `live2-simulation`
3. Image: `559089787622.dkr.ecr.eu-central-1.amazonaws.com/live2-simulation:latest`
4. Job role: `Live2JobRole`
5. vCPUs: 8, Memory: 16384
6. Retry: 2 attempts
7. Timeout: 3600 seconds
8. Create

---

## ✅ Checklist

- [x] S3 bucket utworzony
- [x] ECR repository utworzony
- [x] IAM user utworzony
- [x] Access keys utworzone
- [ ] S3 lifecycle policy skonfigurowana (ręcznie)
- [ ] IAM Role Live2JobRole skonfigurowana (ręcznie)
- [ ] IAM User policy skonfigurowana (ręcznie)
- [ ] Docker image zbudowany i wypushowany
- [ ] AWS Batch compute environment utworzony
- [ ] AWS Batch job queue utworzona
- [ ] AWS Batch job definition utworzona
- [ ] AWS credentials dodane do DO Droplet `.env`

---

**Ostatnia aktualizacja:** 2025-12-27

