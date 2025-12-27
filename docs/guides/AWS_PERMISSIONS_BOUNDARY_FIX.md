---
date: 2025-12-23
label: [guide, aws, troubleshooting]
---

# Fix: Permissions Boundary Blokuje Akcje

**Problem:** Wszystkie akcje są blokowane przez "permissions boundary"

---

## 🔴 Problem

```
An error occurred (AccessDenied) when calling the CreateBucket operation: 
User: arn:aws:iam::559089787622:user/Michal is not authorized to perform: 
s3:CreateBucket on resource: "arn:aws:s3:::live2-artifacts" because 
no permissions boundary allows the s3:CreateBucket action
```

**Przyczyna:** Użytkownik "Michal" ma ustawiony **Permissions Boundary**, który ogranicza wszystkie akcje, nawet jeśli użytkownik ma uprawnienia.

---

## ✅ Rozwiązanie

### Opcja 1: Usuń Permissions Boundary (jeśli masz dostęp)

1. AWS Console → IAM → Users → **Michal**
2. Kliknij zakładkę **Permissions boundary**
3. Kliknij **Remove permissions boundary**
4. Potwierdź

**Lub przez AWS CLI:**
```powershell
aws iam put-user-permissions-boundary --user-name Michal --permissions-boundary ""
```

### Opcja 2: Poproś Administratora AWS

Jeśli nie masz uprawnień do modyfikacji IAM:

**Wyślij administratorowi:**

```
Użytkownik: Michal (arn:aws:iam::559089787622:user/Michal)

Problem: Permissions Boundary blokuje wszystkie akcje

Akcja wymagana:
- Usuń Permissions Boundary dla użytkownika "Michal"
- Lub zaktualizuj Permissions Boundary, aby pozwalał na:
  * s3:CreateBucket, s3:PutBucketLifecycleConfiguration
  * ecr:CreateRepository
  * iam:CreateUser, iam:CreateRole, iam:PutRolePolicy, iam:CreateAccessKey
```

### Opcja 3: Utwórz Ręcznie przez Console

Jeśli nie możesz usunąć Permissions Boundary, utwórz wszystko ręcznie przez AWS Console:

#### 1. S3 Bucket

1. AWS Console → S3 → **Create bucket**
2. Name: `live2-artifacts`
3. Region: `eu-central-1`
4. **Create bucket**
5. Management → **Lifecycle** → **Create lifecycle rule**
   - Name: `DeleteOldArtifacts`
   - Expiration: **90 days**

#### 2. ECR Repository

1. AWS Console → ECR → **Repositories** → **Create repository**
2. Name: `live2-simulation`
3. Region: `eu-central-1`
4. **Create repository**
5. Zapisuj URI: `559089787622.dkr.ecr.eu-central-1.amazonaws.com/live2-simulation`

#### 3. IAM Role: Live2JobRole

1. AWS Console → IAM → **Roles** → **Create role**
2. Trusted entity: **AWS service** → **ECS** → **ECS Task**
3. Permissions: **Create custom policy**

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
5. **Create role**

#### 4. IAM User: live2-do-orchestrator

1. AWS Console → IAM → **Users** → **Create user**
2. Username: `live2-do-orchestrator`
3. **Attach policies directly** → **Create custom policy**

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

4. **Create user**
5. **Security credentials** → **Create access key**
6. **Zapisuj Access Key ID i Secret Access Key** (dla DO Droplet `.env`)

---

## ✅ Weryfikacja

Po usunięciu Permissions Boundary lub utworzeniu ręcznie, sprawdź:

```powershell
# S3
aws s3 ls s3://live2-artifacts

# ECR
aws ecr describe-repositories --repository-names live2-simulation

# IAM Role
aws iam get-role --role-name Live2JobRole

# IAM User
aws iam get-user --user-name live2-do-orchestrator
```

---

## 📋 Checklist

- [ ] Permissions Boundary usunięty (lub zaktualizowany)
- [ ] S3 bucket `live2-artifacts` utworzony
- [ ] ECR repository `live2-simulation` utworzony
- [ ] IAM Role `Live2JobRole` utworzony
- [ ] IAM User `live2-do-orchestrator` utworzony
- [ ] Access Keys dla `live2-do-orchestrator` zapisane

---

**Ostatnia aktualizacja:** 2025-12-23

