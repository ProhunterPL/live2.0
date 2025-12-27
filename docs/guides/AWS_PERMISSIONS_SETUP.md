---
date: 2025-12-23
label: [guide, aws, permissions]
---

# AWS Permissions Setup

**Jak dodać uprawnienia dla użytkownika AWS do tworzenia infrastruktury**

---

## 🔐 Problem

Użytkownik AWS nie ma uprawnień do:
- S3: CreateBucket
- ECR: CreateRepository
- IAM: CreateUser, CreateRole, PutRolePolicy

---

## ✅ Rozwiązanie: Dodaj Uprawnienia

### Opcja 1: Przez AWS Console (Rekomendowane)

#### 1. Zaloguj się do AWS Console
https://console.aws.amazon.com

#### 2. IAM → Users → Michal → Add permissions

**Attach policies directly:**

Wybierz następujące managed policies:
- ✅ **AmazonS3FullAccess** (lub bardziej restrykcyjna: `AmazonS3ReadWriteAccess`)
- ✅ **AmazonEC2ContainerRegistryFullAccess** (lub `AmazonEC2ContainerRegistryPowerUser`)
- ✅ **IAMFullAccess** (lub bardziej restrykcyjna: `IAMReadOnlyAccess` + custom policy)

**Lub utwórz custom policy:**

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:CreateBucket",
        "s3:DeleteBucket",
        "s3:ListBucket",
        "s3:GetBucketLocation",
        "s3:PutBucketLifecycleConfiguration",
        "s3:GetBucketLifecycleConfiguration",
        "s3:PutObject",
        "s3:GetObject",
        "s3:DeleteObject"
      ],
      "Resource": [
        "arn:aws:s3:::live2-artifacts",
        "arn:aws:s3:::live2-artifacts/*"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
        "ecr:CreateRepository",
        "ecr:DescribeRepositories",
        "ecr:ListImages",
        "ecr:PutImage",
        "ecr:GetAuthorizationToken"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "iam:CreateUser",
        "iam:GetUser",
        "iam:ListUsers",
        "iam:CreateRole",
        "iam:GetRole",
        "iam:PutRolePolicy",
        "iam:AttachRolePolicy",
        "iam:CreateAccessKey",
        "iam:ListAccessKeys",
        "iam:DeleteAccessKey"
      ],
      "Resource": [
        "arn:aws:iam::559089787622:user/live2-do-orchestrator",
        "arn:aws:iam::559089787622:role/Live2JobRole"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
        "sts:GetCallerIdentity"
      ],
      "Resource": "*"
    }
  ]
}
```

#### 3. Usuń Permissions Boundary (jeśli istnieje)

Jeśli widzisz błąd "no permissions boundary allows", musisz:
1. IAM → Users → Michal → Permissions boundary
2. Remove permissions boundary (lub poproś administratora)

---

### Opcja 2: Przez AWS CLI (jeśli masz admin access)

```powershell
# Attach managed policies
aws iam attach-user-policy --user-name Michal --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess
aws iam attach-user-policy --user-name Michal --policy-arn arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryFullAccess
aws iam attach-user-policy --user-name Michal --policy-arn arn:aws:iam::aws:policy/IAMFullAccess
```

---

### Opcja 3: Poproś Administratora AWS

Jeśli nie masz uprawnień do modyfikacji IAM:

**Wyślij administratorowi AWS następujące uprawnienia:**

```
Użytkownik: Michal (arn:aws:iam::559089787622:user/Michal)

Wymagane uprawnienia:
1. S3:
   - s3:CreateBucket
   - s3:PutBucketLifecycleConfiguration
   - s3:ListBucket
   - s3:GetBucketLocation

2. ECR:
   - ecr:CreateRepository
   - ecr:DescribeRepositories

3. IAM:
   - iam:CreateUser
   - iam:CreateRole
   - iam:PutRolePolicy
   - iam:CreateAccessKey
   - iam:ListAccessKeys

4. Usunięcie Permissions Boundary (jeśli istnieje)
```

---

## 🔧 Alternatywa: Ręczne Utworzenie przez Console

Jeśli nie możesz dodać uprawnień, możesz utworzyć wszystko ręcznie przez AWS Console:

### 1. S3 Bucket

1. AWS Console → S3 → Create bucket
2. Name: `live2-artifacts`
3. Region: `eu-central-1`
4. Create bucket
5. Management → Lifecycle → Create lifecycle rule
   - Name: `DeleteOldArtifacts`
   - Expiration: 90 days

### 2. ECR Repository

1. AWS Console → ECR → Repositories → Create repository
2. Name: `live2-simulation`
3. Region: `eu-central-1`
4. Create repository
5. Zapisuj URI: `559089787622.dkr.ecr.eu-central-1.amazonaws.com/live2-simulation`

### 3. IAM Role (Live2JobRole)

1. AWS Console → IAM → Roles → Create role
2. Trusted entity: AWS service → ECS → ECS Task
3. Permissions: Create custom policy:

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

### 4. IAM User (live2-do-orchestrator)

1. AWS Console → IAM → Users → Create user
2. Username: `live2-do-orchestrator`
3. Attach policies: Create custom policy:

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

4. Create user
5. Security credentials → Create access key
6. Zapisuj Access Key ID i Secret Access Key (dla DO Droplet)

---

## ✅ Weryfikacja

Po dodaniu uprawnień, uruchom ponownie:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/deploy/setup_aws_infrastructure.ps1
```

**Lub sprawdź ręcznie:**

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

## 🆘 Jeśli Nadal Nie Działa

1. **Sprawdź Permissions Boundary:**
   ```powershell
   aws iam get-user --user-name Michal --query 'User.PermissionsBoundary'
   ```
   Jeśli zwraca coś (nie null), musisz poprosić administratora o usunięcie.

2. **Sprawdź czy masz MFA wymagane:**
   - Niektóre akcje mogą wymagać MFA

3. **Sprawdź Organization SCP:**
   - Jeśli konto jest w AWS Organization, mogą być dodatkowe ograniczenia

---

**Ostatnia aktualizacja:** 2025-12-23

