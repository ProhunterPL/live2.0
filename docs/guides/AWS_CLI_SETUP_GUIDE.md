---
date: 2025-12-27
label: [guide, aws, cli]
---

# AWS Setup przez CLI - Kompletny Przewodnik

**Wszystkie zasoby AWS można utworzyć przez AWS CLI!**

---

## ✅ Co Można Utworzyć przez CLI

- ✅ S3 bucket i lifecycle policy
- ✅ ECR repository
- ✅ IAM Role i policies
- ✅ IAM User i policies
- ✅ AWS Batch Compute Environment
- ✅ AWS Batch Job Queue
- ✅ AWS Batch Job Definition

---

## 🚀 Szybki Start

Uruchom kompleksowy skrypt:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/deploy/setup_aws_complete_cli.ps1
```

**Lub ręcznie krok po kroku:**

---

## 📋 Krok po Kroku (CLI)

### 1. S3 Lifecycle Policy

```powershell
# Utwórz plik lifecycle.json
$lifecycle = @{
    Rules = @(
        @{
            ID = "DeleteOldArtifacts"
            Status = "Enabled"
            Expiration = @{
                Days = 90
            }
        }
    )
} | ConvertTo-Json -Depth 10 -Compress

$lifecycle | Out-File -FilePath lifecycle.json -Encoding utf8

aws s3api put-bucket-lifecycle-configuration `
    --bucket live2-artifacts `
    --lifecycle-configuration file://lifecycle.json
```

### 2. IAM Role: Live2JobRole

```powershell
# Trust policy
$trustPolicy = @{
    Version = "2012-10-17"
    Statement = @(
        @{
            Effect = "Allow"
            Principal = @{
                Service = "ecs-tasks.amazonaws.com"
            }
            Action = "sts:AssumeRole"
        }
    )
} | ConvertTo-Json -Depth 10 -Compress

$trustPolicy | Out-File -FilePath trust-policy.json -Encoding utf8

# Create role
aws iam create-role `
    --role-name Live2JobRole `
    --assume-role-policy-document file://trust-policy.json

# Attach policy
$s3Policy = @{
    Version = "2012-10-17"
    Statement = @(
        @{
            Effect = "Allow"
            Action = @("s3:PutObject", "s3:GetObject")
            Resource = "arn:aws:s3:::live2-artifacts/*"
        },
        @{
            Effect = "Allow"
            Action = @("logs:CreateLogStream", "logs:PutLogEvents")
            Resource = "arn:aws:logs:*:*:log-group:/aws/batch/job"
        }
    )
} | ConvertTo-Json -Depth 10 -Compress

$s3Policy | Out-File -FilePath s3-policy.json -Encoding utf8

aws iam put-role-policy `
    --role-name Live2JobRole `
    --policy-name S3ArtifactsAccess `
    --policy-document file://s3-policy.json
```

### 3. IAM User Policy

```powershell
$batchPolicy = @{
    Version = "2012-10-17"
    Statement = @(
        @{
            Effect = "Allow"
            Action = @("batch:SubmitJob", "batch:DescribeJobs", "batch:CancelJob")
            Resource = "*"
        }
    )
} | ConvertTo-Json -Depth 10 -Compress

$batchPolicy | Out-File -FilePath batch-policy.json -Encoding utf8

aws iam put-user-policy `
    --user-name live2-do-orchestrator `
    --policy-name BatchAccess `
    --policy-document file://batch-policy.json
```

### 4. AWS Batch Compute Environment

```powershell
# Get VPC and subnets
$vpcId = aws ec2 describe-vpcs --filters "Name=isDefault,Values=true" --query "Vpcs[0].VpcId" --output text
$subnets = aws ec2 describe-subnets --filters "Name=vpc-id,Values=$vpcId" --query "Subnets[*].SubnetId" --output text
$securityGroup = aws ec2 describe-security-groups --filters "Name=vpc-id,Values=$vpcId" "Name=group-name,Values=default" --query "SecurityGroups[0].GroupId" --output text

# Create compute environment JSON
$computeEnv = @{
    computeEnvironmentName = "live2-compute-env"
    type = "MANAGED"
    state = "ENABLED"
    computeResources = @{
        type = "SPOT"
        minvCpus = 0
        maxvCpus = 32
        desiredvCpus = 0
        instanceTypes = @("c5.2xlarge", "c5.4xlarge", "c5.9xlarge")
        subnets = $subnets -split "`t"
        securityGroupIds = @($securityGroup)
        bidPercentage = 100
    }
} | ConvertTo-Json -Depth 10

$computeEnv | Out-File -FilePath compute-env.json -Encoding utf8

aws batch create-compute-environment --cli-input-json file://compute-env.json
```

### 5. AWS Batch Job Queue

```powershell
aws batch create-job-queue `
    --job-queue-name live2-job-queue `
    --priority 1 `
    --compute-environment-order order=1,computeEnvironment=live2-compute-env `
    --state ENABLED
```

### 6. AWS Batch Job Definition

```powershell
$accountId = aws sts get-caller-identity --query Account --output text
$region = aws configure get region

$jobDef = @{
    jobDefinitionName = "live2-simulation"
    type = "container"
    containerProperties = @{
        image = "$accountId.dkr.ecr.$region.amazonaws.com/live2-simulation:latest"
        vcpus = 8
        memory = 16384
        jobRoleArn = "arn:aws:iam::${accountId}:role/Live2JobRole"
        environment = @(
            @{ name = "S3_BUCKET"; value = "live2-artifacts" },
            @{ name = "ENV"; value = "prod" }
        )
    }
    retryStrategy = @{ attempts = 2 }
    timeout = @{ attemptDurationSeconds = 3600 }
} | ConvertTo-Json -Depth 10

$jobDef | Out-File -FilePath job-definition.json -Encoding utf8

aws batch register-job-definition --cli-input-json file://job-definition.json
```

---

## ✅ Weryfikacja

```powershell
# S3
aws s3 ls s3://live2-artifacts

# ECR
aws ecr describe-repositories --repository-names live2-simulation

# IAM Role
aws iam get-role --role-name Live2JobRole

# IAM User
aws iam get-user --user-name live2-do-orchestrator

# Batch
aws batch describe-compute-environments
aws batch describe-job-queues
aws batch describe-job-definitions --job-definition-name live2-simulation
```

---

## 🆘 Troubleshooting

### Problem: JSON parsing errors
**Rozwiązanie:** Użyj `-Compress` w `ConvertTo-Json` i zapisz jako UTF-8 bez BOM

### Problem: Permissions boundary
**Rozwiązanie:** Usuń permissions boundary dla użytkownika (wymaga uprawnień admin)

### Problem: Missing service roles
**Rozwiązanie:** Utwórz `AWSBatchServiceRole` i `aws-ec2-spot-fleet-role` przez Console lub użyj istniejących

---

**Ostatnia aktualizacja:** 2025-12-27

