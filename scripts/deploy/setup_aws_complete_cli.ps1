# Complete AWS Setup via CLI
# This script creates all AWS resources needed for Split Deploy
# Run: powershell -ExecutionPolicy Bypass -File setup_aws_complete_cli.ps1

$ErrorActionPreference = "Stop"

Write-Host "========================================" -ForegroundColor Green
Write-Host "Complete AWS Setup via CLI" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

# Get AWS Account ID and Region
$accountId = aws sts get-caller-identity --query Account --output text
$region = aws configure get region
if (-not $region) {
    $region = "eu-central-1"
    Write-Host "Using default region: $region" -ForegroundColor Yellow
}

Write-Host "Account ID: $accountId" -ForegroundColor Cyan
Write-Host "Region: $region" -ForegroundColor Cyan
Write-Host ""

# Helper function to write JSON without BOM
function Write-JsonFile {
    param(
        [string]$Path,
        [object]$Object
    )
    $json = $Object | ConvertTo-Json -Depth 10 -Compress
    $json = $json -replace "`r`n", "`n" -replace "`r", "`n"
    [System.IO.File]::WriteAllText($Path, $json, [System.Text.UTF8Encoding]::new($false))
}

# ============================================
# 1. S3 Lifecycle Policy
# ============================================
Write-Host "1. Setting S3 lifecycle policy..." -ForegroundColor Cyan

# S3 lifecycle requires JSON format (not XML for CLI)
$lifecyclePolicy = @{
    Rules = @(
        @{
            ID = "DeleteOldArtifacts"
            Status = "Enabled"
            Expiration = @{
                Days = 90
            }
        }
    )
}

$lifecycleFile = "$env:TEMP\lifecycle.json"
Write-JsonFile -Path $lifecycleFile -Object $lifecyclePolicy

try {
    aws s3api put-bucket-lifecycle-configuration `
        --bucket live2-artifacts `
        --lifecycle-configuration "file://$lifecycleFile" `
        --region $region 2>&1 | Out-Null
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "   SUCCESS: S3 lifecycle policy set" -ForegroundColor Green
    } else {
        Write-Host "   WARNING: Failed to set lifecycle policy. You can set it manually in Console." -ForegroundColor Yellow
    }
} catch {
    Write-Host "   WARNING: Failed to set lifecycle policy. You can set it manually in Console." -ForegroundColor Yellow
}

Write-Host ""

# ============================================
# 2. IAM Role: Live2JobRole
# ============================================
Write-Host "2. Creating IAM Role: Live2JobRole..." -ForegroundColor Cyan

# Check if role exists
$ErrorActionPreference = "SilentlyContinue"
$null = aws iam get-role --role-name Live2JobRole 2>&1
$roleExists = $LASTEXITCODE -eq 0
$ErrorActionPreference = "Stop"

if ($roleExists) {
    Write-Host "   Role already exists, skipping creation" -ForegroundColor Yellow
} else {
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
    }

    $trustFile = "$env:TEMP\trust-policy.json"
    Write-JsonFile -Path $trustFile -Object $trustPolicy

    try {
        aws iam create-role `
            --role-name Live2JobRole `
            --assume-role-policy-document "file://$trustFile" 2>&1 | Out-Null
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "   SUCCESS: IAM Role created" -ForegroundColor Green
        } else {
            Write-Host "   WARNING: Failed to create role. It may already exist." -ForegroundColor Yellow
        }
    } catch {
        Write-Host "   WARNING: Failed to create role. It may already exist." -ForegroundColor Yellow
    }
}

# Attach policy to role
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
}

$s3PolicyFile = "$env:TEMP\s3-policy.json"
Write-JsonFile -Path $s3PolicyFile -Object $s3Policy

aws iam put-role-policy `
    --role-name Live2JobRole `
    --policy-name S3ArtifactsAccess `
    --policy-document "file://$s3PolicyFile"

Write-Host "   SUCCESS: IAM Role policy attached" -ForegroundColor Green
Write-Host ""

# ============================================
# 3. IAM User Policy: live2-do-orchestrator
# ============================================
Write-Host "3. Creating IAM User policy..." -ForegroundColor Cyan

$batchPolicy = @{
    Version = "2012-10-17"
    Statement = @(
        @{
            Effect = "Allow"
            Action = @("batch:SubmitJob", "batch:DescribeJobs", "batch:CancelJob")
            Resource = "*"
        }
    )
}

$batchPolicyFile = "$env:TEMP\batch-policy.json"
Write-JsonFile -Path $batchPolicyFile -Object $batchPolicy

aws iam put-user-policy `
    --user-name live2-do-orchestrator `
    --policy-name BatchAccess `
    --policy-document "file://$batchPolicyFile"

Write-Host "   SUCCESS: IAM User policy attached" -ForegroundColor Green
Write-Host ""

# ============================================
# 4. AWS Batch Service Roles (if needed)
# ============================================
Write-Host "4. Checking AWS Batch service roles..." -ForegroundColor Cyan

# Check AWSBatchServiceRole
$ErrorActionPreference = "SilentlyContinue"
$null = aws iam get-role --role-name AWSBatchServiceRole 2>&1
$batchServiceRoleExists = $LASTEXITCODE -eq 0
$ErrorActionPreference = "Stop"

if (-not $batchServiceRoleExists) {
    Write-Host "   WARNING: AWSBatchServiceRole not found" -ForegroundColor Yellow
    Write-Host "   You may need to create it manually or use existing service role" -ForegroundColor Yellow
} else {
    Write-Host "   SUCCESS: AWSBatchServiceRole exists" -ForegroundColor Green
}

# Check aws-ec2-spot-fleet-role
$ErrorActionPreference = "SilentlyContinue"
$null = aws iam get-role --role-name aws-ec2-spot-fleet-role 2>&1
$spotFleetRoleExists = $LASTEXITCODE -eq 0
$ErrorActionPreference = "Stop"

if (-not $spotFleetRoleExists) {
    Write-Host "   WARNING: aws-ec2-spot-fleet-role not found" -ForegroundColor Yellow
    Write-Host "   You may need to create it manually" -ForegroundColor Yellow
} else {
    Write-Host "   SUCCESS: aws-ec2-spot-fleet-role exists" -ForegroundColor Green
}

Write-Host ""

# ============================================
# 5. AWS Batch Compute Environment
# ============================================
Write-Host "5. Creating AWS Batch Compute Environment..." -ForegroundColor Cyan

# Check if exists
$ErrorActionPreference = "SilentlyContinue"
$null = aws batch describe-compute-environments --compute-environments live2-compute-env 2>&1
$computeEnvExists = $LASTEXITCODE -eq 0
$ErrorActionPreference = "Stop"

if ($computeEnvExists) {
    Write-Host "   Compute environment already exists, skipping" -ForegroundColor Yellow
} else {
    # Get default VPC and subnets
    $vpcId = aws ec2 describe-vpcs --filters "Name=isDefault,Values=true" --query "Vpcs[0].VpcId" --output text
    $subnets = aws ec2 describe-subnets --filters "Name=vpc-id,Values=$vpcId" --query "Subnets[*].SubnetId" --output text
    
    if (-not $subnets) {
        Write-Host "   WARNING: No subnets found. You may need to specify subnets manually." -ForegroundColor Yellow
        $subnetList = @()
    } else {
        $subnetList = $subnets -split "`t"
    }

    # Get default security group
    $securityGroup = aws ec2 describe-security-groups --filters "Name=vpc-id,Values=$vpcId" "Name=group-name,Values=default" --query "SecurityGroups[0].GroupId" --output text

    $computeResources = @{
        type = "SPOT"
        minvCpus = 0
        maxvCpus = 32
        desiredvCpus = 0
        instanceTypes = @("c5.2xlarge", "c5.4xlarge", "c5.9xlarge")
        bidPercentage = 100
    }

    if ($subnetList.Count -gt 0) {
        $computeResources.subnets = $subnetList
    }
    if ($securityGroup) {
        $computeResources.securityGroupIds = @($securityGroup)
    }

    # Try to get spot fleet role ARN
    $ErrorActionPreference = "SilentlyContinue"
    $spotFleetRoleArn = aws iam get-role --role-name aws-ec2-spot-fleet-role --query "Role.Arn" --output text 2>&1
    $ErrorActionPreference = "Stop"
    if ($spotFleetRoleArn -and $LASTEXITCODE -eq 0 -and $spotFleetRoleArn -ne "None") {
        $computeResources.spotIamFleetRole = $spotFleetRoleArn
    }

    $computeEnv = @{
        computeEnvironmentName = "live2-compute-env"
        type = "MANAGED"
        state = "ENABLED"
        computeResources = $computeResources
    }

    # Try to get service role ARN
    $ErrorActionPreference = "SilentlyContinue"
    $serviceRoleArn = aws iam get-role --role-name AWSBatchServiceRole --query "Role.Arn" --output text 2>&1
    $ErrorActionPreference = "Stop"
    if ($serviceRoleArn -and $LASTEXITCODE -eq 0 -and $serviceRoleArn -ne "None") {
        $computeEnv.serviceRole = $serviceRoleArn
    }

    $computeEnvFile = "$env:TEMP\compute-env.json"
    Write-JsonFile -Path $computeEnvFile -Object $computeEnv

    aws batch create-compute-environment --cli-input-json "file://$computeEnvFile"
    Write-Host "   SUCCESS: Compute environment created" -ForegroundColor Green
}

Write-Host ""

# ============================================
# 6. AWS Batch Job Queue
# ============================================
Write-Host "6. Creating AWS Batch Job Queue..." -ForegroundColor Cyan

$ErrorActionPreference = "SilentlyContinue"
$null = aws batch describe-job-queues --job-queues live2-job-queue 2>&1
$queueExists = $LASTEXITCODE -eq 0
$ErrorActionPreference = "Stop"

if ($queueExists) {
    Write-Host "   Job queue already exists, skipping" -ForegroundColor Yellow
} else {
    $computeEnvOrder = @(
        @{
            order = 1
            computeEnvironment = "live2-compute-env"
        }
    )

    $queueOrderFile = "$env:TEMP\queue-order.json"
    Write-JsonFile -Path $queueOrderFile -Object $computeEnvOrder

    aws batch create-job-queue `
        --job-queue-name live2-job-queue `
        --priority 1 `
        --compute-environment-order "file://$queueOrderFile" `
        --state ENABLED

    Write-Host "   SUCCESS: Job queue created" -ForegroundColor Green
}

Write-Host ""

# ============================================
# 7. AWS Batch Job Definition
# ============================================
Write-Host "7. Creating AWS Batch Job Definition..." -ForegroundColor Cyan

$ecrUri = "$accountId.dkr.ecr.$region.amazonaws.com/live2-simulation"
$jobRoleArn = "arn:aws:iam::${accountId}:role/Live2JobRole"

$containerProps = @{
    image = "$ecrUri:latest"
    vcpus = 8
    memory = 16384
    jobRoleArn = $jobRoleArn
    environment = @(
        @{
            name = "S3_BUCKET"  # Environment variable name in container (kept as S3_BUCKET for compatibility)
            value = "live2-artifacts"
        },
        @{
            name = "ENV"
            value = "prod"
        },
        @{
            name = "AWS_REGION"
            value = $region
        }
    )
}

$jobDef = @{
    jobDefinitionName = "live2-simulation"
    type = "container"
    containerProperties = $containerProps
    retryStrategy = @{
        attempts = 2
    }
    timeout = @{
        attemptDurationSeconds = 3600
    }
}

$jobDefFile = "$env:TEMP\job-definition.json"
Write-JsonFile -Path $jobDefFile -Object $jobDef

aws batch register-job-definition --cli-input-json "file://$jobDefFile"

Write-Host "   SUCCESS: Job definition created" -ForegroundColor Green
Write-Host ""

# ============================================
# Summary
# ============================================
Write-Host "========================================" -ForegroundColor Green
Write-Host "AWS Setup Complete!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Created/Verified:" -ForegroundColor Cyan
Write-Host "  - S3 lifecycle policy" -ForegroundColor Green
Write-Host "  - IAM Role: Live2JobRole" -ForegroundColor Green
Write-Host "  - IAM User policy: live2-do-orchestrator" -ForegroundColor Green
Write-Host "  - Batch Compute Environment: live2-compute-env" -ForegroundColor Green
Write-Host "  - Batch Job Queue: live2-job-queue" -ForegroundColor Green
Write-Host "  - Batch Job Definition: live2-simulation" -ForegroundColor Green
Write-Host ""
Write-Host "Next Steps:" -ForegroundColor Cyan
Write-Host "1. Build and push Docker image to ECR" -ForegroundColor Yellow
Write-Host "2. Add AWS credentials to DO Droplet .env file" -ForegroundColor Yellow
Write-Host "3. Test job submission" -ForegroundColor Yellow
Write-Host ""
Write-Host "Remember: Batch with minvCpus=0 costs NOTHING when idle!" -ForegroundColor Green

