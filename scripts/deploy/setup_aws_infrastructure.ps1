# Setup AWS Infrastructure for Split Deploy
# This script creates S3, ECR, IAM, and Batch infrastructure
# Run: powershell -ExecutionPolicy Bypass -File setup_aws_infrastructure.ps1

Write-Host "Setting up AWS Infrastructure for Live 2.0..." -ForegroundColor Green
Write-Host ""

# Check AWS CLI
$awsInstalled = Get-Command aws -ErrorAction SilentlyContinue
if (-not $awsInstalled) {
    Write-Host "ERROR: AWS CLI not found. Please install it first." -ForegroundColor Red
    exit 1
}

# Get AWS Account ID and Region
Write-Host "Getting AWS Account ID..." -ForegroundColor Cyan
$accountId = aws sts get-caller-identity --query Account --output text
$region = aws configure get region
if (-not $region) {
    $region = "us-east-1"
    Write-Host "No region configured, using default: $region" -ForegroundColor Yellow
}

Write-Host "Account ID: $accountId" -ForegroundColor Green
Write-Host "Region: $region" -ForegroundColor Green
Write-Host ""

# 1. S3 Bucket
Write-Host "1. Creating S3 bucket..." -ForegroundColor Cyan
$bucketName = "live2-artifacts"
$bucketExists = aws s3 ls "s3://$bucketName" 2>$null
if ($LASTEXITCODE -eq 0) {
    Write-Host "   S3 bucket already exists: $bucketName" -ForegroundColor Yellow
} else {
    aws s3 mb "s3://$bucketName" --region $region
    Write-Host "   SUCCESS: S3 bucket created: $bucketName" -ForegroundColor Green
}

# Lifecycle policy
Write-Host "   Setting lifecycle policy (delete after 90 days)..." -ForegroundColor Yellow
$lifecyclePolicy = @{
    Rules = @(
        @{
            Id = "DeleteOldArtifacts"
            Status = "Enabled"
            Expiration = @{
                Days = 90
            }
        }
    )
}

# Write JSON without BOM, compressed, and replace Windows line endings
$json = ($lifecyclePolicy | ConvertTo-Json -Depth 10 -Compress)
$json = $json -replace "`r`n", "`n" -replace "`r", "`n"
[System.IO.File]::WriteAllText("$env:TEMP\lifecycle.json", $json, [System.Text.UTF8Encoding]::new($false))
aws s3api put-bucket-lifecycle-configuration --bucket $bucketName --lifecycle-configuration "file://$env:TEMP\lifecycle.json"
Write-Host "   SUCCESS: Lifecycle policy set" -ForegroundColor Green

# 2. ECR Repository
Write-Host ""
Write-Host "2. Creating ECR repository..." -ForegroundColor Cyan
$ecrRepo = "live2-simulation"
$ecrExists = aws ecr describe-repositories --repository-names $ecrRepo --region $region 2>$null
if ($LASTEXITCODE -eq 0) {
    Write-Host "   ECR repository already exists: $ecrRepo" -ForegroundColor Yellow
} else {
    aws ecr create-repository --repository-name $ecrRepo --region $region
    Write-Host "   SUCCESS: ECR repository created: $ecrRepo" -ForegroundColor Green
}

$ecrUri = "$accountId.dkr.ecr.$region.amazonaws.com/$ecrRepo"
Write-Host "   ECR URI: $ecrUri" -ForegroundColor Green

# 3. IAM Role for Job Execution
Write-Host ""
Write-Host "3. Creating IAM role for job execution..." -ForegroundColor Cyan

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
} | ConvertTo-Json -Depth 10

# Write JSON without BOM, compressed, and fix line endings
$json = ($trustPolicy | ConvertTo-Json -Depth 10 -Compress)
$json = $json -replace "`r`n", "`n" -replace "`r", "`n"
[System.IO.File]::WriteAllText("$env:TEMP\trust-policy.json", $json, [System.Text.UTF8Encoding]::new($false))

# Check if role exists
$roleExists = aws iam get-role --role-name Live2JobRole 2>$null
if ($LASTEXITCODE -eq 0) {
    Write-Host "   IAM role already exists: Live2JobRole" -ForegroundColor Yellow
} else {
    aws iam create-role --role-name Live2JobRole --assume-role-policy-document "file://$env:TEMP\trust-policy.json"
    Write-Host "   SUCCESS: IAM role created: Live2JobRole" -ForegroundColor Green
}

# S3 and CloudWatch policy
$s3Policy = @{
    Version = "2012-10-17"
    Statement = @(
        @{
            Effect = "Allow"
            Action = @("s3:PutObject", "s3:GetObject")
            Resource = "arn:aws:s3:::$bucketName/*"
        },
        @{
            Effect = "Allow"
            Action = @("logs:CreateLogStream", "logs:PutLogEvents")
            Resource = "arn:aws:logs:*:*:log-group:/aws/batch/job"
        }
    )
} | ConvertTo-Json -Depth 10

# Write JSON without BOM, compressed, and fix line endings
$json = ($s3Policy | ConvertTo-Json -Depth 10 -Compress)
$json = $json -replace "`r`n", "`n" -replace "`r", "`n"
[System.IO.File]::WriteAllText("$env:TEMP\s3-policy.json", $json, [System.Text.UTF8Encoding]::new($false))
aws iam put-role-policy --role-name Live2JobRole --policy-name S3ArtifactsAccess --policy-document "file://$env:TEMP\s3-policy.json"
Write-Host "   SUCCESS: S3 and CloudWatch policies attached" -ForegroundColor Green

# 4. IAM User for DO Orchestrator
Write-Host ""
Write-Host "4. Creating IAM user for DO orchestrator..." -ForegroundColor Cyan

$userExists = aws iam get-user --user-name live2-do-orchestrator 2>$null
if ($LASTEXITCODE -eq 0) {
    Write-Host "   IAM user already exists: live2-do-orchestrator" -ForegroundColor Yellow
} else {
    aws iam create-user --user-name live2-do-orchestrator
    Write-Host "   SUCCESS: IAM user created: live2-do-orchestrator" -ForegroundColor Green
}

# Batch access policy
$batchPolicy = @{
    Version = "2012-10-17"
    Statement = @(
        @{
            Effect = "Allow"
            Action = @("batch:SubmitJob", "batch:DescribeJobs", "batch:CancelJob")
            Resource = "*"
        }
    )
} | ConvertTo-Json -Depth 10

# Write JSON without BOM, compressed, and fix line endings
$json = ($batchPolicy | ConvertTo-Json -Depth 10 -Compress)
$json = $json -replace "`r`n", "`n" -replace "`r", "`n"
[System.IO.File]::WriteAllText("$env:TEMP\batch-policy.json", $json, [System.Text.UTF8Encoding]::new($false))
aws iam put-user-policy --user-name live2-do-orchestrator --policy-name BatchAccess --policy-document "file://$env:TEMP\batch-policy.json"
Write-Host "   SUCCESS: Batch access policy attached" -ForegroundColor Green

# Check for existing access keys
Write-Host ""
Write-Host "   Checking for existing access keys..." -ForegroundColor Yellow
$existingKeys = aws iam list-access-keys --user-name live2-do-orchestrator --query 'AccessKeyMetadata[*].AccessKeyId' --output text
if ($existingKeys) {
    Write-Host "   WARNING: Access keys already exist for this user" -ForegroundColor Yellow
    Write-Host "   Existing keys: $existingKeys" -ForegroundColor Yellow
    Write-Host "   If you need new keys, delete old ones first or use existing ones" -ForegroundColor Yellow
} else {
    Write-Host "   Creating new access key..." -ForegroundColor Yellow
    $newKey = aws iam create-access-key --user-name live2-do-orchestrator
    Write-Host ""
    Write-Host "   ========================================" -ForegroundColor Red
    Write-Host "   IMPORTANT: Save these credentials!" -ForegroundColor Red
    Write-Host "   ========================================" -ForegroundColor Red
    Write-Host $newKey -ForegroundColor Yellow
    Write-Host "   ========================================" -ForegroundColor Red
    Write-Host "   Add these to your DO Droplet .env file:" -ForegroundColor Yellow
    Write-Host "   AWS_ACCESS_KEY_ID=<AccessKeyId>" -ForegroundColor Yellow
    Write-Host "   AWS_SECRET_ACCESS_KEY=<SecretAccessKey>" -ForegroundColor Yellow
    Write-Host ""
}

# Summary
Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host "AWS Infrastructure Setup Complete!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Created/Verified:" -ForegroundColor Cyan
Write-Host "  - S3 Bucket: $bucketName" -ForegroundColor Green
Write-Host "  - ECR Repository: $ecrRepo" -ForegroundColor Green
Write-Host "  - ECR URI: $ecrUri" -ForegroundColor Green
Write-Host "  - IAM Role: Live2JobRole" -ForegroundColor Green
Write-Host "  - IAM User: live2-do-orchestrator" -ForegroundColor Green
Write-Host ""
Write-Host "Next Steps:" -ForegroundColor Cyan
Write-Host "1. Build and push Docker image to ECR" -ForegroundColor Yellow
Write-Host "2. Create Batch Compute Environment (minvCpus=0 for zero cost)" -ForegroundColor Yellow
Write-Host "3. Create Batch Job Queue" -ForegroundColor Yellow
Write-Host "4. Create Batch Job Definition" -ForegroundColor Yellow
Write-Host "5. Add AWS credentials to DO Droplet .env file" -ForegroundColor Yellow
Write-Host ""
Write-Host "Remember: Batch with minvCpus=0 costs NOTHING when idle!" -ForegroundColor Green

