# Create AWS Batch Compute Environment and Job Queue
# Run: powershell -ExecutionPolicy Bypass -File scripts/deploy/create_batch_queue.ps1

$ErrorActionPreference = "Stop"

Write-Host "========================================" -ForegroundColor Green
Write-Host "Creating AWS Batch Queue" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

# Get AWS Account ID and Region
$accountId = (aws sts get-caller-identity --query Account --output text).Trim()
$region = (aws configure get region).Trim()
if (-not $region) {
    $region = "eu-central-1"
}

Write-Host "Account ID: $accountId" -ForegroundColor Cyan
Write-Host "Region: $region" -ForegroundColor Cyan
Write-Host ""

# Function to write JSON without BOM
function Write-JsonFile {
    param(
        [string]$Path,
        [object]$Object
    )
    $json = $Object | ConvertTo-Json -Depth 10 -Compress
    [System.IO.File]::WriteAllText($Path, $json, [System.Text.UTF8Encoding]::new($false))
}

# ============================================
# 1. Check/Create Compute Environment
# ============================================
Write-Host "1. Checking Compute Environment..." -ForegroundColor Cyan

$ErrorActionPreference = "SilentlyContinue"
$null = aws batch describe-compute-environments --compute-environments live2-compute-env --region $region 2>&1
$computeEnvExists = $LASTEXITCODE -eq 0
$ErrorActionPreference = "Stop"

if ($computeEnvExists) {
    Write-Host "   Compute environment exists" -ForegroundColor Green
} else {
    Write-Host "   Creating Compute Environment..." -ForegroundColor Yellow
    
    # Get VPC and subnets
    $vpcId = aws ec2 describe-vpcs --filters "Name=isDefault,Values=true" --query "Vpcs[0].VpcId" --output text --region $region
    $subnets = aws ec2 describe-subnets --filters "Name=vpc-id,Values=$vpcId" --query "Subnets[*].SubnetId" --output text --region $region
    $securityGroup = aws ec2 describe-security-groups --filters "Name=vpc-id,Values=$vpcId" "Name=group-name,Values=default" --query "SecurityGroups[0].GroupId" --output text --region $region
    
    $subnetList = $subnets -split "`t" | Where-Object { $_ }
    
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
    
    $computeEnv = @{
        computeEnvironmentName = "live2-compute-env"
        type = "MANAGED"
        state = "ENABLED"
        computeResources = $computeResources
    }
    
    $computeEnvFile = "$env:TEMP\compute-env.json"
    Write-JsonFile -Path $computeEnvFile -Object $computeEnv
    
    try {
        aws batch create-compute-environment --cli-input-json "file://$computeEnvFile" --region $region
        Write-Host "   SUCCESS: Compute environment created" -ForegroundColor Green
        Write-Host "   Waiting for environment to become VALID..." -ForegroundColor Yellow
        Start-Sleep -Seconds 10
    } catch {
        Write-Host "   ERROR: Failed to create compute environment: $_" -ForegroundColor Red
        exit 1
    }
}

Write-Host ""

# ============================================
# 2. Check/Create Job Queue
# ============================================
Write-Host "2. Checking Job Queue..." -ForegroundColor Cyan

$ErrorActionPreference = "SilentlyContinue"
$null = aws batch describe-job-queues --job-queues live2-job-queue --region $region 2>&1
$queueExists = $LASTEXITCODE -eq 0
$ErrorActionPreference = "Stop"

if ($queueExists) {
    Write-Host "   Job queue exists" -ForegroundColor Green
} else {
    Write-Host "   Creating Job Queue..." -ForegroundColor Yellow
    
    $computeEnvOrder = @(
        @{
            order = 1
            computeEnvironment = "live2-compute-env"
        }
    )
    
    $queueOrderFile = "$env:TEMP\queue-order.json"
    Write-JsonFile -Path $queueOrderFile -Object $computeEnvOrder
    
    try {
        aws batch create-job-queue `
            --job-queue-name live2-job-queue `
            --priority 1 `
            --compute-environment-order "file://$queueOrderFile" `
            --state ENABLED `
            --region $region
        
        Write-Host "   SUCCESS: Job queue created" -ForegroundColor Green
    } catch {
        Write-Host "   ERROR: Failed to create job queue: $_" -ForegroundColor Red
        exit 1
    }
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host "Batch Queue Setup Complete!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Compute Environment: live2-compute-env" -ForegroundColor Cyan
Write-Host "Job Queue: live2-job-queue" -ForegroundColor Cyan
Write-Host ""

