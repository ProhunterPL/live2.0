# Create AWS Batch Job Definition
# Run: powershell -ExecutionPolicy Bypass -File scripts/deploy/create_job_definition.ps1

$ErrorActionPreference = "Stop"

Write-Host "========================================" -ForegroundColor Green
Write-Host "Creating AWS Batch Job Definition" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

# Get AWS Account ID and Region
$accountId = aws sts get-caller-identity --query Account --output text
$region = aws configure get region
if (-not $region) {
    $region = "eu-central-1"
}

$ecrUri = "$accountId.dkr.ecr.$region.amazonaws.com/live2-simulation"
$jobRoleArn = "arn:aws:iam::${accountId}:role/Live2JobRole"

Write-Host "Account ID: $accountId" -ForegroundColor Cyan
Write-Host "Region: $region" -ForegroundColor Cyan
Write-Host "ECR URI: $ecrUri:latest" -ForegroundColor Cyan
Write-Host "Job Role: $jobRoleArn" -ForegroundColor Cyan
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

# Create Job Definition
Write-Host "Creating Job Definition..." -ForegroundColor Cyan

$containerProps = @{
    image = "${ecrUri}:latest"
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

try {
    aws batch register-job-definition --cli-input-json "file://$jobDefFile"
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "   SUCCESS: Job definition created" -ForegroundColor Green
    } else {
        Write-Host "   ERROR: Failed to create job definition" -ForegroundColor Red
        exit 1
    }
} catch {
    Write-Host "   ERROR: Failed to create job definition: $_" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host "Job Definition Created!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Job Definition: live2-simulation" -ForegroundColor Cyan
Write-Host "Image: $ecrUri:latest" -ForegroundColor Cyan
Write-Host ""
Write-Host "Next Steps:" -ForegroundColor Cyan
Write-Host "1. Test job submission" -ForegroundColor Yellow
Write-Host "2. Monitor job execution in AWS Console" -ForegroundColor Yellow
Write-Host ""

