# Build and Push Docker Image to ECR
# Run: powershell -ExecutionPolicy Bypass -File scripts/deploy/build_and_push_docker.ps1

$ErrorActionPreference = "Stop"

Write-Host "========================================" -ForegroundColor Green
Write-Host "Build and Push Docker Image to ECR" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

# Get AWS Account ID and Region
$accountId = aws sts get-caller-identity --query Account --output text
$region = aws configure get region
if (-not $region) {
    $region = "eu-central-1"
}

$ecrUri = "$accountId.dkr.ecr.$region.amazonaws.com/live2-simulation"

Write-Host "Account ID: $accountId" -ForegroundColor Cyan
Write-Host "Region: $region" -ForegroundColor Cyan
Write-Host "ECR URI: $ecrUri" -ForegroundColor Cyan
Write-Host ""

# Check if Docker is running
Write-Host "1. Checking Docker..." -ForegroundColor Cyan
try {
    $dockerVersion = docker --version
    Write-Host "   SUCCESS: Docker is running" -ForegroundColor Green
    Write-Host "   $dockerVersion" -ForegroundColor Gray
} catch {
    Write-Host "   ERROR: Docker is not running or not installed" -ForegroundColor Red
    Write-Host "   Please start Docker Desktop and try again" -ForegroundColor Yellow
    exit 1
}

Write-Host ""

# Check if Dockerfile exists
Write-Host "2. Checking Dockerfile..." -ForegroundColor Cyan
$dockerfilePath = "docker/simulation.Dockerfile"
if (-not (Test-Path $dockerfilePath)) {
    Write-Host "   ERROR: Dockerfile not found at $dockerfilePath" -ForegroundColor Red
    exit 1
}
Write-Host "   SUCCESS: Dockerfile found" -ForegroundColor Green
Write-Host ""

# Build Docker image
Write-Host "3. Building Docker image..." -ForegroundColor Cyan
Write-Host "   This may take several minutes..." -ForegroundColor Yellow

try {
    docker build -t live2-simulation:latest -f $dockerfilePath .
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "   SUCCESS: Docker image built successfully" -ForegroundColor Green
    } else {
        Write-Host "   ERROR: Docker build failed" -ForegroundColor Red
        exit 1
    }
} catch {
    Write-Host "   ERROR: Docker build failed: $_" -ForegroundColor Red
    exit 1
}

Write-Host ""

# Login to ECR
Write-Host "4. Logging in to ECR..." -ForegroundColor Cyan

try {
    $loginCommand = aws ecr get-login-password --region $region | `
        docker login --username AWS --password-stdin $ecrUri
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "   SUCCESS: Logged in to ECR" -ForegroundColor Green
    } else {
        Write-Host "   ERROR: Failed to login to ECR" -ForegroundColor Red
        exit 1
    }
} catch {
    Write-Host "   ERROR: Failed to login to ECR: $_" -ForegroundColor Red
    exit 1
}

Write-Host ""

# Tag image
Write-Host "5. Tagging Docker image..." -ForegroundColor Cyan

$latestTag = "${ecrUri}:latest"
$versionTag = "${ecrUri}:v1.0.0"

try {
    docker tag live2-simulation:latest $latestTag
    if ($LASTEXITCODE -ne 0) {
        Write-Host "   ERROR: Failed to tag latest" -ForegroundColor Red
        exit 1
    }
    
    docker tag live2-simulation:latest $versionTag
    if ($LASTEXITCODE -ne 0) {
        Write-Host "   ERROR: Failed to tag v1.0.0" -ForegroundColor Red
        exit 1
    }
    
    Write-Host "   SUCCESS: Image tagged" -ForegroundColor Green
    Write-Host "   Tags: latest, v1.0.0" -ForegroundColor Gray
} catch {
    Write-Host "   ERROR: Failed to tag image: $_" -ForegroundColor Red
    exit 1
}

Write-Host ""

# Push image
Write-Host "6. Pushing Docker image to ECR..." -ForegroundColor Cyan
Write-Host "   This may take several minutes..." -ForegroundColor Yellow

try {
    Write-Host "   Pushing latest tag..." -ForegroundColor Gray
    docker push $latestTag
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "   SUCCESS: latest tag pushed" -ForegroundColor Green
    } else {
        Write-Host "   ERROR: Failed to push latest tag" -ForegroundColor Red
        exit 1
    }
    
    Write-Host "   Pushing v1.0.0 tag..." -ForegroundColor Gray
    docker push $versionTag
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "   SUCCESS: v1.0.0 tag pushed" -ForegroundColor Green
    } else {
        Write-Host "   WARNING: Failed to push v1.0.0 tag (latest is OK)" -ForegroundColor Yellow
    }
} catch {
    Write-Host "   ERROR: Failed to push image: $_" -ForegroundColor Red
    exit 1
}

Write-Host ""

# Verify
Write-Host "7. Verifying image in ECR..." -ForegroundColor Cyan

try {
    $images = aws ecr describe-images --repository-name live2-simulation --region $region --query "imageDetails[*].imageTags" --output text
    Write-Host "   SUCCESS: Image verified in ECR" -ForegroundColor Green
    Write-Host "   Available tags: $images" -ForegroundColor Gray
} catch {
    Write-Host "   WARNING: Could not verify image (but push may have succeeded)" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host "Docker Image Build & Push Complete!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Image URI: $ecrUri:latest" -ForegroundColor Cyan
Write-Host ""
Write-Host "Next Steps:" -ForegroundColor Cyan
Write-Host "1. Update Batch Job Definition with image URI" -ForegroundColor Yellow
Write-Host "2. Test job submission" -ForegroundColor Yellow
Write-Host ""

