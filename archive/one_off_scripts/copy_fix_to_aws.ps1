# Deploy Phase 2B Cluster Fix to AWS
# ===================================

$ErrorActionPreference = "Stop"

Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host "🚀 DEPLOYING CLUSTER FIX TO AWS" -ForegroundColor Cyan
Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host ""

# Get AWS IP
$awsIp = "ip-172-31-0-42"
Write-Host "Target AWS instance: $awsIp" -ForegroundColor Yellow
Write-Host ""

# Confirm
$response = Read-Host "Deploy fix to AWS? (y/n)"
if ($response -ne "y") {
    Write-Host "Cancelled." -ForegroundColor Red
    exit
}

Write-Host ""
Write-Host "📦 Copying files to AWS..." -ForegroundColor Green
Write-Host ""

# Copy all fix files
$files = @(
    "aws_test/DEPLOY_FIX_NOW.sh",
    "aws_test/EMERGENCY_SUMMARY.md",
    "aws_test/CLUSTER_FIX_INSTRUCTIONS.md",
    "aws_test/scripts/kill_stuck_simulations.sh",
    "aws_test/scripts/apply_cluster_fix.sh",
    "aws_test/scripts/restart_phase2b_safe.sh",
    "aws_test/scripts/check_actual_progress.py",
    "aws_test/scripts/monitor_by_filesize.py",
    "aws_test/configs/phase2_miller_urey_extended_SAFER.yaml",
    "aws_test/patches/disable_cluster_detection.patch"
)

foreach ($file in $files) {
    Write-Host "  → $file" -ForegroundColor Gray
    scp $file "ubuntu@${awsIp}:~/live2.0/$file"
    if ($LASTEXITCODE -ne 0) {
        Write-Host "    ❌ Failed to copy $file" -ForegroundColor Red
        exit 1
    }
}

Write-Host ""
Write-Host "================================================================================" -ForegroundColor Green
Write-Host "✅ ALL FILES DEPLOYED SUCCESSFULLY" -ForegroundColor Green
Write-Host "================================================================================" -ForegroundColor Green
Write-Host ""
Write-Host "Next steps on AWS:" -ForegroundColor Yellow
Write-Host ""
Write-Host "  ssh ubuntu@$awsIp" -ForegroundColor Cyan
Write-Host "  cd ~/live2.0" -ForegroundColor Cyan
Write-Host "  bash aws_test/DEPLOY_FIX_NOW.sh" -ForegroundColor Cyan
Write-Host ""
Write-Host "Or read instructions first:" -ForegroundColor Yellow
Write-Host "  cat aws_test/EMERGENCY_SUMMARY.md" -ForegroundColor Cyan
Write-Host ""

