# Skrypt do skopiowania plikow na AWS
# Uzycie: .\copy_to_aws.ps1

# EDYTUJ TO:
# Use environment variables or command-line arguments for security
param(
    [Parameter(Mandatory=$false)]
    [string]$AwsIp = $env:AWS_IP,
    [Parameter(Mandatory=$false)]
    [string]$KeyPath = $env:AWS_SSH_KEY_PATH,
    [Parameter(Mandatory=$false)]
    [bool]$UseKey = $true
)

# Validate required parameters
if (-not $AwsIp) {
    Write-Host "ERROR: AWS IP not provided. Set AWS_IP environment variable or use -AwsIp parameter" -ForegroundColor Red
    exit 1
}

if ($UseKey -and -not $KeyPath) {
    Write-Host "ERROR: SSH key path not provided. Set AWS_SSH_KEY_PATH environment variable or use -KeyPath parameter" -ForegroundColor Red
    exit 1
}

$AWS_IP = $AwsIp
$KEY_PATH = $KeyPath
$USE_KEY = $UseKey

Write-Host "Kopiowanie plikow na AWS..." -ForegroundColor Green

if ($USE_KEY) {
    # Z kluczem SSH
    scp -i $KEY_PATH diagnose_round1.sh ubuntu@${AWS_IP}:~/live2.0/
    scp -i $KEY_PATH AWS_EMERGENCY_FIX.txt ubuntu@${AWS_IP}:~/live2.0/
    scp -i $KEY_PATH AWS_RECOMMENDED_ACTION.txt ubuntu@${AWS_IP}:~/live2.0/
    scp -i $KEY_PATH aws_start_missing_9.sh ubuntu@${AWS_IP}:~/live2.0/
} else {
    # Bez klucza (z haslem)
    scp diagnose_round1.sh ubuntu@${AWS_IP}:~/live2.0/
    scp AWS_EMERGENCY_FIX.txt ubuntu@${AWS_IP}:~/live2.0/
    scp AWS_RECOMMENDED_ACTION.txt ubuntu@${AWS_IP}:~/live2.0/
    scp aws_start_missing_9.sh ubuntu@${AWS_IP}:~/live2.0/
}

Write-Host "Gotowe! Teraz polacz sie przez SSH i uruchom:" -ForegroundColor Green
Write-Host "ssh ubuntu@$AWS_IP" -ForegroundColor Yellow
Write-Host "cd live2.0" -ForegroundColor Yellow
Write-Host "bash diagnose_round1.sh" -ForegroundColor Yellow

