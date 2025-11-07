# Skrypt do skopiowania plikow na AWS
# Uzycie: .\copy_to_aws.ps1

# EDYTUJ TO:
$AWS_IP = "35.157.92.39"  # np. "3.15.123.45"
$KEY_PATH = "C:\Users\klawi\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem"  # Jesli masz klucz
$USE_KEY = $true  # Zmien na $true jesli uzywasz klucza

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

