---
date: 2025-12-23
label: [guide, installation]
---

# Instalacja Narzędzi - doctl i AWS CLI

**Przewodnik instalacji narzędzi potrzebnych do wdrożenia**

---

## 🎯 Co Instalujemy

- **doctl** - DigitalOcean CLI (do zarządzania Dropletami)
- **AWS CLI** - AWS Command Line Interface (do zarządzania AWS)

---

## 🪟 Windows

### Opcja 1: Automatyczna Instalacja (Rekomendowane)

```powershell
# Uruchom PowerShell jako Administrator
# (Right-click → Run as Administrator)

# Przejdź do katalogu projektu
cd C:\Users\klawi\live2.0\live2.0

# Uruchom skrypt instalacyjny
powershell -ExecutionPolicy Bypass -File scripts/deploy/install_tools_windows.ps1
```

### Opcja 2: Ręczna Instalacja

#### doctl (DigitalOcean CLI)

**Przez Chocolatey:**
```powershell
# Jeśli nie masz Chocolatey, zainstaluj:
Set-ExecutionPolicy Bypass -Scope Process -Force
[System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072
iex ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))

# Zainstaluj doctl
choco install doctl -y
```

**Lub pobierz ręcznie:**
1. Idź na https://github.com/digitalocean/doctl/releases
2. Pobierz `doctl-X.X.X-windows-amd64.zip`
3. Rozpakuj i dodaj do PATH

#### AWS CLI

**Przez MSI Installer:**
1. Pobierz: https://awscli.amazonaws.com/AWSCLIV2.msi
2. Uruchom installer
3. Postępuj zgodnie z instrukcjami

**Lub przez PowerShell:**
```powershell
# Download
$awsCliUrl = "https://awscli.amazonaws.com/AWSCLIV2.msi"
$awsCliInstaller = "$env:TEMP\AWSCLIV2.msi"
Invoke-WebRequest -Uri $awsCliUrl -OutFile $awsCliInstaller

# Install
Start-Process msiexec.exe -Wait -ArgumentList "/i $awsCliInstaller /quiet /norestart"

# Cleanup
Remove-Item $awsCliInstaller
```

---

## 🍎 macOS

```bash
# doctl
brew install doctl

# AWS CLI
brew install awscli
```

---

## 🐧 Linux

```bash
# doctl
cd ~
wget https://github.com/digitalocean/doctl/releases/download/v1.104.0/doctl-1.104.0-linux-amd64.tar.gz
tar xf doctl-1.104.0-linux-amd64.tar.gz
sudo mv doctl /usr/local/bin

# AWS CLI
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

---

## ✅ Weryfikacja Instalacji

```powershell
# Sprawdź doctl
doctl version

# Sprawdź AWS CLI
aws --version
```

---

## 🔐 Konfiguracja

### doctl (DigitalOcean)

```bash
# Login do DigitalOcean
doctl auth init

# Wprowadź token z DigitalOcean Dashboard
# (Settings → API → Generate New Token)
```

### AWS CLI

```bash
# Konfiguracja
aws configure

# Wprowadź:
# - AWS Access Key ID
# - AWS Secret Access Key
# - Default region (np. us-east-1)
# - Default output format (json)
```

**Gdzie znaleźć AWS Credentials:**
1. AWS Console → IAM → Users
2. Wybierz użytkownika (lub utwórz nowego)
3. Security credentials → Create access key

---

## 💡 Ważne: AWS Batch Scale to Zero

**Kluczowa informacja:** AWS Batch automatycznie scale to zero!

- ✅ **minvCpus: 0** → Nie płacisz gdy brak jobów
- ✅ **Koszt pojawia się tylko gdy job się wykonuje**
- ✅ **Po zakończeniu joba → automatycznie scale down do 0**

**Możesz przygotować całą infrastrukturę AWS BEZ KOSZTÓW:**
- S3 bucket (płacisz tylko za storage, nie za bucket)
- ECR repository (płacisz tylko za storage obrazów)
- IAM roles/users (zawsze darmowe)
- Batch compute environment (gdy minvCpus=0 → zero cost)
- Batch job queue (zawsze darmowe)
- Batch job definition (zawsze darmowe)

**Koszt przygotowania:** ~$0 ✅

**Koszt pojawia się tylko gdy:**
- Klient uruchamia job → ~$0.10-0.50/job
- Artefakty są przechowywane → ~$0.023/GB/mo

---

## 📚 Następne Kroki

Po zainstalowaniu narzędzi:

1. **Skonfiguruj doctl:** `doctl auth init`
2. **Skonfiguruj AWS CLI:** `aws configure`
3. **Przeczytaj:** [`COST_OPTIMIZATION.md`](COST_OPTIMIZATION.md) - szczegóły o kosztach
4. **Zacznij wdrożenie:** [`QUICK_START_DEPLOY.md`](QUICK_START_DEPLOY.md)

---

**Ostatnia aktualizacja:** 2025-12-23

