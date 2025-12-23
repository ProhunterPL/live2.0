# Install doctl and AWS CLI on Windows
# Run as Administrator: powershell -ExecutionPolicy Bypass -File install_tools_windows.ps1

Write-Host "Installing deployment tools..." -ForegroundColor Green

# Check if running as Administrator
$isAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
if (-not $isAdmin) {
    Write-Host "WARNING: This script should be run as Administrator" -ForegroundColor Yellow
    Write-Host "Right-click PowerShell -> Run as Administrator" -ForegroundColor Yellow
    exit 1
}

# 1. Install doctl (DigitalOcean CLI)
Write-Host ""
Write-Host "Installing doctl (DigitalOcean CLI)..." -ForegroundColor Cyan

# Check if doctl is already installed
$doctlInstalled = Get-Command doctl -ErrorAction SilentlyContinue
if ($doctlInstalled) {
    Write-Host "SUCCESS: doctl is already installed: $(doctl version)" -ForegroundColor Green
} else {
    Write-Host "Installing doctl from GitHub releases..." -ForegroundColor Yellow
    
    # Download doctl from GitHub releases
    $doctlVersion = "1.104.0"
    $doctlUrl = "https://github.com/digitalocean/doctl/releases/download/v$doctlVersion/doctl-$doctlVersion-windows-amd64.zip"
    $doctlZip = "$env:TEMP\doctl.zip"
    $doctlDir = "$env:ProgramFiles\doctl"
    
    Write-Host "Downloading doctl v$doctlVersion..." -ForegroundColor Yellow
    try {
        Invoke-WebRequest -Uri $doctlUrl -OutFile $doctlZip -UseBasicParsing
        
        # Create directory
        if (-not (Test-Path $doctlDir)) {
            New-Item -ItemType Directory -Path $doctlDir -Force | Out-Null
        }
        
        # Extract
        Write-Host "Extracting doctl..." -ForegroundColor Yellow
        Expand-Archive -Path $doctlZip -DestinationPath $doctlDir -Force
        
        # Add to PATH (system-wide)
        $currentPath = [System.Environment]::GetEnvironmentVariable("Path","Machine")
        if ($currentPath -notlike "*$doctlDir*") {
            [System.Environment]::SetEnvironmentVariable("Path", "$currentPath;$doctlDir", "Machine")
        }
        
        # Refresh PATH for current session
        $env:Path = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
        
        # Cleanup
        Remove-Item $doctlZip -ErrorAction SilentlyContinue
        
        # Verify
        Start-Sleep -Seconds 1
        $doctlInstalled = Get-Command doctl -ErrorAction SilentlyContinue
        if ($doctlInstalled) {
            Write-Host "SUCCESS: doctl installed successfully: $(doctl version)" -ForegroundColor Green
        } else {
            Write-Host "WARNING: doctl installed but not in PATH. Please restart PowerShell." -ForegroundColor Yellow
        }
    } catch {
        Write-Host "ERROR: Failed to install doctl: $_" -ForegroundColor Red
        Write-Host "You can install manually from: https://github.com/digitalocean/doctl/releases" -ForegroundColor Yellow
    }
}

# 2. Install AWS CLI
Write-Host ""
Write-Host "Installing AWS CLI..." -ForegroundColor Cyan

$awsInstalled = Get-Command aws -ErrorAction SilentlyContinue
if ($awsInstalled) {
    Write-Host "SUCCESS: AWS CLI is already installed: $(aws --version)" -ForegroundColor Green
} else {
    Write-Host "Installing AWS CLI..." -ForegroundColor Yellow
    
    # Download AWS CLI MSI installer
    $awsCliUrl = "https://awscli.amazonaws.com/AWSCLIV2.msi"
    $awsCliInstaller = "$env:TEMP\AWSCLIV2.msi"
    
    Write-Host "Downloading AWS CLI installer..." -ForegroundColor Yellow
    Invoke-WebRequest -Uri $awsCliUrl -OutFile $awsCliInstaller
    
    Write-Host "Installing AWS CLI (this may take a few minutes)..." -ForegroundColor Yellow
    Start-Process msiexec.exe -Wait -ArgumentList "/i $awsCliInstaller /quiet /norestart"
    
    # Cleanup
    Remove-Item $awsCliInstaller -ErrorAction SilentlyContinue
    
    # Refresh PATH
    $machinePath = [System.Environment]::GetEnvironmentVariable("Path","Machine")
    $userPath = [System.Environment]::GetEnvironmentVariable("Path","User")
    $env:Path = "$machinePath;$userPath"
    
    # Verify installation
    Start-Sleep -Seconds 2
    $awsInstalled = Get-Command aws -ErrorAction SilentlyContinue
    if ($awsInstalled) {
        Write-Host "SUCCESS: AWS CLI installed successfully: $(aws --version)" -ForegroundColor Green
    } else {
        Write-Host "WARNING: AWS CLI installed but not in PATH. Please restart PowerShell." -ForegroundColor Yellow
    }
}

# 3. Verify installations
Write-Host ""
Write-Host "Installation Summary:" -ForegroundColor Green
Write-Host "====================" -ForegroundColor Green

if (Get-Command doctl -ErrorAction SilentlyContinue) {
    Write-Host "SUCCESS: doctl: $(doctl version)" -ForegroundColor Green
} else {
    Write-Host "ERROR: doctl: Not found (restart PowerShell)" -ForegroundColor Red
}

if (Get-Command aws -ErrorAction SilentlyContinue) {
    Write-Host "SUCCESS: AWS CLI: $(aws --version)" -ForegroundColor Green
} else {
    Write-Host "ERROR: AWS CLI: Not found (restart PowerShell)" -ForegroundColor Red
}

Write-Host ""
Write-Host "Next Steps:" -ForegroundColor Cyan
Write-Host "1. Restart PowerShell if tools are not found" -ForegroundColor Yellow
Write-Host "2. Configure doctl: doctl auth init" -ForegroundColor Yellow
Write-Host "3. Configure AWS CLI: aws configure" -ForegroundColor Yellow
Write-Host ""
Write-Host "Remember: AWS Batch scales to ZERO - no cost when no jobs running!" -ForegroundColor Green
