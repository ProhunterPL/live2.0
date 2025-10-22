# AWS / Cloud Deployment

Dokumentacja wdrożenia Live 2.0 na AWS i inne platformy cloud.

---

## 📄 Quick Start Guides

### [QUICK_AWS_COMMANDS.md](QUICK_AWS_COMMANDS.md)
Najczęściej używane komendy AWS:
- Setup instance
- Start/stop symulacji
- Monitoring
- Download results

### [QUICK_START_AWS_PIPELINE.md](QUICK_START_AWS_PIPELINE.md)
Pełny pipeline AWS:
- Krok po kroku setup
- Konfiguracja EC2
- Storage (S3)
- Networking

### [AWS_QUICK_START.md](AWS_QUICK_START.md)
Szybki start:
- Minimalna konfiguracja
- Przykładowe komendy
- Troubleshooting

---

## 🚀 Architektura AWS

### Komponenty:

```
┌─────────────────────────────────────────┐
│  EC2 Instance (Ubuntu 22.04)           │
│  ├─ Backend (FastAPI + Taichi)         │
│  ├─ Frontend (React + Vite)            │
│  └─ Nginx (reverse proxy)              │
└─────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────┐
│  S3 Bucket                              │
│  ├─ Snapshots                           │
│  ├─ Results                             │
│  └─ Logs                                │
└─────────────────────────────────────────┘
```

### Wymagania:
- **Instance type:** t3.xlarge lub większy (dla Taichi GPU)
- **Storage:** 100GB EBS
- **Memory:** 16GB RAM minimum
- **GPU:** Optional (CUDA support)

---

## 📋 Setup Checklist

### 1. Pre-requisites
```bash
- [ ] AWS Account
- [ ] SSH Key pair
- [ ] Security group configured
- [ ] S3 bucket created
```

### 2. Instance Setup
```bash
- [ ] Launch EC2 instance
- [ ] Install dependencies
- [ ] Clone repository
- [ ] Setup virtual environment
```

### 3. Application Setup
```bash
- [ ] Install Python packages
- [ ] Install Node.js packages
- [ ] Configure Nginx
- [ ] Setup SSL (optional)
```

### 4. Testing
```bash
- [ ] Test backend API
- [ ] Test frontend
- [ ] Test simulation run
- [ ] Verify S3 upload
```

---

## 🔧 Common Tasks

### Start Simulation
```bash
ssh -i key.pem ubuntu@ec2-xxx.compute.amazonaws.com
cd live2.0
./start_backend.ps1  # lub bash equivalent
```

### Monitor Progress
```bash
# Logs
tail -f logs/logs.txt

# Results
aws s3 ls s3://live2-results/

# Metrics
curl http://localhost:8000/simulation/current/metrics
```

### Download Results
```bash
aws s3 sync s3://live2-results/ ./results/
```

---

## 📊 Cost Estimation

### Typical Run (24h):
- EC2 t3.xlarge: ~$0.17/hour × 24 = $4.08
- Storage (100GB): ~$10/month
- S3 storage: ~$0.023/GB
- Data transfer: Variable

**Total: ~$5-10 per day**

---

## 🐛 Troubleshooting

### Instance Won't Start
- Check security group (ports 8000, 5173, 80, 443)
- Verify SSH key
- Check instance limits

### Backend Crashes
- Check memory (htop)
- Review logs (journalctl -u live2-backend)
- Verify Taichi installation

### Frontend Not Loading
- Check Nginx config
- Verify backend is running
- Check CORS settings

---

## 🔗 Zobacz Też

- [Cloud Deployment Guide](../../CLOUD_DEPLOYMENT_GUIDE.md) - Pełna dokumentacja
- [AWS Results Pipeline](../../AWS_RESULTS_PIPELINE.md) - Pipeline results
- [Session 2024-10-16](../../sessions/) - AWS integration session

