---
date: 2025-12-23
label: [guide]
---
# Rekomendacja Instancji AWS - Live 2.0 z Monetyzacją

## 🎯 Na Początek (MVP / Staging)

### **Opcja 1: t3.medium (ZALECANA dla startu)**

**Specyfikacja:**
- **vCPUs:** 2
- **RAM:** 4 GB
- **Network:** Up to 5 Gbps
- **Koszt:** ~$0.04/godz = ~$1/dzień = ~$30/miesiąc

**Dlaczego:**
- ✅ Wystarczająca dla API + Frontend + Billing
- ✅ Niski koszt (idealne do testów)
- ✅ Można skalować w górę gdy potrzeba
- ✅ Burstable performance (T3) - OK dla API

**Dla czego wystarczy:**
- FastAPI backend (API v1 + billing)
- React frontend (statyczne pliki)
- PostgreSQL (external - Supabase)
- Redis (external - Redis Labs)
- Stripe webhooks

**Limitacje:**
- ❌ Nie uruchomisz symulacji (za mało CPU/RAM)
- ⚠️ Burstable credits - może być wolniejsze przy ciągłym obciążeniu

---

### **Opcja 2: t3.large (Jeśli oczekujesz większego ruchu)**

**Specyfikacja:**
- **vCPUs:** 2
- **RAM:** 8 GB
- **Network:** Up to 5 Gbps
- **Koszt:** ~$0.08/godz = ~$2/dzień = ~$60/miesiąc

**Dlaczego:**
- ✅ Więcej RAM (lepsze dla wielu równoczesnych requestów)
- ✅ Większy baseline performance
- ✅ Nadal tanie

---

### **Opcja 3: t3.small (Minimalne - tylko testy)**

**Specyfikacja:**
- **vCPUs:** 2
- **RAM:** 2 GB
- **Network:** Up to 5 Gbps
- **Koszt:** ~$0.02/godz = ~$0.50/dzień = ~$15/miesiąc

**Dlaczego:**
- ✅ Najtańsze
- ✅ OK dla development/staging

**Limitacje:**
- ⚠️ Tylko 2 GB RAM - może być ciasno
- ⚠️ Nie dla production z realnym ruchem

---

## 🚀 Production (Gdy masz użytkowników)

### **Opcja 1: t3.xlarge (ZALECANA dla production)**

**Specyfikacja:**
- **vCPUs:** 4
- **RAM:** 16 GB
- **Network:** Up to 5 Gbps
- **Koszt:** ~$0.17/godz = ~$4/dzień = ~$120/miesiąc

**Dlaczego:**
- ✅ Wystarczająca dla production API
- ✅ Może obsłużyć setki requestów/min
- ✅ Więcej RAM dla cache/sessions
- ✅ Nadal reasonable cost

---

### **Opcja 2: m5.large (Jeśli potrzebujesz więcej CPU)**

**Specyfikacja:**
- **vCPUs:** 2
- **RAM:** 8 GB
- **Network:** Up to 10 Gbps
- **Koszt:** ~$0.10/godz = ~$2.40/dzień = ~$72/miesiąc

**Dlaczego:**
- ✅ General purpose (nie burstable)
- ✅ Stała wydajność
- ✅ Lepsze dla ciągłego obciążenia

---

## 🔬 Jeśli Planujesz Uruchamiać Symulacje

### **c6i.16xlarge (Compute-Optimized)**

**Specyfikacja:**
- **vCPUs:** 64
- **RAM:** 128 GB
- **Network:** Up to 50 Gbps
- **Koszt:** ~$2.72/godz = ~$65/dzień = ~$1950/miesiąc

**Dlaczego:**
- ✅ Idealne dla symulacji (Phase 2B)
- ✅ 16-32 równoległych symulacji
- ✅ Najszybsze dla compute-heavy workloads

**Uwaga:**
- 💰 Drogie - używaj tylko gdy uruchamiasz symulacje
- ⏱️ Można użyć Spot Instances (oszczędność 50-90%)

---

## 📊 Porównanie

| Instancja | vCPU | RAM | Koszt/godz | Koszt/mies | Użycie |
|-----------|------|-----|------------|------------|--------|
| **t3.small** | 2 | 2 GB | $0.02 | $15 | Dev/Test |
| **t3.medium** | 2 | 4 GB | $0.04 | $30 | **MVP/Staging** ⭐ |
| **t3.large** | 2 | 8 GB | $0.08 | $60 | Większy ruch |
| **t3.xlarge** | 4 | 16 GB | $0.17 | $120 | **Production** ⭐ |
| **m5.large** | 2 | 8 GB | $0.10 | $72 | Stała wydajność |
| **c6i.16xlarge** | 64 | 128 GB | $2.72 | $1950 | Symulacje |

---

## 🎯 Moja Rekomendacja

### **Na Początek: t3.medium**

**Dlaczego:**
1. ✅ Wystarczająca dla MVP (API + Frontend + Billing)
2. ✅ Niski koszt (~$30/miesiąc)
3. ✅ Można łatwo skalować w górę (resize instance)
4. ✅ External services (Supabase DB, Redis Labs) - nie obciążają instancji

**Setup:**
```bash
# AWS Console → EC2 → Launch Instance
Instance type: t3.medium
OS: Ubuntu 22.04 LTS
Storage: 20 GB (wystarczy)
Security Group: Porty 22 (SSH), 80 (HTTP), 443 (HTTPS), 8000 (API - opcjonalnie)
```

**Kiedy skalować:**
- → **t3.large** gdy: >100 requestów/min, RAM >80%
- → **t3.xlarge** gdy: >500 requestów/min, production traffic
- → **c6i.16xlarge** gdy: chcesz uruchamiać symulacje

---

## 💡 Tips

### 1. **Start Small, Scale Up**
- Zaczynaj od t3.medium
- Monitoruj CloudWatch metrics
- Skaluj w górę gdy potrzeba (resize instance)

### 2. **Use Spot Instances dla Symulacji**
- Jeśli uruchamiasz symulacje → użyj Spot Instances
- Oszczędność 50-90%
- OK dla batch jobs (symulacje)

### 3. **External Services**
- PostgreSQL → Supabase (free tier OK dla startu)
- Redis → Redis Labs (free tier OK dla startu)
- Nie obciążają instancji EC2

### 4. **Auto Scaling (Future)**
- Gdy masz production traffic → setup Auto Scaling Group
- Automatycznie dodaje instancje przy obciążeniu

---

## 📝 Quick Start Command

```bash
# Launch t3.medium instance
aws ec2 run-instances \
  --image-id ami-0c55b159cbfafe1f0 \
  --instance-type t3.medium \
  --key-name your-key \
  --security-group-ids sg-xxxxx \
  --subnet-id subnet-xxxxx \
  --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":20}}]' \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=live2-mvp}]'
```

---

## ✅ Checklist

- [ ] Wybrano instancję (t3.medium recommended)
- [ ] Security Group skonfigurowany (22, 80, 443)
- [ ] Key pair utworzony
- [ ] External services setup (Supabase, Redis Labs)
- [ ] Environment variables przygotowane
- [ ] Domain name + SSL certificate (Let's Encrypt)

---

**Ostatnia aktualizacja:** 2025-12-23

