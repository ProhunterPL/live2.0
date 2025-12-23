---
date: 2025-12-23
label: [guide, cost-optimization]
---

# Cost Optimization - Split Deploy

**Jak minimalizować koszty w architekturze Split Deploy**

---

## 💰 Model Kosztowy

### DigitalOcean (Always-On)
- **Droplet:** $24/mo (2 vCPU, 4GB RAM)
- **Koszt stały:** $24/mo (płacisz zawsze, nawet bez klientów)

### AWS (On-Demand Compute)
- **AWS Batch:** Scale to ZERO (0 vCPUs gdy brak jobów)
- **Koszt:** $0 gdy brak jobów ✅
- **Koszt per job:** ~$0.10-0.50/job (zależnie od czasu trwania)

### Supabase
- **Free tier:** $0/mo (do 500MB DB, 2GB bandwidth)
- **Pro tier:** $25/mo (jeśli potrzebne)

### Redis
- **Redis Labs Free:** $0/mo (30MB)
- **Lub własny na DO:** wliczone w Droplet

### **Total Minimum:** ~$24/mo (tylko DO Droplet)

---

## 🎯 Strategia: Zero Cost When Idle

### AWS Batch - Automatyczne Scale to Zero

**Kluczowa funkcja:** AWS Batch automatycznie:
- ✅ Scale down do 0 vCPUs gdy brak jobów
- ✅ Scale up tylko gdy job jest submitted
- ✅ Scale down po zakończeniu joba

**Konfiguracja:**
```json
{
  "computeResources": {
    "minvCpus": 0,    // ← ZERO! Nie płacisz gdy brak jobów
    "maxvCpus": 32,   // Maksymalna pojemność
    "desiredvCpus": 0 // ZERO gdy idle
  }
}
```

**Koszt:**
- **0 jobów:** $0/mo ✅
- **1 job (30min):** ~$0.10
- **10 jobów/miesiąc:** ~$1-5/mo

---

## 📊 Przykładowe Scenariusze Kosztowe

### Scenariusz 1: Brak Klientów (MVP Start)
- **DO Droplet:** $24/mo
- **AWS Batch:** $0 (0 jobów)
- **Supabase:** $0 (free tier)
- **Redis:** $0 (free tier)
- **Total:** $24/mo

### Scenariusz 2: 1 Klient, 5 Jobów/Miesiąc
- **DO Droplet:** $24/mo
- **AWS Batch:** ~$0.50 (5 jobów × $0.10)
- **Supabase:** $0 (free tier)
- **Redis:** $0 (free tier)
- **Total:** ~$24.50/mo

### Scenariusz 3: 10 Klientów, 50 Jobów/Miesiąc
- **DO Droplet:** $24/mo
- **AWS Batch:** ~$5 (50 jobów × $0.10)
- **Supabase:** $0 (free tier) lub $25 (pro)
- **Redis:** $0 (free tier)
- **Total:** ~$29-54/mo

### Scenariusz 4: 100 Klientów, 500 Jobów/Miesiąc
- **DO Droplet:** $48/mo (upgrade do 4 vCPU)
- **AWS Batch:** ~$50 (500 jobów)
- **Supabase:** $25 (pro tier)
- **Redis:** $0 (free tier)
- **Total:** ~$123/mo

---

## ✅ Co Można Zrobić BEZ Kosztów AWS

### Przygotowanie Infrastruktury (Zero Cost)

Możesz przygotować całą infrastrukturę AWS **bez uruchamiania jobów**:

1. ✅ **S3 Bucket** - $0 (płacisz tylko za storage, nie za bucket)
2. ✅ **ECR Repository** - $0 (płacisz tylko za storage obrazów)
3. ✅ **IAM Roles/Users** - $0 (zawsze darmowe)
4. ✅ **Batch Compute Environment** - $0 (gdy minvCpus=0, nie ma kosztów)
5. ✅ **Batch Job Queue** - $0 (zawsze darmowe)
6. ✅ **Batch Job Definition** - $0 (zawsze darmowe)

**Koszt przygotowania:** $0 ✅

### Koszty Pojawiają Się Tylko Gdy:
- 📦 **S3 Storage:** ~$0.023/GB/mo (płacisz za przechowywane artefakty)
- 🐳 **ECR Storage:** ~$0.10/GB/mo (płacisz za przechowywane obrazy Docker)
- ⚡ **Batch Compute:** ~$0.10-0.50/job (płacisz tylko gdy job się wykonuje)

---

## 🛡️ Guardrails Kosztowe

### 1. AWS Budgets
```bash
# Utwórz budżet $10/mo dla AWS
aws budgets create-budget \
  --account-id YOUR_ACCOUNT \
  --budget file://budget.json
```

**Alerty:**
- 50% → Email
- 80% → Email + Slack
- 100% → Email + Slack + Auto-stop (opcjonalnie)

### 2. Batch Compute Limits
```json
{
  "maxvCpus": 32,  // Maksymalna pojemność
  "minvCpus": 0    // Zawsze zero gdy idle
}
```

### 3. S3 Lifecycle Policy
```json
{
  "Rules": [{
    "Expiration": {"Days": 90}  // Automatyczne usuwanie starych artefaktów
  }]
}
```

### 4. Rate Limiting (DO)
- Ogranicza liczbę jobów per user
- Zapobiega nadużyciom
- Implementowane w Redis (darmowe)

---

## 📋 Checklist: Przygotowanie Bez Kosztów

Możesz wykonać wszystkie te kroki **BEZ PŁACENIA**:

- [ ] Utwórz S3 bucket (zero cost)
- [ ] Utwórz ECR repository (zero cost)
- [ ] Utwórz IAM roles/users (zero cost)
- [ ] Utwórz Batch compute environment (minvCpus=0 → zero cost)
- [ ] Utwórz Batch job queue (zero cost)
- [ ] Utwórz Batch job definition (zero cost)
- [ ] Build & push Docker image (mały koszt storage: ~$0.10/mo)
- [ ] Skonfiguruj DO backend (już masz Droplet)

**Total:** ~$0.10/mo (tylko storage obrazu Docker)

---

## 🚀 Kiedy Pojawiają Się Koszty

### Koszty Pojawiają Się Tylko Gdy:

1. **Klient subskrybuje** → Stripe transaction fee (2.9% + $0.30)
2. **Klient uruchamia job** → AWS Batch compute (~$0.10-0.50/job)
3. **Job generuje artefakty** → S3 storage (~$0.023/GB/mo)

### Przykład: Pierwszy Płatny Klient

**Scenariusz:**
- Klient subskrybuje Hobby tier ($10/mo)
- Uruchamia 3 joby w pierwszym miesiącu

**Koszty:**
- DO Droplet: $24/mo (już płacisz)
- AWS Batch: ~$0.30 (3 joby)
- S3 Storage: ~$0.05 (artefakty)
- Stripe fee: $0.59 (2.9% + $0.30)

**Total dodatkowe koszty:** ~$0.94

**Przychód:** $10/mo
**Zysk:** $9.06/mo (po pierwszym kliencie)

---

## 💡 Best Practices

### 1. Zawsze Używaj Spot Instances
```json
{
  "type": "SPOT",           // 70% oszczędności
  "bidPercentage": 100      // Max savings
}
```

### 2. Auto-Terminate Jobs
```json
{
  "timeout": {
    "attemptDurationSeconds": 3600  // Max 1h per job
  }
}
```

### 3. Cleanup Stare Artefakty
```json
{
  "lifecycle": {
    "expiration": {"Days": 90}  // Auto-delete po 90 dniach
  }
}
```

### 4. Monitoruj Koszty
```bash
# Sprawdź koszty AWS
aws ce get-cost-and-usage \
  --time-period Start=2025-12-01,End=2025-12-31 \
  --granularity MONTHLY \
  --metrics BlendedCost
```

---

## ✅ Podsumowanie

### ✅ Możesz Przygotować Wszystko Bez Kosztów
- S3, ECR, IAM, Batch (gdy minvCpus=0) → $0
- Tylko storage obrazu Docker → ~$0.10/mo

### ✅ Koszty Pojawiają Się Tylko Gdy:
- Klient uruchamia job → ~$0.10-0.50/job
- Artefakty są przechowywane → ~$0.023/GB/mo

### ✅ AWS Batch Scale to Zero
- Automatycznie: 0 vCPUs gdy brak jobów
- Automatycznie: Scale up gdy job submitted
- Automatycznie: Scale down po zakończeniu

**Total minimum:** ~$24/mo (tylko DO Droplet) ✅

---

**Ostatnia aktualizacja:** 2025-12-23

