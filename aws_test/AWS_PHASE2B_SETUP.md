# AWS Instance dla Phase 2B - Rekomendacje

**Data**: 24 października 2025  
**Cel**: Uruchomienie Phase 2B (30 dodatkowych symulacji + debug formamide)

---

## 🖥️ **Rekomendowana Instancja AWS**

### **Opcja 1: c6i.16xlarge (Zalecana)**
- **vCPUs**: 64
- **RAM**: 128 GB
- **Storage**: 100 GB SSD
- **Koszt**: ~$2.50/godzina
- **Uzasadnienie**: 
  - Wystarczająca moc dla 30 symulacji równolegle
  - 128 GB RAM dla długich symulacji (500K kroków)
  - Stabilna wydajność

### **Opcja 2: c6i.8xlarge (Ekonomiczna)**
- **vCPUs**: 32
- **RAM**: 64 GB
- **Storage**: 100 GB SSD
- **Koszt**: ~$1.25/godzina
- **Uzasadnienie**:
  - Wystarczająca dla 15-20 symulacji równolegle
  - Niższy koszt
  - Może wymagać uruchamiania w partiach

### **Opcja 3: c6i.4xlarge (Minimalna)**
- **vCPUs**: 16
- **RAM**: 32 GB
- **Storage**: 100 GB SSD
- **Koszt**: ~$0.63/godzina
- **Uzasadnienie**:
  - Najniższy koszt
  - Uruchamianie po 5-10 symulacji
  - Dłuższy czas wykonania

---

## 🚀 **Setup Instancji**

### **1. Uruchom Instancję**
```bash
# AWS Console lub CLI
aws ec2 run-instances \
  --image-id ami-0c02fb55956c7d316 \
  --instance-type c6i.16xlarge \
  --key-name your-key-name \
  --security-group-ids sg-xxxxxxxxx \
  --subnet-id subnet-xxxxxxxxx \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=live2-phase2b}]'
```

### **2. Połącz się z Instancją**
```bash
ssh -i your-key.pem ubuntu@<instance-ip>
```

### **3. Uruchom Setup**
```bash
# Pobierz i uruchom setup script
wget https://raw.githubusercontent.com/ProhunterPL/live2.0/main/setup_aws_instance.sh
bash setup_aws_instance.sh
```

---

## 📋 **Phase 2B na AWS - Kompletny Plan**

### **Krok 1: Upload Phase 2B Files**
```bash
# Na lokalnej maszynie
scp -r aws_test/ ubuntu@<instance-ip>:~/live2.0/
```

### **Krok 2: Uruchom na AWS**
```bash
# Na instancji AWS
cd live2.0/aws_test

# Uruchom wszystko
python run_phase2b_master.py --mode all
```

### **Krok 3: Monitoring**
```bash
# W drugim terminalu
python scripts/monitor_runs.py --results-dir results/phase2b_additional
```

---

## 💰 **Szacowany Koszt**

### **c6i.16xlarge**:
- **Czas**: 3-4 dni (72-96 godzin)
- **Koszt**: $180-240
- **Zalety**: Najszybsze wykonanie, stabilność

### **c6i.8xlarge**:
- **Czas**: 5-6 dni (120-144 godzin)
- **Koszt**: $150-180
- **Zalety**: Niższy koszt, nadal szybko

### **c6i.4xlarge**:
- **Czas**: 8-10 dni (192-240 godzin)
- **Koszt**: $120-150
- **Zalety**: Najniższy koszt, długi czas

---

## 🔧 **Konfiguracja dla AWS**

### **Security Group**:
- **SSH (22)**: Twój IP
- **HTTP (80)**: 0.0.0.0/0 (opcjonalne)
- **HTTPS (443)**: 0.0.0.0/0 (opcjonalne)

### **Storage**:
- **Root Volume**: 100 GB gp3
- **Dodatkowy Volume**: 500 GB gp3 (opcjonalne dla wyników)

### **Network**:
- **VPC**: Default VPC
- **Subnet**: Public subnet
- **Elastic IP**: Zalecane dla stabilności

---

## 📊 **Monitoring i Logi**

### **System Monitoring**:
```bash
# CPU i Memory
htop

# Disk usage
df -h

# Running processes
ps aux | grep python
```

### **Application Monitoring**:
```bash
# Progress tracking
tail -f results/phase2b_additional/logs/phase2b_runner.log

# Real-time dashboard
python scripts/monitor_runs.py --results-dir results/phase2b_additional
```

---

## 🎯 **Rekomendacja**

### **Dla Phase 2B zalecam: c6i.16xlarge**

**Uzasadnienie**:
1. **Wydajność**: 64 vCPUs wystarczy dla 30 symulacji równolegle
2. **Pamięć**: 128 GB RAM dla długich symulacji (500K kroków)
3. **Stabilność**: Mniej prawdopodobne zawieszenia
4. **Czas**: 3-4 dni vs 8-10 dni na mniejszej instancji
5. **Koszt**: $180-240 vs $120-150 (różnica $60-90 za oszczędność 4-6 dni)

### **Alternatywnie: c6i.8xlarge**
Jeśli budżet jest ograniczony, c6i.8xlarge też będzie działać dobrze, ale:
- Uruchamiaj po 15-20 symulacji
- Monitoruj pamięć
- Może wymagać restartów

---

## 🚀 **Gotowe do Uruchomienia**

**Po uruchomieniu instancji**:
1. Uruchom setup: `bash setup_aws_instance.sh`
2. Upload Phase 2B files
3. Uruchom: `python run_phase2b_master.py --mode all`
4. Monitoruj postęp

**Szacowany czas**: 3-4 dni  
**Szacowany koszt**: $180-240  
**Oczekiwane wyniki**: 100+ molekuł, 10+ cykli autokatalitycznych
