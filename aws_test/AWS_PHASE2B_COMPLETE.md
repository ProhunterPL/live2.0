# AWS Phase 2B - Kompletny Plan Uruchomienia

**Data**: 24 października 2025  
**Status**: ✅ **GOTOWE DO URUCHOMIENIA NA AWS**

---

## 🖥️ **Rekomendowana Instancja AWS**

### **c6i.16xlarge (Zalecana)**
- **vCPUs**: 64
- **RAM**: 128 GB  
- **Storage**: 100 GB SSD
- **Koszt**: ~$2.50/godzina
- **Szacowany koszt**: $180-240 (3-4 dni)

### **Alternatywy**:
- **c6i.8xlarge**: 32 vCPUs, 64 GB RAM, ~$1.25/godzina
- **c6i.4xlarge**: 16 vCPUs, 32 GB RAM, ~$0.63/godzina

---

## 🚀 **Kroki Uruchomienia**

### **1. Uruchom Instancję AWS**
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

### **3. Setup Instancji**
```bash
# Pobierz i uruchom setup
wget https://raw.githubusercontent.com/ProhunterPL/live2.0/main/setup_aws_instance.sh
bash setup_aws_instance.sh
```

### **4. Upload Phase 2B Files**
```bash
# Na lokalnej maszynie
scp -r aws_test/ ubuntu@<instance-ip>:~/live2.0/
```

### **5. Uruchom Phase 2B**
```bash
# Na instancji AWS
cd live2.0/aws_test
bash run_phase2b_aws.sh
```

---

## 📊 **Co Zostanie Uruchomione**

### **Debug Formamide (9 testów)**:
- 3 testy krótkie (10K kroków)
- 3 testy średnie (50K kroków)
- 3 testy długie (100K kroków)

### **30 Dodatkowych Symulacji**:
- **10 Miller-Urey Extended** (500K kroków, seeds 100-109)
- **10 Hydrothermal Extended** (500K kroków, seeds 110-119)
- **10 Formamide Extended** (500K kroków, seeds 120-129)

---

## 🔍 **Monitoring**

### **Real-time Monitoring**:
```bash
# Na instancji AWS
python scripts/monitor_runs.py --results-dir results/phase2b_additional
```

### **System Monitoring**:
```bash
# CPU i Memory
htop

# Disk usage
df -h

# Running processes
ps aux | grep python
```

---

## 📥 **Pobieranie Wyników**

### **Sprawdź Status**:
```bash
# Na lokalnej maszynie
python scripts/download_phase2b_results.py --host <instance-ip> --key <key.pem> --status-only
```

### **Pobierz Wyniki**:
```bash
# Na lokalnej maszynie
python scripts/download_phase2b_results.py --host <instance-ip> --key <key.pem> --local-dir results/phase2b_local
```

### **Analizuj Wyniki**:
```bash
# Na lokalnej maszynie
python scripts/download_phase2b_results.py --analyze-only --local-dir results/phase2b_local
```

---

## 📁 **Struktura Wyników**

```
results/phase2b_additional/
├── miller_urey_extended/
│   ├── run_1/
│   ├── run_2/
│   └── ... (10 runs)
├── hydrothermal_extended/
│   ├── run_1/
│   ├── run_2/
│   └── ... (10 runs)
├── formamide_extended/
│   ├── run_1/
│   ├── run_2/
│   └── ... (10 runs)
├── formamide_debug/
│   ├── short_test/
│   ├── medium_test/
│   └── long_test/
├── logs/
│   ├── phase2b_runner.log
│   └── phase2b_*.log
├── phase2b_summary_report.md
├── phase2b_analysis_report.md
└── formamide_debug_report.md
```

---

## 🎯 **Oczekiwane Wyniki**

### **Po Ukończeniu**:
- **Total molecules**: 50-150 (vs obecne 11)
- **Autocatalytic cycles**: 5-20 (vs obecne 0)
- **Formamide active**: 10-30 molekuł (vs obecne 0)
- **Completion rate**: ≥90%

### **Dla Publikacji**:
- ✅ Wystarczające dane do napisania papera
- ✅ Solidne wyniki naukowe
- ✅ Kompletne Phase 2

---

## ⏱️ **Timeline**

| Etap | Czas | Opis |
|------|------|------|
| Setup instancji | 30 min | Instalacja i konfiguracja |
| Debug formamide | 2-4 godziny | 9 krótkich testów |
| 30 symulacji | 60-90 godzin | 500K kroków każda |
| Analiza | 1-2 godziny | Automatyczna analiza |
| **Total** | **3-4 dni** | Ciągłe uruchomienie |

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

## 🛠️ **Troubleshooting**

### **Problem: Instancja się zawiesza**
```bash
# Sprawdź system
htop
free -h
df -h

# Restart jeśli potrzeba
sudo reboot
```

### **Problem: Symulacje się zawieszają**
```bash
# Sprawdź logi
tail -f results/phase2b_additional/logs/phase2b_runner.log

# Sprawdź procesy
ps aux | grep python
```

### **Problem: Brak miejsca na dysku**
```bash
# Sprawdź miejsce
df -h

# Usuń stare logi
rm -rf logs/*.log.old
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

---

## 🚀 **Gotowe do Uruchomienia**

**Po uruchomieniu instancji**:
1. ✅ Uruchom setup: `bash setup_aws_instance.sh`
2. ✅ Upload Phase 2B files
3. ✅ Uruchom: `bash run_phase2b_aws.sh`
4. ✅ Monitoruj postęp

**Szacowany czas**: 3-4 dni  
**Szacowany koszt**: $180-240  
**Oczekiwane wyniki**: 100+ molekuł, 10+ cykli autokatalitycznych

---

## 📞 **Wsparcie**

### **W przypadku problemów**:
1. Sprawdź logi w `results/phase2b_additional/logs/`
2. Uruchom monitoring: `python scripts/monitor_runs.py`
3. Sprawdź system: `htop`, `free -h`, `df -h`
4. Skontaktuj się z zespołem

### **Przydatne komendy**:
```bash
# Sprawdź status
python scripts/download_phase2b_results.py --host <ip> --key <key> --status-only

# Sprawdź postęp
tail -f results/phase2b_additional/logs/phase2b_runner.log

# Sprawdź system
htop
```

---

**Gotowe do uruchomienia na AWS!** 🚀

**Rekomendacja**: c6i.16xlarge, 3-4 dni, $180-240
