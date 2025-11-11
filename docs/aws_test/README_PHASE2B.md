# Phase 2B: Dodatkowe Uruchomienia - Instrukcje

**Status**: Gotowy do uruchomienia  
**Cel**: Uruchomienie 30 dodatkowych symulacji dla osiągnięcia celów Phase 2  
**Timeline**: 2-3 tygodnie

---

## 🚀 **Szybki Start**

### **Uruchomienie Wszystkiego (Zalecane)**
```bash
cd aws_test
python run_phase2b_master.py --mode all
```

### **Uruchomienie Krok po Kroku**
```bash
# 1. Debug formamide
python run_phase2b_master.py --mode debug

# 2. Uruchom 30 symulacji
python run_phase2b_master.py --mode run

# 3. Monitoruj postęp (opcjonalne)
python run_phase2b_master.py --mode monitor

# 4. Analizuj wyniki
python run_phase2b_master.py --mode analyze
```

---

## 📋 **Co Zostanie Uruchomione**

### **30 Dodatkowych Symulacji**:
- **10 Miller-Urey Extended** (500K kroków, seeds 100-109)
- **10 Hydrothermal Extended** (500K kroków, seeds 110-119)  
- **10 Formamide Extended** (500K kroków, seeds 120-129)

### **Debug Formamide**:
- **3 testy krótkie** (10K kroków)
- **3 testy średnie** (50K kroków)
- **3 testy długie** (100K kroków)

---

## 📁 **Struktura Plików**

```
aws_test/
├── configs/
│   ├── phase2_miller_urey_extended.yaml
│   ├── phase2_hydrothermal_extended.yaml
│   ├── phase2_formamide_debug.yaml
│   └── phase2_formamide_extended.yaml
├── scripts/
│   ├── run_phase2b_additional.py
│   ├── debug_formamide.py
│   ├── monitor_runs.py
│   └── analyze_additional_results.py
├── run_phase2b_master.py
└── results/phase2b_additional/
    ├── miller_urey_extended/
    ├── hydrothermal_extended/
    ├── formamide_extended/
    ├── formamide_debug/
    └── logs/
```

---

## 🎯 **Cele i Kryteria Sukcesu**

### **Główne Cele**:
- **Różnorodność molekularna**: 11 → 100+ unikalnych molekuł
- **Autocatalytic cycles**: 0 → 10+ cykli autokatalitycznych
- **Per-scenario diversity**: 8-9 → 30+ molekuł na scenariusz
- **Debug formamide**: 0 → aktywny scenariusz

### **Kryteria Sukcesu**:
- **Minimum Success**: 50+ molekuł, formamide aktywny, ≥90% completion
- **Optimal Success**: 100+ molekuł, 10+ cykli, 30+ na scenariusz, ≥95% completion

---

## ⏱️ **Szacowany Czas**

| Etap | Czas | Opis |
|------|------|------|
| Debug formamide | 2-4 godziny | 9 krótkich testów |
| 30 symulacji | 60-90 godzin | 500K kroków każda |
| Analiza | 1-2 godziny | Automatyczna analiza |
| **Total** | **3-4 dni** | Ciągłe uruchomienie |

---

## 🔍 **Monitoring i Logi**

### **Real-time Monitoring**:
```bash
python scripts/monitor_runs.py --results-dir results/phase2b_additional
```

### **Logi**:
- `results/phase2b_additional/logs/phase2b_runner.log` - Główny log
- `results/phase2b_additional/logs/phase2b_*.log` - Logi scenariuszy
- `results/phase2b_additional/formamide_debug/logs/` - Logi debug

### **Progress Tracking**:
- `results/phase2b_additional/phase2b_results.json` - Postęp w czasie rzeczywistym
- `results/phase2b_additional/phase2b_summary_report.md` - Raport końcowy

---

## 📊 **Analiza Wyników**

### **Automatyczna Analiza**:
```bash
python scripts/analyze_additional_results.py
```

### **Generowane Raporty**:
- `phase2b_analysis_report.md` - Kompletna analiza
- `phase2b_analysis_results.json` - Dane JSON
- `formamide_debug_report.md` - Raport debug formamide

### **Porównanie z Phase 2A**:
- Różnorodność molekularna (11 vs nowe)
- Aktywność scenariuszy
- Stabilność symulacji
- Rekomendacje

---

## 🛠️ **Troubleshooting**

### **Problem: Formamide nadal nieaktywny**
```bash
# Sprawdź debug report
cat results/phase2b_additional/formamide_debug/formamide_debug_report.md

# Uruchom dodatkowy debug
python scripts/debug_formamide.py --test long_test
```

### **Problem: Symulacje się zawieszają**
```bash
# Sprawdź logi
tail -f results/phase2b_additional/logs/phase2b_runner.log

# Sprawdź system
python scripts/monitor_runs.py --report-only
```

### **Problem: Niska różnorodność molekularna**
```bash
# Sprawdź konfiguracje
cat configs/phase2_*_extended.yaml

# Uruchom analizę
python scripts/analyze_additional_results.py
```

---

## 📈 **Oczekiwane Wyniki**

### **Po Ukończeniu**:
- **Total molecules**: 50-150 (vs obecne 11)
- **Autocatalytic cycles**: 5-20 (vs obecne 0)
- **Formamide active**: 10-30 molekuł (vs obecne 0)
- **Completion rate**: ≥90%

### **Dla Publikacji**:
- Wystarczające dane do napisania papera
- Solidne wyniki naukowe
- Kompletne Phase 2

---

## 🎯 **Następne Kroki**

### **Po Ukończeniu Phase 2B**:
1. **Przejrzyj raporty** analizy
2. **Wygeneruj figury** publikacyjne
3. **Przejdź do Phase 3** (Paper Writing)
4. **Przygotuj manuskrypt** do publikacji

### **Timeline do Publikacji**:
- **Tydzień 1-2**: Phase 2B uruchomienia
- **Tydzień 3**: Analiza i figury
- **Tydzień 4-6**: Pisanie papera
- **Tydzień 7**: Submission

---

## 📞 **Wsparcie**

### **W przypadku problemów**:
1. Sprawdź logi w `results/phase2b_additional/logs/`
2. Uruchom monitoring: `python scripts/monitor_runs.py`
3. Sprawdź raporty debug
4. Skontaktuj się z zespołem

### **Przydatne komendy**:
```bash
# Sprawdź status
python scripts/monitor_runs.py --report-only

# Sprawdź postęp
cat results/phase2b_additional/phase2b_results.json

# Sprawdź logi
tail -f results/phase2b_additional/logs/phase2b_runner.log
```

---

**Gotowy do uruchomienia!** 🚀

Uruchom: `python run_phase2b_master.py --mode all`
