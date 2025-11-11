# 📥 Pobieranie i Analiza Wyników Phase 2B z AWS

**Status**: ✅ Phase 2B zakończone na AWS  
**Data**: 5 listopada 2025  
**Czas wykonania**: ~1-2 dni (SUPER FAST MODE)

---

## 🎉 **Podsumowanie Uruchomienia**

### ✅ **Ukończone:**
- ✅ Formamide debug (9 testów)
- ✅ 30 dodatkowych symulacji (500K kroków każda)
- ✅ Analiza wyników
- ✅ Raporty wygenerowane

### 📊 **Raporty na AWS:**
- `results/phase2b_additional/phase2b_summary_report.md`
- `results/phase2b_additional/phase2b_analysis_report.md`
- `results/phase2b_additional/formamide_debug/formamide_debug_report.md`

---

## 📥 **Krok 1: Pobierz Wyniki z AWS**

### **Opcja A: Użyj skryptu downloader (zalecane)**

```bash
# Na lokalnej maszynie
python aws_test/scripts/download_phase2b_results.py \
    --host <AWS-INSTANCE-IP> \
    --key <path-to-key.pem> \
    --local-dir results/phase2b_aws_results
```

### **Opcja B: Ręczne pobieranie przez SCP**

```bash
# Pobierz wszystkie wyniki
scp -r -i <key.pem> \
    ubuntu@<AWS-INSTANCE-IP>:~/live2.0/aws_test/results/phase2b_additional \
    results/phase2b_aws_results

# Lub tylko raporty (szybsze)
scp -r -i <key.pem> \
    ubuntu@<AWS-INSTANCE-IP>:~/live2.0/aws_test/results/phase2b_additional/*.md \
    results/phase2b_aws_results/
```

### **Opcja C: Sprawdź status przed pobraniem**

```bash
# Sprawdź status bez pobierania
python aws_test/scripts/download_phase2b_results.py \
    --host <AWS-INSTANCE-IP> \
    --key <path-to-key.pem> \
    --status-only
```

---

## 📊 **Krok 2: Przeczytaj Raporty**

Po pobraniu sprawdź raporty:

```bash
# Summary report
cat results/phase2b_aws_results/phase2b_summary_report.md

# Analysis report
cat results/phase2b_aws_results/phase2b_analysis_report.md

# Formamide debug report
cat results/phase2b_aws_results/formamide_debug/formamide_debug_report.md
```

---

## 🔍 **Krok 3: Analiza Offline - Novelty Detection**

Wszystkie symulacje były w SUPER FAST MODE, więc novelty detection był wyłączony.  
**Musisz uruchomić analizę offline na snapshotach:**

### **Automatyczna analiza dla wszystkich uruchomień:**

```bash
# Dla każdego scenariusza
for scenario in miller_urey_extended hydrothermal_extended formamide_extended; do
    echo "Analizuję $scenario..."
    for run_dir in results/phase2b_aws_results/$scenario/run_*; do
        if [ -d "$run_dir/snapshots" ]; then
            echo "  -> $run_dir"
            python scripts/post_detect_batch.py --input "$run_dir" --parallel 4
        fi
    done
done
```

### **Lub pojedynczo:**

```bash
# Analiza jednego uruchomienia
python scripts/post_detect_batch.py \
    --input results/phase2b_aws_results/miller_urey_extended/run_01 \
    --parallel 4

# Analiza wszystkich snapshotów w katalogu
python scripts/post_detect_batch.py \
    --input results/phase2b_aws_results/miller_urey_extended \
    --parallel 8
```

---

## 📈 **Krok 4: Agregacja Wyników**

Po analizie offline możesz zaggregować wyniki:

```bash
# Jeśli masz skrypt agregujący
python scripts/analyze_batch_results.py \
    --input-dir results/phase2b_aws_results \
    --output results/phase2b_final_analysis.json
```

---

## 📊 **Krok 5: Sprawdź Strukturę Wyników**

```bash
# Sprawdź strukturę
tree -L 3 results/phase2b_aws_results/ | head -50

# Lub
find results/phase2b_aws_results -name "*.json" -type f | wc -l
find results/phase2b_aws_results -name "snapshots" -type d | wc -l
```

---

## 🎯 **Oczekiwane Wyniki**

### **Po analizie offline powinieneś mieć:**

| Scenariusz | Uruchomień | Oczekiwane molekuły |
|------------|-----------|---------------------|
| Miller-Urey | 10 | 30-50 na run |
| Hydrothermal | 10 | 20-40 na run |
| Formamide | 10 | 10-30 na run |
| **Total** | **30** | **600-1200 unikalnych** |

### **Autocatalytic cycles:**
- Oczekiwane: 5-20 cykli na scenariusz
- Total: 15-60 cykli

---

## 🔧 **Rozwiązywanie Problemów**

### **Problem: Wyniki nie zostały pobrane**

```bash
# Sprawdź połączenie
ssh -i <key.pem> ubuntu@<AWS-INSTANCE-IP> "ls -la ~/live2.0/aws_test/results/phase2b_additional"

# Sprawdź czy pliki istnieją
ssh -i <key.pem> ubuntu@<AWS-INSTANCE-IP> "find ~/live2.0/aws_test/results/phase2b_additional -name '*.md'"
```

### **Problem: Analiza offline nie działa**

```bash
# Sprawdź czy snapshoty istnieją
find results/phase2b_aws_results -name "snapshots" -type d

# Sprawdź zawartość snapshotów
ls results/phase2b_aws_results/miller_urey_extended/run_01/snapshots/ | head -10
```

### **Problem: Brak pamięci podczas analizy**

```bash
# Uruchom mniej równoległych procesów
python scripts/post_detect_batch.py --input <dir> --parallel 2

# Lub analizuj pojedynczo
python scripts/post_detect_batch.py --input <dir> --parallel 1
```

---

## 📋 **Checklist Po Pobraniu**

- [ ] Wyniki pobrane z AWS
- [ ] Raporty przeczytane
- [ ] Analiza offline uruchomiona dla wszystkich uruchomień
- [ ] Wyniki zaggregowane
- [ ] Statystyki obliczone
- [ ] Wykresy wygenerowane (jeśli potrzebne)
- [ ] Instancja AWS zatrzymana (aby oszczędzić koszty)

---

## 💰 **Koszt AWS**

### **Szacowany koszt:**
- Instancja c6i.16xlarge: ~$2.50/godzinę
- Czas wykonania: 1-2 dni (24-48 godzin)
- **Total koszt: $60-120** ⬇️ (znacznie niższy niż poprzednio!)

### **Zatrzymaj instancję po pobraniu wyników:**

```bash
# Sprawdź ID instancji
aws ec2 describe-instances --filters "Name=tag:Name,Values=live2-phase2b-optimized"

# Zatrzymaj instancję (oszczędza koszty)
aws ec2 stop-instances --instance-ids <instance-id>

# Lub usuń jeśli nie potrzebujesz
aws ec2 terminate-instances --instance-ids <instance-id>
```

---

## 🎉 **Następne Kroki**

1. ✅ Pobierz wyniki z AWS
2. ✅ Uruchom analizę offline na snapshotach
3. ✅ Zaggreguj wyniki
4. ✅ Wygeneruj raporty końcowe
5. ✅ Przejdź do Phase 3 (Paper Writing)

---

**Gotowe do pobrania wyników!** 🚀

