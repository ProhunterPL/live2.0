---
date: 2025-11-25
label: guide
---

# Formamide - Gotowość do uruchomienia na AWS

## ✅ Status: GOTOWE DO URUCHOMIENIA

### Sprawdzone elementy:

1. **✅ Konfiguracja**: `phase2_formamide_extended_SUPER_FAST.yaml` istnieje
2. **✅ Skrypt kolejkowania**: `auto_queue_restart_formamide.sh` stworzony
3. **✅ Launcher**: `start_formamide_queue_aws.sh` stworzony
4. **✅ Skrypty przesłane na AWS**: Gotowe

---

## 🚀 Szybkie uruchomienie (1 komenda)

```bash
ssh -i "C:\Users\klawi\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem" \
    ubuntu@63.178.224.65 \
    "cd ~/live2.0 && bash aws_test/scripts/start_formamide_queue_aws.sh"
```

**To uruchomi:**
- 8 runów formamide (run_1 do run_8)
- 4 równolegle (max)
- Automatyczne kolejkowanie
- Auto-restart przy zakończeniu

---

## 📋 Szczegóły konfiguracji

### Parametry:
- **Scenariusz**: Formamide Extended
- **Liczba runów**: 8
- **Max równoległe**: 4
- **Steps per run**: 500,000
- **Seeds**: 200-207 (run_1=200, run_2=201, ...)
- **CPU threads per run**: 16 (optymalne dla 4 równoległych na 64-core)
- **Config**: `aws_test/configs/phase2_formamide_extended_SUPER_FAST.yaml`

### Oczekiwany czas:
- **Per run**: 6-8 godzin
- **Total (4 równoległe)**: ~12-16 godzin
  - Batch 1 (run_1-4): ~6-8h
  - Batch 2 (run_5-8): ~6-8h (startuje po zakończeniu batch 1)

---

## 📊 Monitorowanie

### Podstawowe komendy:

```bash
# Status kolejki
ssh -i "..." ubuntu@63.178.224.65 "tail -20 ~/live2.0/logs/auto_restart_formamide_main.log"

# Sprawdź działające procesy
ssh -i "..." ubuntu@63.178.224.65 "ps aux | grep formamide | grep run_phase2_full"

# Sprawdź postęp
ssh -i "..." ubuntu@63.178.224.65 "cd ~/live2.0 && python3 aws_test/scripts/check_actual_progress.py"
```

### Logi:
- **Główny log**: `logs/auto_restart_formamide_main.log`
- **Postęp**: `logs/auto_restart_formamide_progress.log`
- **Per run**: `logs/formamide_run_X.log`

---

## ⚙️ Konfiguracja systemu

### Auto-restart system:
- **Check interval**: 5 minut (300s)
- **Max parallel**: 4
- **Auto-restart**: Tak (automatycznie startuje następne gdy slot wolny)
- **Stuck detection**: Tak (wykrywa procesy bez aktualizacji >2h)

### CPU & Memory:
- **CPU threads per run**: 16
- **Total CPU usage**: ~64 cores (4 × 16)
- **Memory per run**: ~4-5 GB
- **Total memory**: ~16-20 GB

---

## 🔄 Workflow

1. **Uruchomienie**:
   ```bash
   bash aws_test/scripts/start_formamide_queue_aws.sh
   ```

2. **System automatycznie**:
   - Startuje run_1-4 (4 równolegle)
   - Monitoruje postęp co 5 minut
   - Gdy run kończy → startuje następny z kolejki
   - Kontynuuje aż wszystkie 8 runów zakończone

3. **Zakończenie**:
   - Wszystkie 8 runów zakończone
   - Log: "ALL FORMAMIDE RUNS COMPLETED!"

---

## ✅ Weryfikacja przed uruchomieniem

### Sprawdź na AWS:

```bash
# 1. Czy config istnieje
ls -lh ~/live2.0/aws_test/configs/phase2_formamide_extended_SUPER_FAST.yaml

# 2. Czy skrypty są wykonywalne
ls -lh ~/live2.0/aws_test/scripts/auto_queue_restart_formamide.sh
ls -lh ~/live2.0/aws_test/scripts/start_formamide_queue_aws.sh

# 3. Czy katalog wyników istnieje
ls -ld ~/live2.0/results/phase2b_additional/formamide_extended

# 4. Sprawdź dostępne zasoby
free -h
nproc
```

---

## 🎯 Po zakończeniu

### Pobierz wyniki:

```bash
# Z lokalnej maszyny
scp -i "C:\Users\klawi\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem" \
    -r ubuntu@63.178.224.65:~/live2.0/results/phase2b_additional/formamide_extended \
    results/phase2b_additional/
```

### Analiza:

```bash
# Ekstrakcja cząsteczek
python scripts/fix_run1_molecules.py --scenario formamide_extended

# Pełna analiza
python scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data
```

---

## ⚠️ Uwagi

1. **Nie uruchamiaj równocześnie z innymi scenariuszami** - może przeciążyć system
2. **Monitoruj pierwsze 10 minut** - upewnij się, że wszystko działa
3. **Sprawdź logi** - jeśli błędy, przerwij i zdiagnozuj

---

## 📝 Quick Reference

### Start:
```bash
bash aws_test/scripts/start_formamide_queue_aws.sh
```

### Stop:
```bash
pkill -f auto_queue_restart_formamide.sh
pkill -f "formamide.*run_phase2_full"
```

### Status:
```bash
tail -f logs/auto_restart_formamide_main.log
```

---

**Status**: ✅ GOTOWE  
**Ostatnia aktualizacja**: 2025-11-25  
**Czas uruchomienia**: ~12-16 godzin (8 runów, 4 równolegle)

