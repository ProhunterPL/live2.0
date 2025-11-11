# 🔧 Instrukcje Pobierania Wyników Phase 2B z AWS

## ❌ Problem

Skrypt `download_phase2b_results.py` **musi być uruchomiony na lokalnej maszynie Windows**, nie na AWS!

Na AWS wyniki już są gotowe w:
```
~/live2.0/aws_test/results/phase2b_additional/
```

---

## ✅ Rozwiązanie

### **Opcja 1: Użyj skryptu downloader (z lokalnej maszyny Windows)**

```powershell
# Na lokalnej maszynie Windows (PowerShell)
python aws_test\scripts\download_phase2b_results.py `
    --host 63.178.224.65 `
    --key "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    --status-only
```

Jeśli status OK, pobierz wyniki:

```powershell
python aws_test\scripts\download_phase2b_results.py `
    --host 63.178.224.65 `
    --key "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    --local-dir results\phase2b_aws_results
```

### **Opcja 2: Bezpośrednie pobieranie przez SCP (z lokalnej maszyny Windows)**

```powershell
# Pobierz wszystkie wyniki
scp -r -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    ubuntu@63.178.224.65:~/live2.0/aws_test/results/phase2b_additional `
    results\phase2b_aws_results
```

### **Opcja 3: Tylko raporty (szybsze, mniejsze)**

```powershell
# Najpierw stwórz katalog
mkdir results\phase2b_aws_results

# Pobierz tylko pliki MD (raporty)
scp -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    ubuntu@63.178.224.65:~/live2.0/aws_test/results/phase2b_additional/*.md `
    results\phase2b_aws_results\

# Pobierz strukturę katalogów (bez zawartości snapshotów)
scp -r -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    ubuntu@63.178.224.65:~/live2.0/aws_test/results/phase2b_additional/*/run_*/summary.txt `
    results\phase2b_aws_results\
```

---

## 🔍 Sprawdzenie Statusu na AWS

Jeśli chcesz sprawdzić status na AWS, użyj:

```bash
# Na AWS (SSH)
cd ~/live2.0/aws_test/results/phase2b_additional
ls -la

# Sprawdź raporty
cat phase2b_summary_report.md | head -50

# Sprawdź liczbę ukończonych symulacji
find . -name "summary.txt" -type f | wc -l
find . -name "snapshots" -type d | wc -l
```

---

## 📊 Co Pobierać

### **Minimum (tylko raporty):**
- `phase2b_summary_report.md` - podsumowanie
- `phase2b_analysis_report.md` - analiza
- `formamide_debug_report.md` - debug formamide

### **Pełne (dla analizy offline):**
- Wszystkie katalogi `run_*/` z snapshotami
- Wszystkie pliki JSON z wynikami
- Logi symulacji

**Uwaga**: Snapshoty mogą być duże (kilka GB), więc pobieraj tylko jeśli potrzebujesz analizy offline.

---

## 🚀 Po Pobraniu - Analiza Offline

Jeśli pobrałeś snapshoty, uruchom analizę offline:

```powershell
# Dla każdego scenariusza
$scenarios = @("miller_urey_extended", "hydrothermal_extended", "formamide_extended")
foreach ($scenario in $scenarios) {
    $runs = Get-ChildItem "results\phase2b_aws_results\$scenario\run_*" -Directory
    foreach ($run in $runs) {
        Write-Host "Analizuję $($run.Name)..."
        python scripts\post_detect_batch.py --input "$($run.FullName)" --parallel 4
    }
}
```

---

## 💡 Szybkie Sprawdzenie Wyników

Najszybsze podejście - sprawdź tylko raporty na AWS:

```bash
# Na AWS (SSH)
cat ~/live2.0/aws_test/results/phase2b_additional/phase2b_summary_report.md
cat ~/live2.0/aws_test/results/phase2b_additional/phase2b_analysis_report.md | head -100
```

Następnie zdecyduj czy pobierać pełne wyniki czy tylko raporty.

