# 🔧 Rozwiązywanie Problemu z Pobieraniem Wyników

## Problem

Skrypt `download_phase2b_results.py` pokazuje:
- `Completed runs: 0`
- `Running processes: 0`
- `⚠️ Phase 2B may have issues`

Ale wcześniej widzieliśmy że Phase 2B zostało zakończone i raporty zostały wygenerowane.

---

## 🔍 Diagnostyka

### **Sprawdź ręcznie na AWS co rzeczywiście jest:**

```bash
# Na AWS (SSH)
ssh -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" ubuntu@63.178.224.65

# Sprawdź strukturę katalogów
cd ~/live2.0/aws_test/results/phase2b_additional
ls -la

# Sprawdź raporty
ls -la *.md

# Sprawdź run directories
find . -type d -name "run_*" | wc -l

# Sprawdź pliki results.json
find . -name "results.json" | wc -l

# Sprawdź pliki summary.txt
find . -name "summary.txt" | wc -l
```

---

## ✅ Rozwiązanie 1: Pobierz Bezpośrednio Raporty

Jeśli raporty istnieją, pobierz je bezpośrednio:

```powershell
# Na lokalnej maszynie Windows
scp -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    ubuntu@63.178.224.65:~/live2.0/aws_test/results/phase2b_additional/*.md `
    results\phase2b_aws_results\
```

---

## ✅ Rozwiązanie 2: Pobierz Wszystko Pomimo Statusu

Skrypt został poprawiony żeby szukał też `results.json` i raportów MD.  
Spróbuj pobrać mimo statusu "partial":

```powershell
# Na lokalnej maszynie Windows
python aws_test\scripts\download_phase2b_results.py `
    --host 63.178.224.65 `
    --key "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    --local-dir results\phase2b_aws_results
```

Skrypt spróbuje pobrać wszystko co znajdzie, nawet jeśli status pokazuje problemy.

---

## ✅ Rozwiązanie 3: Bezpośrednie Pobieranie przez SCP

```powershell
# Pobierz cały katalog results/phase2b_additional
scp -r -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    ubuntu@63.178.224.65:~/live2.0/aws_test/results/phase2b_additional `
    results\phase2b_aws_results
```

To pobierze wszystko co jest w tym katalogu, niezależnie od struktury.

---

## 📊 Sprawdź Co Jest Na AWS

Najpierw sprawdź co rzeczywiście jest na AWS:

```bash
# Na AWS (SSH)
cd ~/live2.0/aws_test/results/phase2b_additional

# Lista wszystkiego
find . -type f -name "*.md" -o -name "*.json" -o -name "summary.txt" | head -20

# Sprawdź strukturę
tree -L 3 -d . | head -50

# Sprawdź rozmiary
du -sh */run_* 2>/dev/null | head -20
```

---

## 💡 Najszybsze Rozwiązanie

Jeśli Phase 2B faktycznie się zakończyło (widziałeś raporty), najszybsze będzie:

```powershell
# 1. Pobierz tylko raporty (szybkie)
mkdir results\phase2b_aws_results
scp -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    ubuntu@63.178.224.65:~/live2.0/aws_test/results/phase2b_additional/*.md `
    results\phase2b_aws_results\

# 2. Przeczytaj raporty żeby zobaczyć co jest
cat results\phase2b_aws_results\phase2b_summary_report.md
```

Jeśli raporty mówią że wszystko jest OK, pobierz pełne wyniki:

```powershell
# 3. Pobierz pełne wyniki
scp -r -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    ubuntu@63.178.224.65:~/live2.0/aws_test/results/phase2b_additional `
    results\phase2b_aws_results
```

---

**Spróbuj najpierw sprawdzić co jest na AWS, a potem pobierz bezpośrednio przez SCP.**

