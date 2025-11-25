---
date: 2025-11-25
label: guide
---

# Quick Start: Generowanie sieci reakcji na AWS

## 🚀 Szybkie uruchomienie (1 komenda)

```bash
ssh -i "C:\Users\klawi\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem" \
    ubuntu@63.178.224.65 \
    "cd ~/live2.0 && bash scripts/aws_build_reaction_networks.sh"
```

**Czas**: ~5-15 minut (z równoległością 4 workers)

---

## 📋 Krok po kroku

### Krok 1: Połącz się z AWS
```bash
ssh -i "C:\Users\klawi\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem" \
    ubuntu@63.178.224.65
```

### Krok 2: Przejdź do katalogu projektu
```bash
cd ~/live2.0
```

### Krok 3: Uruchom skrypt
```bash
bash scripts/aws_build_reaction_networks.sh
```

LUB bezpośrednio:
```bash
python3 scripts/build_reaction_networks_batch.py \
    --scenario hydrothermal_extended \
    --base-dir results/phase2b_additional \
    --parallel 4
```

---

## ✅ Weryfikacja wyników

Po zakończeniu sprawdź:
```bash
# Liczba wygenerowanych plików
find results/phase2b_additional/hydrothermal_extended/run_*/reaction_network.json | wc -l
# Powinno być: 17

# Sprawdź przykładowy plik
cat results/phase2b_additional/hydrothermal_extended/run_1/reaction_network.json | head -20
```

---

## 🔄 Następne kroki

Po wygenerowaniu sieci reakcji:

1. **Uruchom detektor autokatalityczny**:
```bash
python3 scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data
```

2. **Pobierz wyniki lokalnie** (opcjonalnie):
```bash
# Z lokalnej maszyny
scp -i "C:\Users\klawi\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem" \
    -r ubuntu@63.178.224.65:~/live2.0/results/phase2b_additional/hydrothermal_extended/run_*/reaction_network.json \
    results/phase2b_additional/hydrothermal_extended/
```

---

## ⚡ Szybka komenda (wszystko w jednej linii)

```bash
ssh -i "C:\Users\klawi\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem" ubuntu@63.178.224.65 "cd ~/live2.0 && python3 scripts/build_reaction_networks_batch.py --scenario hydrothermal_extended --base-dir results/phase2b_additional --parallel 4 && echo 'Done! Check results with: find results/phase2b_additional/hydrothermal_extended/run_*/reaction_network.json | wc -l'"
```

---

**Status**: Gotowe do uruchomienia  
**Czas**: ~5-15 minut  
**Ostatnia aktualizacja**: 2025-11-25

