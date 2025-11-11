# 🔧 Naprawa Timeout - Phase 2B

## ❌ Problem

Symulacje są zabijane przez timeout po 6 godzinach:
```
⏰ miller_urey_extended run 1 timed out after 6 hours
```

**Przyczyna**: Timeout był ustawiony na 6 godzin, ale symulacje 500K kroków potrzebują **10-14 godzin**.

## ✅ Rozwiązanie

Zwiększono timeout z **6 godzin** do **24 godzin** w `aws_test/scripts/run_phase2b_additional.py`.

## 🚀 Co Zrobić Na AWS

### Opcja 1: Zatrzymaj i uruchom ponownie (zalecane)

```bash
# 1. Zatrzymaj obecne symulacje
pkill -f run_phase2_full.py
pkill -f run_phase2b_additional.py

# 2. Zaktualizuj kod
cd ~/live2.0
git pull

# 3. Uruchom ponownie
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run
```

### Opcja 2: Pozostaw obecne i uruchom nowe (szybsze)

```bash
# 1. Zaktualizuj kod (nie zatrzymuj obecnych symulacji)
cd ~/live2.0
git pull

# 2. Obecne symulacje (run_5, run_6) mogą się nie ukończyć (stary timeout)
# Ale nowe symulacje będą miały 24h timeout

# 3. Jeśli obecne się nie ukończą, uruchom ponownie tylko te które failed:
cd ~/live2.0/aws_test
python3 scripts/run_phase2b_additional.py \
  --output-dir results/phase2b_additional \
  --max-parallel 2 \
  --scenario miller_urey_extended
```

## 📊 Status Obecnych Symulacji

Z diagnostyki:
- **run_1**: Zatrzymana po 6h (Step 214,000/500,000 = 42.8%)
- **run_2**: Zatrzymana po 6h (Step 72,000/500,000 = 14.4%)
- **run_3**: Zatrzymana po 6h (Step 185,000/500,000 = 37.0%)
- **run_4**: Zatrzymana po 6h (Step 185,000/500,000 = 37.0%)
- **run_5**: Działa (Step 72,000/500,000 = 14.4%) - może się nie ukończyć (stary timeout)
- **run_6**: Działa (Step 72,000/500,000 = 14.4%) - może się nie ukończyć (stary timeout)

## 💡 Rekomendacja

**Zatrzymaj obecne symulacje i uruchom ponownie** - obecne (run_5, run_6) prawdopodobnie nie ukończą się przed timeoutem 6h.

Po zaktualizowaniu kodu wszystkie nowe symulacje będą miały 24h timeout, co wystarczy na ukończenie 500K kroków.

## ⏱️ Szacowany Czas

- **Jedna symulacja**: ~10-14 godzin (500K kroków)
- **Z timeoutem 24h**: Wystarczająco dużo czasu
- **30 symulacji ÷ 2 równolegle**: ~150-210 godzin (~6-9 dni)

