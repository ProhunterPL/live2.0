---
date: 2025-11-25
label: guide
---

# AWS vs Lokalnie: Generowanie sieci reakcji

## 📊 Analiza wymagań

### Dane do przetworzenia:
- **17 runów** hydrothermal
- **170 snapshotów** (10 na run)
- **Rozmiar snapshotu**: ~100-500 KB JSON
- **Łączny rozmiar**: ~17-85 MB danych

### Operacje:
1. **Czytanie JSON** (I/O-bound)
2. **Budowa grafów** (connected components) - O(n) gdzie n = liczba wiązań
3. **Porównywanie temporalne** (O(m²) gdzie m = liczba typów cząsteczek)
4. **Zapisywanie wyników** (I/O-bound)

### Wymagania obliczeniowe:
- **CPU**: Niskie (głównie I/O i proste operacje)
- **RAM**: ~100-500 MB na run
- **Czas**: ~1-5 minut na run (szacunkowo)
- **Łączny czas**: ~20-85 minut dla wszystkich runów

---

## 🤔 AWS vs Lokalnie

### ✅ Użyj AWS jeśli:

1. **Dane już są na AWS**
   - Nie trzeba pobierać 17-85 MB
   - Szybsze (brak transferu sieciowego)

2. **Równoległe przetwarzanie**
   - Można przetwarzać wiele runów jednocześnie
   - Przyspieszenie 4-8x (zależnie od liczby CPU cores)

3. **Lokalna maszyna jest wolna**
   - Stary komputer
   - Mało RAM
   - Wolny dysk

4. **Chcesz wykorzystać instancję**
   - Instancja już płatna
   - Nie chcesz obciążać lokalnej maszyny

### ✅ Użyj lokalnie jeśli:

1. **Dane już pobrane**
   - Mamy wszystkie snapshoty lokalnie
   - Nie trzeba transferować

2. **Szybka lokalna maszyna**
   - Nowoczesny CPU
   - SSD
   - Wystarczająco RAM

3. **Prostsze debugowanie**
   - Łatwiejszy dostęp do logów
   - Szybsze iteracje

4. **Oszczędność kosztów**
   - AWS kosztuje (nawet jeśli mała instancja)
   - Lokalnie = darmowe

---

## 💡 Rekomendacja

### **Opcja 1: Lokalnie (REKOMENDOWANE)**

**Dlaczego:**
- Dane już są lokalnie (pobrane wcześniej)
- To nie jest bardzo wymagające obliczeniowo
- Prostsze debugowanie
- Szybsze iteracje przy testowaniu

**Czas**: ~20-85 minut (można uruchomić w tle)

**Jak:**
```bash
# Stwórz batch script
python scripts/build_reaction_networks_batch.py \
    --scenario hydrothermal_extended \
    --base-dir results/phase2b_additional
```

### **Opcja 2: AWS (jeśli chcesz równoległe przetwarzanie)**

**Dlaczego:**
- Równoległe przetwarzanie (4-8 runów jednocześnie)
- Szybsze jeśli dane są już na AWS
- Wykorzystanie instancji

**Czas**: ~5-15 minut (z równoległością)

**Jak:**
```bash
# Na AWS
cd ~/live2.0
git pull origin main

# Stwórz batch script z równoległością
python scripts/build_reaction_networks_batch.py \
    --scenario hydrothermal_extended \
    --base-dir results/phase2b_additional \
    --parallel 4
```

### **Opcja 3: Hybrydowa (NAJLEPSZA)**

**Dlaczego:**
- Test lokalnie na 1-2 runach
- Pełne przetwarzanie na AWS z równoległością

**Jak:**
1. **Lokalnie**: Przetestuj na run_1
   ```bash
   python scripts/build_reaction_network_from_snapshots.py \
       --run results/phase2b_additional/hydrothermal_extended/run_1
   ```

2. **AWS**: Przetwórz wszystkie runy równolegle
   ```bash
   # Na AWS
   python scripts/build_reaction_networks_batch.py \
       --scenario hydrothermal_extended \
       --parallel 4
   ```

---

## 🚀 Plan działania

### Jeśli wybierasz LOKALNIE:

1. **Stwórz batch script** (15 min)
   ```python
   # scripts/build_reaction_networks_batch.py
   # Przetwarza wszystkie runy sekwencyjnie
   ```

2. **Uruchom** (20-85 min)
   ```bash
   python scripts/build_reaction_networks_batch.py \
       --scenario hydrothermal_extended
   ```

3. **Weryfikuj wyniki** (10 min)
   - Sprawdź wygenerowane `reaction_network.json`
   - Uruchom detektor autokatalityczny

### Jeśli wybierasz AWS:

1. **Prześlij skrypt na AWS** (5 min)
   ```bash
   scp scripts/build_reaction_networks_batch.py ubuntu@63.178.224.65:~/live2.0/scripts/
   ```

2. **Uruchom na AWS** (5-15 min z równoległością)
   ```bash
   ssh ubuntu@63.178.224.65
   cd ~/live2.0
   python scripts/build_reaction_networks_batch.py \
       --scenario hydrothermal_extended \
       --parallel 4
   ```

3. **Pobierz wyniki** (5 min)
   ```bash
   scp -r ubuntu@63.178.224.65:~/live2.0/results/phase2b_additional/hydrothermal_extended/run_*/reaction_network.json \
       results/phase2b_additional/hydrothermal_extended/
   ```

---

## 📊 Porównanie czasów

| Metoda | Czas | Zalety | Wady |
|--------|------|--------|------|
| **Lokalnie (sekwencyjnie)** | 20-85 min | Proste, darmowe | Wolniejsze |
| **Lokalnie (równolegle)** | 5-20 min | Szybkie, darmowe | Wymaga multi-core CPU |
| **AWS (sekwencyjnie)** | 20-85 min | Dane już tam | Koszt, transfer |
| **AWS (równolegle)** | 5-15 min | Najszybsze | Koszt, transfer |

---

## ✅ Moja rekomendacja

**Użyj AWS z równoległością**, jeśli:
- ✅ Instancja już jest płatna
- ✅ Chcesz szybkie wyniki (5-15 min)
- ✅ Masz doświadczenie z AWS

**Użyj lokalnie**, jeśli:
- ✅ Chcesz prostsze debugowanie
- ✅ Masz szybką lokalną maszynę
- ✅ Nie zależy Ci na czasie (można w tle)

**Hybrydowa** (najlepsza):
- ✅ Test lokalnie (run_1)
- ✅ Pełne przetwarzanie na AWS z równoległością

---

**Status**: Gotowe do implementacji  
**Czas implementacji**: 15-30 minut (batch script)  
**Czas wykonania**: 5-85 minut (zależnie od metody)

