# ⚡ FAST MODE - Phase 2B Production Run

## 🎯 Nowe Podejście

**Problem**: Novelty detection jest zbyt wolny (10 min/krok)  
**Rozwiązanie**: Wyłącz detektor podczas symulacji, uruchom offline batch

---

## ✅ Jak To Działa

### 1. Symulacja (FAST - ~1-2 godziny)
- ✅ Novelty detection: **KOMPLETNIE WYŁĄCZONE**
- ✅ Tylko fizyka symulacji
- ✅ Zapis snapshotów co 10K kroków
- ⚡ **100x szybciej** niż z detektorem

### 2. Post-Process Batch (offline)
- Uruchom analizę offline po symulacji
- Można uruchomić równolegle na wielu CPU cores
- Można uruchomić na taniej maszynie CPU
- Można uruchomić w tle gdy komputer nie jest używany

---

## 🚀 Uruchomienie

### Krok 1: Uruchom Symulację (FAST MODE)

```powershell
# Uruchom z FAST MODE (bez detektora)
python scripts/run_phase2_full.py `
  --config aws_test/configs/phase2_miller_urey_extended_FAST.yaml `
  --output results/phase2b_local/miller_urey/run_01 `
  --steps 500000 `
  --seed 100
```

**Czas**: ~1-2 godziny (vs 35 dni z detektorem!)

### Krok 2: Po Ukończeniu - Batch Analysis

```powershell
# Uruchom offline analizę na wszystkich snapshotach
python scripts/post_detect_batch.py `
  --dir results/phase2b_local/miller_urey/run_01 `
  --parallel 16
```

**Czas**: ~30 minut (16 równoległych workerów)

### Krok 3: Aggreguj Wyniki

```powershell
# Stwórz plik molecules.json z wszystkich wyników
python scripts/aggregate_post_detect.py `
  --dir results/phase2b_local/miller_urey/run_01/post_detect `
  --output results/phase2b_local/miller_urey/run_01/molecules.json
```

---

## 📊 Porównanie

| Tryb | Novelty Detection | Czas 500K kroków | Detekcja |
|------|------------------|------------------|----------|
| **Original** | Co 500 kroków | 35 dni | ❌ Zbyt wolne |
| **Optimized** | Co 10K kroków | 35 dni | ❌ Wciąż wolne |
| **FAST** | Wyłączone | **1-2 godziny** | ✅ Offline batch |

---

## 🎯 Zalety FAST MODE

### ✅ Szybkość
- 500K kroków: **1-2 godziny** (vs 35 dni)
- Batch 10 symulacji: **10-20 godzin** (vs 350 dni!)
- Batch 30 symulacji: **30-60 godzin** (1-2.5 dni)

### ✅ Elastyczność
- Możesz uruchomić analizę offline później
- Możesz użyć więcej CPU cores
- Możesz użyć taniej maszyny CPU
- Możesz uruchomić w tle

### ✅ Skalowalność
- Równoległe przetwarzanie (16-32 workers)
- Możesz przerwać i wznowić
- Możesz użyć AWS dla batch analysis

---

## 🔧 Batch Analysis na Jeden Raz

```powershell
# 1. Uruchom wszystkie 30 symulacji (bez detektora)
python run_phase2b_local.py --all --runs 10

# 2. Po ukończeniu - batch analysis wszystkich wyników
Get-ChildItem -Recurse results/phase2b_local -Filter "snapshots" | ForEach-Object {
    python scripts/post_detect_batch.py --dir $_.FullName --parallel 16
}
```

**Total**: 30-60 godzin symulacji + ~2-4 godziny batch analysis = **1-3 dni**

---

## 💡 Zaawansowane: GNU Parallel (Linux/Mac)

```bash
# Najpierw stwórz listę snapshotów
find results/phase2b_local -name "state_*.npz" > snapshot_list.txt

# Batch process z GNU parallel
cat snapshot_list.txt | \
  parallel -j 32 \
  'python scripts/post_detect_batch.py --input {} --output {}.json'
```

**32 równoległych workers** = szybka analiza!

---

## 🎉 Podsumowanie

### ✅ NOWE: FAST MODE
- ⚡ **100x szybciej** niż z detektorem
- 📊 **Te same wyniki** (offline analysis)
- 🔧 **Elastyczne** (możesz uruchomić batch później)
- 💰 **Tanie** (użyj lokalnego CPU dla batch)

### 🚀 Uruchom

```powershell
# 1. FAST MODE symulacja
python scripts/run_phase2_full.py --config aws_test/configs/phase2_miller_urey_extended_FAST.yaml --output results/phase2b_local/miller_urey/run_01 --steps 500000 --seed 100

# 2. Batch analysis po ukończeniu
python scripts/post_detect_batch.py --dir results/phase2b_local/miller_urey/run_01 --parallel 16
```

**Szacowany czas**: 1-2 godziny symulacji + 30 minut batch analysis = **2-3 godziny total** 🎉

---

*Czas to pieniądz - FAST MODE to jest to!* ⚡

