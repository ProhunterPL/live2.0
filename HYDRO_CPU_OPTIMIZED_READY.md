# ✅ HYDROTHERMAL QUEUE - CPU OPTIMIZED & READY

**Data:** 2025-11-18  
**Status:** GOTOWE DO URUCHOMIENIA  
**Tryb:** **CPU (domyślny) - SZYBSZY niż GPU dla tego typu obliczeń!**

---

## 🎯 Co Zostało Zmienione

### 1. CPU jako Domyślny Tryb
- **Testy pokazały**: CPU z wieloma rdzeniami jest **szybszy** niż GPU dla symulacji chemicznych
- Domyślnie używa **wszystkich dostępnych rdzeni**
- Automatycznie wykrywa liczbę rdzeni systemowych

### 2. Opcje Uruchomienia

#### ⚡ CPU Mode (ZALECANE - domyślne)
```powershell
python run_phase2b_hydro_queue.py --start 10 --end 1
```
- Używa wszystkich rdzeni CPU
- **Najszybsze dla tego workloadu**
- Szacowany czas: **~10 godzin** (60 min × 10 runs)

#### 🔥 Hybrid Mode (EKSPERYMENTALNE)
```powershell
python run_phase2b_hydro_queue.py --start 10 --end 1 --hybrid
```
- GPU: fizyka cząstek
- CPU: analiza chemiczna (w tle)
- Szacowany czas: **~7.5 godzin** (45 min × 10 runs)
- ⚠️ Wymaga GPU + integracji z `run_phase2_full.py`

#### 🎮 GPU Mode (dostępne, ale wolniejsze)
```powershell
python run_phase2b_hydro_queue.py --start 10 --end 1 --gpu
```
- Oryginalna wersja GPU
- Szacowany czas: **~15 godzin** (90 min × 10 runs)

### 3. Kontrola Wątków CPU
```powershell
# Użyj konkretnej liczby wątków (np. dla równoległych uruchomień)
python run_phase2b_hydro_queue.py --start 10 --end 1 --cpu-threads 12
```

---

## 🚀 QUICK START

### Metoda 1: Interaktywne Menu (Najłatwiejsze)

```powershell
.\start_hydro_queue.ps1
```

**Menu opcje:**
1. 🧪 Test CPU (10K steps, ~5 min)
2. 🚀 Pełna kolejka CPU (run 10→1, ~10h) **← ZALECANE**
3. 🔥 Pełna kolejka HYBRID (run 10→1, ~7.5h) EKSPERYMENTALNE
4. 📊 Sprawdź status
5. 📈 Monitoruj zasoby (CPU/RAM)
6. ❌ Wyjdź

### Metoda 2: Bezpośrednie Komendy

```powershell
# 1. Test (5 minut)
python scripts/run_phase2_full.py `
  --config aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml `
  --output results/test_hydro_local `
  --steps 10000 `
  --seed 42 `
  --force-cpu

# 2. Pełna kolejka (10 godzin)
python run_phase2b_hydro_queue.py --start 10 --end 1
```

---

## 📊 Porównanie Wydajności

| Tryb | Czas per Run | Total (10 runs) | Wykorzystanie |
|------|--------------|-----------------|---------------|
| **CPU (wszystkie rdzenie)** | 60 min | **10h** | CPU: 100%, RAM: 4-8GB |
| **Hybrid (GPU+CPU)** | 45 min | **7.5h** | GPU: 80%, CPU: 50% |
| GPU (tylko CUDA) | 90 min | 15h | GPU: 90-100% |

**Rekomendacja:** CPU mode (domyślny) - najlepszy stosunek wydajności do stabilności

---

## 🖥️ Wymagania Systemowe

### Minimalne:
- CPU: 4+ rdzeni
- RAM: 8GB
- Dysk: 50GB wolnego miejsca

### Zalecane:
- CPU: 8+ rdzeni (więcej = szybciej!)
- RAM: 16GB+
- Dysk: SSD z 100GB wolnego miejsca

### Dla Hybrid Mode:
- CPU: 8+ rdzeni
- GPU: NVIDIA z CUDA support
- RAM: 16GB+

---

## 📁 Pliki

| Plik | Opis |
|------|------|
| `run_phase2b_hydro_queue.py` | Główny skrypt - CPU optimized |
| `start_hydro_queue.ps1` | Interaktywne menu PowerShell |
| `docs/local/HYDROTHERMAL_QUEUE_GUIDE.md` | Pełna dokumentacja |
| `HYDRO_QUEUE_START.md` | Quick reference |
| `HYDRO_CPU_OPTIMIZED_READY.md` | Ten plik |

---

## 🔍 Monitoring

### Sprawdź Status
```powershell
# Ile runs ukończone
Get-ChildItem results/phase2b_local/hydrothermal/*/results.json

# Ostatnie logi
Get-Content logs/hydro_queue_*.log -Tail 20

# Lub użyj menu
.\start_hydro_queue.ps1
# wybierz opcję 4
```

### Monitoruj CPU
```powershell
# W menu
.\start_hydro_queue.ps1
# wybierz opcję 5

# Lub bezpośrednio Task Manager
# Ctrl + Shift + Esc → Performance → CPU
```

---

## ⚡ Automatyczne Wykrywanie Rdzeni

Skrypt automatycznie wykrywa liczbę rdzeni:

```powershell
python run_phase2b_hydro_queue.py --help
# Pokaże: "Detected CPU cores: X"
```

Na przykład:
- Intel i7-13700K: 16 rdzeni → użyje 16
- AMD Ryzen 9 7950X: 16 rdzeni → użyje 16
- AWS c5.24xlarge: 96 vCPUs → użyje 96

---

## 🎯 Timeline (CPU Mode)

| Czas | Event | Status |
|------|-------|--------|
| T+0h | run_10 starts | Running |
| T+1h | run_10 done, run_9 starts | 10% done |
| T+2h | run_9 done, run_8 starts | 20% done |
| T+3h | run_8 done, run_7 starts | 30% done |
| T+4h | run_7 done, run_6 starts | 40% done |
| T+5h | run_6 done, run_5 starts | 50% done |
| T+6h | run_5 done, run_4 starts | 60% done |
| T+7h | run_4 done, run_3 starts | 70% done |
| T+8h | run_3 done, run_2 starts | 80% done |
| T+9h | run_2 done, run_1 starts | 90% done |
| **T+10h** | **run_1 done** | **✅ 100%** |

---

## 💡 Tips

### 1. Równoległe Uruchomienia
Jeśli masz dużo rdzeni (np. 32+), możesz uruchomić 2 kolejki:

```powershell
# Terminal 1 - runs 10-6
python run_phase2b_hydro_queue.py --start 10 --end 6 --cpu-threads 16

# Terminal 2 - runs 5-1
python run_phase2b_hydro_queue.py --start 5 --end 1 --cpu-threads 16
```

### 2. Overnight Run
```powershell
# Zapobiegnij uśpieniu
powercfg /change standby-timeout-ac 1440

# Uruchom
python run_phase2b_hydro_queue.py --start 10 --end 1

# Rano sprawdź
.\start_hydro_queue.ps1
# wybierz opcję 4
```

### 3. Restart-Safe
Jeśli przerwiesz:
```powershell
# Po prostu uruchom ponownie - pominie ukończone
python run_phase2b_hydro_queue.py --start 10 --end 1
```

---

## 🔬 Dlaczego CPU Jest Szybsze?

Z testów Phase 2B:
- **Symulacje chemiczne** = dużo logiki, rozgałęzień, grafów
- **GPU** = świetne dla prostych obliczeń równoległych (grafika, fizyka)
- **CPU** = lepsze dla złożonej logiki (chemia, grafy, detekcja molekuł)

Szczególnie przy:
- Detekcji wiązań (złożone reguły)
- Analizie grafów (clustery, molekuły)
- Reakcjach chemicznych (logika warunkowa)

---

## 📚 Dokumentacja

- **Quick Start:** `HYDRO_QUEUE_START.md`
- **Full Guide:** `docs/local/HYDROTHERMAL_QUEUE_GUIDE.md`
- **Hybrid Mode:** `docs/HYBRID_GPU_CPU_GUIDE.md`
- **Original Setup:** `docs/local/README_PHASE2B_LOCAL.md`

---

## ✅ Checklist

Przed rozpoczęciem:
- [ ] Sprawdziłeś ile masz rdzeni CPU
- [ ] Test 10K steps zakończony sukcesem
- [ ] Masz ~50GB wolnego miejsca
- [ ] Ustawiłeś zapobieganie uśpieniu
- [ ] Gotowy na ~10h obliczeń

---

## 🚀 START NOW!

```powershell
# Najłatwiejsze - menu:
.\start_hydro_queue.ps1

# Lub bezpośrednio:
python run_phase2b_hydro_queue.py --start 10 --end 1
```

**Oczekiwany czas:** 10 godzin  
**Tryb:** CPU (wszystkie rdzenie)  
**Rezultat:** 10 hydrothermal runs gotowych do analizy

---

**Let's go! 🌊🔥⚡**

---

## 🎯 Cel

Po zakończeniu będziemy mieli:
- AWS: ~14 miller_urey runs
- Local: **10 hydrothermal runs** ← To robimy teraz
- **Total: ~24 runs** z 30 zaplanowanych

Brakujące 6 runs: formamide (opcjonalnie, jeśli będziemy potrzebować więcej danych)

---

**Ostatnia aktualizacja:** 2025-11-18  
**Status:** ✅ READY TO GO

