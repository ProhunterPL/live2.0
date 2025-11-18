# 🚀 Quick Start - Hydrothermal Queue (CPU OPTIMIZED)

## TL;DR

```powershell
# 1. Test CPU (5 min)
.\start_hydro_queue.ps1
# Wybierz opcję 1

# 2. Jeśli test OK, uruchom pełną kolejkę CPU (10h)
python run_phase2b_hydro_queue.py --start 10 --end 1
```

**⚡ CPU jest SZYBSZE niż GPU dla tego workloadu!**

## Co To Robi?

- Uruchamia symulacje **hydrothermal** lokalnie na **CPU**
- Kolejność: **run_10 → run_9 → ... → run_1** (od końca)
- Równolegle z AWS (które robi miller_urey)
- **Cel**: Przyspieszyć progres Phase 2B
- **Mode**: CPU z wszystkimi rdzeniami (szybsze niż GPU!)

## Pliki

| Plik | Opis |
|------|------|
| `start_hydro_queue.ps1` | Menu startowe (CPU/Hybrid/GPU opcje) |
| `run_phase2b_hydro_queue.py` | Główny skrypt - CPU optimized |
| `HYDRO_CPU_OPTIMIZED_READY.md` | Setup summary & performance info |
| `docs/local/HYDROTHERMAL_QUEUE_GUIDE.md` | Pełna dokumentacja |

## Szybkie Komendy

```powershell
# Test CPU
python scripts/run_phase2_full.py --config aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml --output results/test_hydro_local --steps 10000 --seed 42 --force-cpu

# Pełna kolejka CPU (domyślne, najszybsze)
python run_phase2b_hydro_queue.py --start 10 --end 1

# Hybrid mode (eksperymentalne, jeszcze szybsze)
python run_phase2b_hydro_queue.py --start 10 --end 1 --hybrid

# GPU mode (jeśli chcesz porównać)
python run_phase2b_hydro_queue.py --start 10 --end 1 --gpu

# Status
Get-ChildItem results/phase2b_local/hydrothermal/*/results.json

# Sprawdź CPU
Get-Process | Where-Object {$_.ProcessName -eq "python"} | Select-Object CPU
```

## Timeline

### CPU Mode (domyślny)
| Run | Czas | Łącznie |
|-----|------|---------|
| run_10 | 60 min | 1h |
| run_9 | 60 min | 2h |
| run_8 | 60 min | 3h |
| ... | ... | ... |
| run_1 | 60 min | **10h** |

### Hybrid Mode (eksperymentalne)
| Run | Czas | Łącznie |
|-----|------|---------|
| run_10 | 45 min | 0.75h |
| ... | ... | ... |
| run_1 | 45 min | **7.5h** |

## Monitoring

```powershell
# Logi w czasie rzeczywistym
Get-Content logs/hydro_queue_*.log -Wait -Tail 20

# CPU usage
Get-Process | Where-Object {$_.ProcessName -eq "python"} | Select-Object CPU, WorkingSet

# Lub użyj menu
.\start_hydro_queue.ps1
# wybierz opcję 5
```

## Bezpieczeństwo

- **Restart-safe**: Jeśli się przerwie, uruchom ponownie - pominie ukończone
- **Ctrl+C**: Bezpieczne przerwanie
- **Auto-skip**: Pomija już ukończone runs

## Wyniki

```
results/phase2b_local/hydrothermal/
├── run_10/
│   ├── results.json          ✅ Główne wyniki
│   ├── molecules.json        ✅ Molekuły
│   └── snapshots/            ✅ 10 snapshots
├── run_09/
└── ... (10 runs total)
```

## Tryby Uruchomienia

| Tryb | Komenda | Czas | Zalecenie |
|------|---------|------|-----------|
| **CPU** | `python run_phase2b_hydro_queue.py --start 10 --end 1` | 10h | ⚡ **ZALECANE** |
| **Hybrid** | `python run_phase2b_hydro_queue.py --start 10 --end 1 --hybrid` | 7.5h | 🔥 Eksperymentalne |
| GPU | `python run_phase2b_hydro_queue.py --start 10 --end 1 --gpu` | 15h | Wolniejsze |

## Pełna Dokumentacja

- **Setup Summary:** `HYDRO_CPU_OPTIMIZED_READY.md`
- **Full Guide:** `docs/local/HYDROTHERMAL_QUEUE_GUIDE.md`
- **Hybrid Mode:** `docs/HYBRID_GPU_CPU_GUIDE.md`

---

**Ready? Let's go! 🌊🔥⚡**

```powershell
.\start_hydro_queue.ps1
```

**💡 Tip: CPU jest szybsze niż GPU dla chemii!**

