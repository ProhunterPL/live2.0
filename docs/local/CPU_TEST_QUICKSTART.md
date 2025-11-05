# 🔧 CPU Test - Szybki Start

## 🛑 Krok 1: Zatrzymaj Obecną Symulację

W terminalu gdzie symulacja działa:
```powershell
# Ctrl+C
```

---

## 🚀 Krok 2: Uruchom CPU Test (10K kroków)

```powershell
python run_phase2_cpu_test.py `
  --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml `
  --output results/test_cpu_10k `
  --steps 10000 `
  --seed 42
```

**Czas**: ~10-20 minut  
**Cel**: Sprawdzić czy CPU jest szybszy niż GPU

---

## 📊 Krok 3: Porównaj Wyniki

Po ukończeniu testu sprawdź logi:

```powershell
# CPU test
Get-Content results/test_cpu_10k/simulation.log -Tail 5

# GPU test (jeśli miałeś)
Get-Content results/test_gpu_perf/simulation.log -Tail 5
```

---

## ✅ Jeśli CPU Jest Szybszy

Użyj CPU dla pełnej symulacji:

```powershell
python run_phase2_cpu_test.py `
  --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml `
  --output results/phase2b_local/miller_urey/cpu_run_01 `
  --steps 500000 `
  --seed 100
```

---

## ⚠️ Jeśli CPU Jest Wolniejszy

Może problem jest gdzie indziej:
- Zbyt wiele particles
- Zbyt częste operacje
- Inne bottlenecki

W takim przypadku rozważ:
- Zmniejszyć particles jeszcze bardziej (500 zamiast 1000)
- Zwiększyć timestep (0.05 zamiast 0.01)
- Użyć AWS z większą instancją

---

**Zacznij od CPU test - to tylko 10-20 minut!** ⚡

