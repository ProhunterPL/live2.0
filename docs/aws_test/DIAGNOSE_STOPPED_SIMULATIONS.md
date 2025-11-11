# 🔍 Diagnostyka Zatrzymanych Symulacji Phase 2B

Symulacje się zatrzymały - brak procesów Python. Sprawdź co się stało:

## ✅ Polecenia do Uruchomienia na AWS

### 1. Użyj prostego skryptu diagnostycznego (najlepsze):
```bash
cd ~/live2.0
python3 aws_test/scripts/quick_diagnose.py
```

### 2. Sprawdź ostatnie kroki w logach (kiedy się zatrzymały):
```bash
# Sprawdź run_1
tail -20 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_1/simulation.log | grep -E "Step|ERROR|CRITICAL|Exception"

# Sprawdź run_2  
tail -20 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_2/simulation.log | grep -E "Step|ERROR|CRITICAL|Exception"
```

### 3. Sprawdź kiedy ostatnio były aktualizowane logi:
```bash
ls -lh ~/live2.0/results/phase2b_additional/miller_urey_extended/run_*/simulation.log
```

### 4. Sprawdź czy są pliki results.json (czy się ukończyły):
```bash
find ~/live2.0/results/phase2b_additional -name "results.json"
```

### 5. Sprawdź status w phase2b_results.json:
```bash
cd ~/live2.0/results/phase2b_additional
python3 -c "import json; d=json.load(open('phase2b_results.json')); print('Completed:', d.get('completed_runs', 0)); print('Failed:', d.get('failed_runs', 0)); print('Total:', d.get('total_runs', 0))"
```

### 6. Sprawdź szczegóły statusu każdego run:
```bash
cd ~/live2.0/results/phase2b_additional
python3 -c "import json; d=json.load(open('phase2b_results.json')); [print(f\"{s['scenario']} run_{r['run_id']}: {r.get('status', 'unknown')}\") for s in d.get('scenarios', {}).values() for r in s.get('runs', [])]"
```

### 7. Sprawdź logi master runnera:
```bash
tail -50 ~/live2.0/results/phase2b_additional/logs/phase2b_runner.log
```

### 8. Sprawdź czy były błędy w logach:
```bash
grep -i "error\|exception\|failed\|crashed" ~/live2.0/results/phase2b_additional/miller_urey_extended/run_*/simulation.log | tail -20
```

---

## 🎯 Najszybsze Sprawdzenie (Jedna Komenda)

```bash
cd ~/live2.0/results/phase2b_additional && \
echo "=== Ostatnie kroki run_1 ===" && \
tail -5 miller_urey_extended/run_1/simulation.log | grep "Step" && \
echo "" && \
echo "=== Ostatnie kroki run_2 ===" && \
tail -5 miller_urey_extended/run_2/simulation.log | grep "Step" && \
echo "" && \
echo "=== Kiedy ostatnio aktualizowane ===" && \
ls -lh miller_urey_extended/run_*/simulation.log && \
echo "" && \
echo "=== Pliki results.json ===" && \
find . -name "results.json" | wc -l && \
echo "" && \
echo "=== Status z JSON ===" && \
python3 -c "import json; d=json.load(open('phase2b_results.json')); print(f\"Completed: {d.get('completed_runs', 0)}/{d.get('total_runs', 0)}\"); print(f\"Failed: {d.get('failed_runs', 0)}/{d.get('total_runs', 0)}\")"
```

---

## 💡 Możliwe Scenariusze

1. **Symulacje się zakończyły** - sprawdź czy są pliki `results.json`
2. **Symulacje się zatrzymały** - sprawdź ostatnie kroki w logach i czy były błędy
3. **Symulacje się crashowały** - sprawdź logi pod kątem błędów
4. **Proces został zabity** - sprawdź czy były błędy OOM (Out of Memory) lub inne

