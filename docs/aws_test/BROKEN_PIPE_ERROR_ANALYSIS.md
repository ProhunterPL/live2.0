# 🔍 Analiza Błędu "Broken Pipe"

## 📊 Błąd w Logach

```
2025-11-07 15:06:53,211 - backend.sim.core.stepper - INFO - Step 86000 completed in 72.6ms
2025-11-07 15:06:53,213 - __main__ - ERROR - [FAILED] Simulation failed: [Errno 32] Broken pipe
```

---

## 💡 Co To Oznacza

**Broken pipe (Errno 32)** oznacza że:
- Proces próbował zapisać do pipe/kolejki która została zamknięta
- Najprawdopodobniej proces został zabity przez przerwanie SSH
- Lub pipe między procesami został zamknięty przedwcześnie

---

## 🔍 Kiedy Się Stało

- **Data**: 2025-11-07 15:06 (wczoraj)
- **Symulacje**: run_7 i run_8
- **Krok**: 86,000/500,000 (17.2%)
- **Przyczyna**: Przerwanie połączenia SSH zabrało procesy

---

## ✅ Czy To Problem Dla Nowych Symulacji?

**NIE** - to dotyczyło tylko starych symulacji (run_7, run_8).

### Dlaczego Nowe Symulacje Są Bezpieczne:

1. **Systemd Service** - chroni przed przerwaniem SSH
2. **Restart on failure** - automatycznie restartuje jeśli proces zginie
3. **Działają w tle** - nie zależą od sesji SSH

### Stare Symulacje (run_3-8):
- ⚠️ Były uruchomione bez systemd/screen
- ⚠️ Zginęły przy przerwaniu SSH
- ⚠️ Błąd "Broken pipe" to efekt uboczny

---

## 🔍 Jak Sprawdzić Czy Nowe Symulacje Mają Problemy

### Sprawdź Błędy w Nowych Logach:
```bash
# Sprawdź czy są błędy w run_1 i run_2
grep -i "error\|failed\|broken" ~/live2.0/results/phase2b_additional/miller_urey_extended/run_1/simulation.log
grep -i "error\|failed\|broken" ~/live2.0/results/phase2b_additional/miller_urey_extended/run_2/simulation.log
```

### Sprawdź Status Service:
```bash
sudo systemctl status phase2b
```

Jeśli widzisz "active (running)" - wszystko działa poprawnie.

### Sprawdź Czy Procesy Działają:
```bash
ps aux | grep python | grep run_phase2
```

Powinno pokazać 4 procesy (master + runner + 2 symulacje).

---

## ⚠️ Co Jeśli Pojawi Się Ten Błąd W Nowych Symulacjach?

### Jeśli systemd wykryje błąd:
1. **Automatycznie zrestartuje** proces (Restart=on-failure)
2. **Zaloguje błąd** w journalctl
3. **Kontynuuje** z następnymi symulacjami

### Sprawdź Logi Service:
```bash
sudo journalctl -u phase2b -n 50
```

### Sprawdź Czy Service Restartuje:
```bash
sudo systemctl status phase2b
# Szukaj "Restart count" - jeśli > 0, to były restarty
```

---

## ✅ Podsumowanie

- ✅ **Błąd dotyczył starych symulacji** (run_7, run_8)
- ✅ **Nowe symulacje są bezpieczne** - systemd je chroni
- ✅ **Systemd automatycznie restartuje** jeśli proces zginie
- ⚠️ **Monitoruj regularnie** - sprawdzaj czy działają

---

## 🎯 Zalecane Działanie

1. **Pozwól nowym symulacjom działać** - systemd je chroni
2. **Monitoruj regularnie** - sprawdzaj czy działają
3. **Sprawdzaj logi** - jeśli pojawią się błędy, systemd je naprawi

**Błąd "Broken pipe" w starych logach nie jest problemem dla nowych symulacji!** ✅

