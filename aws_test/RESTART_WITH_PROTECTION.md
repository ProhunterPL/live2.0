# ✅ Potwierdzona Diagnoza: SSH Connection Dropped

## 📊 Wyniki Diagnostyki

- ✅ **Brak OOM kills** - pamięć była OK
- ✅ **Brak błędów w logach** - symulacje działały normalnie
- ⚠️ **Procesy zostały zabite** - najprawdopodobniej przez przerwanie SSH

## 💡 Przyczyna

**SSH Connection Dropped** - połączenie zostało przerwane i procesy otrzymały sygnał SIGHUP, który je zabił.

---

## 🚀 Rozwiązanie: Uruchom Ponownie z Ochroną

### Krok 1: Sprawdź Uptime (Opcjonalnie)

```bash
uptime
last reboot
```

To potwierdzi że system się nie restartował.

### Krok 2: Uruchom w Screen (Zalecane)

```bash
# Uruchom screen
screen -S phase2b

# W screen, uruchom symulacje
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run

# Odłącz: Ctrl+A, potem D
# Podłącz ponownie: screen -r phase2b
```

**Screen chroni przed rozłączeniem SSH** - procesy będą działać nawet po zamknięciu terminala.

### Krok 3: Monitoruj Postęp

W osobnym terminalu (lub po odłączeniu od screen):

```bash
# Sprawdź czy procesy działają
ps aux | grep python | grep -v grep

# Sprawdź postęp
cd ~/live2.0
python3 aws_test/scripts/quick_diagnose.py

# Lub użyj watch do ciągłego monitorowania
watch -n 60 'ps aux | grep python | grep -v grep | wc -l'
```

---

## 📋 Pełna Instrukcja Screen

### Podstawowe Komendy Screen:

```bash
# Uruchom nowy screen
screen -S phase2b

# Lista aktywnych screenów
screen -ls

# Podłącz do screenu
screen -r phase2b

# Jeśli screen jest "attached", wymuś podłączenie
screen -d -r phase2b

# Odłącz od screenu (wewnątrz screenu)
Ctrl+A, potem D

# Wyjście z screenu (zabije procesy!)
Ctrl+A, potem K, potem Y
```

### Wewnątrz Screenu:

- **Ctrl+A, D** - odłącz (procesy działają dalej)
- **Ctrl+A, C** - nowe okno
- **Ctrl+A, N** - następne okno
- **Ctrl+A, P** - poprzednie okno
- **Ctrl+A, "** - lista okien

---

## 🔄 Alternatywa: Tmux

Jeśli wolisz tmux:

```bash
# Uruchom tmux
tmux new -s phase2b

# W tmux, uruchom symulacje
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run

# Odłącz: Ctrl+B, potem D
# Podłącz ponownie: tmux attach -t phase2b
```

---

## ⚠️ Ważne Uwagi

1. **Zawsze używaj screen/tmux** - bez tego procesy zginą przy rozłączeniu SSH
2. **Nie zamykaj terminala** - użyj Ctrl+A, D żeby odłączyć się od screenu
3. **Sprawdzaj regularnie** - używaj `screen -r phase2b` żeby sprawdzić postęp
4. **Backup logów** - przed restartem zapisz logi dla analizy

---

## 🎯 Szybki Start

```bash
# 1. Uruchom screen
screen -S phase2b

# 2. W screen, uruchom symulacje
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run

# 3. Odłącz: Ctrl+A, potem D

# 4. Sprawdź status (w osobnym terminalu)
screen -r phase2b  # Podłącz żeby zobaczyć output
# Lub
ps aux | grep python | grep -v grep  # Sprawdź czy działają
```

---

## 📊 Oczekiwany Czas

- **500K kroków** × **~9.5 kroków/sekundę** = **~14-15 godzin** na symulację
- **30 symulacji** równolegle = **~14-15 godzin** total (jeśli wszystkie działają równolegle)

**Monitoruj regularnie** żeby upewnić się że działają!

