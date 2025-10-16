# ⚡ SZYBKIE KOMENDY AWS - Przepisz po kolei

## 📝 UWAGA: Włącz wklejanie w PowerShell
- **Kliknij prawym przyciskiem myszy** w oknie = wkleja tekst
- Lub użyj: `Shift + Insert`

---

## 🔧 INSTALACJA (Krok po kroku)

### 1️⃣ Zainstaluj narzędzia (już masz update ✓)
```bash
sudo apt install -y python3.11 python3-pip git htop
```

### 2️⃣ Pobierz kod
```bash
git clone https://github.com/ProhunterPL/live2.0.git
cd live2.0
```

### 3️⃣ Zainstaluj Python packages
```bash
pip3 install --upgrade pip
pip3 install -r requirements.txt
```

### 4️⃣ Sprawdź instalację
```bash
python3 -c "import taichi; print('OK')"
```

---

## ✅ TEST WYDAJNOŚCI

```bash
python3 scripts/run_phase2_full.py --config configs/phase2_quick_test.yaml --output results/test1 --steps 1000 --seed 42
```

Sprawdź wynik (powinna być prędkość 4-6 steps/s):
```bash
cat results/test1/summary.txt
```

---

## 🚀 URUCHOMIENIE PRODUKCYJNE

### Dla c6i.16xlarge (64 CPU) - 16 równoległych zadań:

```bash
nohup python3 scripts/phase2_master_1M.py --mode full --scenarios all --max-parallel 16 > production.log 2>&1 &
```

### Monitoruj:
```bash
tail -f production.log
```

### Lub w drugim terminalu (Ctrl+C aby wyjść):
```bash
htop
```

---

## 📊 SPRAWDŹ POSTĘP

```bash
# Ile symulacji zakończonych
find results -name "summary.txt" | wc -l

# Ostatnie 20 linii logu
tail -20 production.log

# Zużycie dysku
du -sh results/
```

---

## 💾 POBIERANIE WYNIKÓW (z lokalnego Windows)

W nowym PowerShell (nie na AWS):
```powershell
scp -i twoj-klucz.pem -r ubuntu@TWOJ-IP:~/live2.0/results ./aws_results
```

---

## ⏹️ ZATRZYMAJ SYMULACJE (jeśli potrzeba)

```bash
pkill -f phase2_master
```

---

## 🎯 SUPER SZYBKA WERSJA (wszystko naraz)

Jeśli uda Ci się wkleić, to całość w jednym kawałku:

```bash
sudo apt install -y python3.11 python3-pip git htop && \
git clone https://github.com/ProhunterPL/live2.0.git && \
cd live2.0 && \
pip3 install --upgrade pip && \
pip3 install -r requirements.txt && \
python3 -c "import taichi; print('Installation OK!')" && \
echo "Ready to run tests!"
```

Potem test:
```bash
cd ~/live2.0
python3 scripts/run_phase2_full.py --config configs/phase2_quick_test.yaml --output results/test1 --steps 1000 --seed 42
```

Potem produkcja:
```bash
cd ~/live2.0
nohup python3 scripts/phase2_master_1M.py --mode full --scenarios all --max-parallel 16 > prod.log 2>&1 &
tail -f prod.log
```

---

## 🆘 PROBLEMY?

### Nie działa wklejanie w PowerShell?
1. Spróbuj Windows Terminal: https://aka.ms/terminal
2. Lub PuTTY: https://putty.org
3. Lub użyj VS Code z Remote SSH

### Błąd przy instalacji?
```bash
# Sprawdź logi
cat /var/log/apt/term.log
```

### Nie ma phase2_master_1M.py?
```bash
cd ~/live2.0
git pull
ls scripts/phase2_master_1M.py
```

---

## ⏱️ CZASY (dla c6i.16xlarge)

- Instalacja: **5-10 minut**
- Test (1000 kroków): **3-5 minut**  
- Jedna symulacja (200k): **5-7 godzin**
- 150 symulacji (16 równoległych): **~2.4 dni**
- Koszt: **~$157**

---

**Powodzenia! 🚀**

Jeśli masz problemy, napisz gdzie się zatrzymałeś!

