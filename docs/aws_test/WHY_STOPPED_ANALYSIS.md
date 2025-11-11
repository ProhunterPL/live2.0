# 🔍 Analiza: Dlaczego Symulacje Się Zatrzymały

## 📊 Diagnoza

Z logów wynika, że:

1. ✅ **Brak błędów w logach** - symulacje działały normalnie
2. ⚠️ **Logi kończą się nagle** - procesy zostały zabite
3. ⏰ **Zatrzymanie o ~15:00** - oba procesy zatrzymały się w tym samym czasie

## 💡 Najbardziej Prawdopodobne Przyczyny

### 1. **SSH Connection Dropped** (Najbardziej prawdopodobne)
- Połączenie SSH zostało przerwane
- Procesy zostały zabite przez system (SIGHUP)
- **Rozwiązanie**: Użyj `screen`, `tmux` lub `nohup`

### 2. **Out of Memory (OOM) Killer**
- System zabrał pamięć
- OOM killer zabił procesy Python
- **Sprawdź**: `dmesg | grep -i oom`

### 3. **System Reboot/Instance Restart**
- Instancja AWS została zrestartowana
- **Sprawdź**: `uptime`, `last reboot`

### 4. **Resource Limits**
- Procesy przekroczyły limity CPU/pamięci
- **Sprawdź**: `ulimit -a`

---

## ✅ Rozwiązanie: Uruchom Ponownie z Ochroną

### Opcja 1: Screen (Najprostsze)

```bash
# Uruchom screen
screen -S phase2b

# W screen, uruchom symulacje
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run

# Odłącz: Ctrl+A, potem D
# Podłącz ponownie: screen -r phase2b
```

### Opcja 2: Tmux

```bash
# Uruchom tmux
tmux new -s phase2b

# W tmux, uruchom symulacje
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run

# Odłącz: Ctrl+B, potem D
# Podłącz ponownie: tmux attach -t phase2b
```

### Opcja 3: Nohup (Background)

```bash
cd ~/live2.0/aws_test
nohup python3 run_phase2b_master.py --mode run > phase2b_run.log 2>&1 &

# Sprawdź status
tail -f phase2b_run.log
```

### Opcja 4: Systemd Service (Najlepsze dla produkcji)

Stwórz plik `/etc/systemd/system/phase2b.service`:

```ini
[Unit]
Description=Phase 2B Simulations
After=network.target

[Service]
Type=simple
User=ubuntu
WorkingDirectory=/home/ubuntu/live2.0/aws_test
ExecStart=/usr/bin/python3 run_phase2b_master.py --mode run
Restart=on-failure
RestartSec=10
StandardOutput=append:/home/ubuntu/live2.0/aws_test/phase2b_service.log
StandardError=append:/home/ubuntu/live2.0/aws_test/phase2b_service.log

[Install]
WantedBy=multi-user.target
```

Potem:
```bash
sudo systemctl daemon-reload
sudo systemctl enable phase2b
sudo systemctl start phase2b
sudo systemctl status phase2b
```

---

## 🔍 Sprawdź Co Się Stało

```bash
# Sprawdź OOM kills
dmesg | grep -i "killed process" | tail -20

# Sprawdź kiedy system się restartował
uptime
last reboot

# Sprawdź limity zasobów
ulimit -a

# Sprawdź użycie pamięci
free -h
```

---

## 🎯 Zalecane Działanie

1. **Sprawdź przyczynę**:
```bash
cd ~/live2.0
python3 aws_test/scripts/check_system_kills.py
```

2. **Uruchom ponownie w screen**:
```bash
screen -S phase2b
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode run
```

3. **Monitoruj**:
```bash
# W osobnym terminalu
watch -n 60 'ps aux | grep python | grep -v grep | wc -l'
```

---

## ⚠️ Ważne

- **Screen/tmux są konieczne** - bez nich procesy zginą przy rozłączeniu SSH
- **Sprawdź pamięć** - jeśli OOM, zwiększ rozmiar instancji AWS
- **Monitoruj** - sprawdzaj czy procesy działają

