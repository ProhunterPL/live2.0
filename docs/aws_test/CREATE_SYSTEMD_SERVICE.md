# 🔧 Jak Stworzyć Systemd Service dla Phase 2B

## ⚠️ Problem

Wklejenie zawartości pliku bezpośrednio do terminala nie działa - trzeba stworzyć plik.

---

## ✅ Rozwiązanie: Użyj Skryptu (Najprostsze)

```bash
cd ~/live2.0/aws_test
chmod +x scripts/create_phase2b_service.sh
sudo ./scripts/create_phase2b_service.sh
```

To automatycznie stworzy plik service i pokaże następne kroki.

---

## 🔧 Alternatywa: Ręczne Tworzenie

### Krok 1: Stwórz Plik Service

```bash
sudo nano /etc/systemd/system/phase2b.service
```

### Krok 2: Wklej Zawartość

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

### Krok 3: Zapisz i Wyjdź

- **Nano**: Ctrl+X, potem Y, potem Enter
- **Vim**: Esc, potem `:wq`, potem Enter

### Krok 4: Aktywuj Service

```bash
sudo systemctl daemon-reload
sudo systemctl enable phase2b
sudo systemctl start phase2b
sudo systemctl status phase2b
```

---

## 📋 Komendy do Zarządzania Service

```bash
# Sprawdź status
sudo systemctl status phase2b

# Zatrzymaj
sudo systemctl stop phase2b

# Uruchom
sudo systemctl start phase2b

# Restart
sudo systemctl restart phase2b

# Wyłącz auto-start
sudo systemctl disable phase2b

# Włącz auto-start
sudo systemctl enable phase2b

# Zobacz logi
sudo journalctl -u phase2b -f

# Lub sprawdź plik log
tail -f ~/live2.0/aws_test/phase2b_service.log
```

---

## ✅ Po Uruchomieniu

Sprawdź czy działa:

```bash
# Sprawdź status service
sudo systemctl status phase2b

# Sprawdź procesy
ps aux | grep python | grep run_phase2b

# Sprawdź logi
tail -f ~/live2.0/aws_test/phase2b_service.log
```

---

## 🎯 Zalecane: Użyj Skryptu

Najprostsze rozwiązanie:

```bash
cd ~/live2.0/aws_test
chmod +x scripts/create_phase2b_service.sh
sudo ./scripts/create_phase2b_service.sh

# Potem
sudo systemctl daemon-reload
sudo systemctl enable phase2b
sudo systemctl start phase2b
sudo systemctl status phase2b
```

---

## ⚠️ Ważne

1. **Service automatycznie restartuje** jeśli proces zginie
2. **Działa nawet po rozłączeniu SSH** - nie potrzebujesz screen
3. **Logi są w pliku** - sprawdzaj `phase2b_service.log`
4. **Auto-start przy boot** - jeśli włączysz `enable`

