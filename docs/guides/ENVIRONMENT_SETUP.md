# 🚀 Środowisko Pracy Live 2.0 - Instrukcja

## ✅ Co zostało przygotowane

Twoje środowisko pracy zostało w pełni skonfigurowane:

- **Conda Environment**: `live` z Python 3.11
- **Backend**: Wszystkie zależności zainstalowane (Taichi, FastAPI, SciPy, itp.)
- **Frontend**: Node.js dependencies zainstalowane
- **Skrypty**: Automatyzacja uruchamiania i zarządzania

## 🎯 Sposób Pracy

### Aktywacja środowiska (jednorazowo)

```powershell
.\activate_live_env.ps1
```

### Uruchomienie aplikacji (2 terminale)

**Terminal 1 - Backend:**
```powershell
.\start_backend.ps1
```

**Terminal 2 - Frontend:**
```powershell
.\start_frontend.ps1
```

### Dostęp do aplikacji

- **Frontend**: http://localhost:3000
- **Backend API**: http://localhost:8000
- **API Documentation**: http://localhost:8000/docs

## 📁 Lokalizacje Środowiska

- **Conda installation**: `D:\conda`
- **Environment path**: `D:\conda_envs\live`
- **Python**: `D:\conda_envs\live\python.exe`
- **Pip**: `D:\conda_envs\live\Scripts\pip.exe`

## 🛠️ Przydatne Komendy

### Sprawdzenie instalacji
```powershell
D:\conda_envs\live\python.exe -c "import taichi; print('Taichi:', taichi.__version__)"
```

### Instalacja dodatkowych pakietów
```powershell
D:\conda_envs\live\python.exe -m pip install <package_name>
```

### Aktualizowanie requirements
```powershell
D:\conda_envs\live\python.exe -m pip install -r backend\requirements.txt
```

### Sprawdzenie portów
```powershell
netstat -ano | findstr :8000  # Backend
netstat -ano | findstr :3000  # Frontend
```

## 🔧 Rozwiązywanie Problemów

### Problem: Conda nie rozpoznawane
**Rozwiązanie**: Uruchom nową sesję PowerShell po wykonaniu `activate_live_env.ps1`

### Problem: Port zajęty
**Rozwiązanie**: 
```powershell
netstat -ano | findstr :8000
taskkill /PID <PID> /F
```

### Problem: Zależności nie zainstalowane
**Rozwiązanie**:
```powershell
D:\conda_envs\live\python.exe -m pip install -r backend\requirements.txt
```

### Problem: Frontend nie uruchamia się
**Rozwiązanie**:
```powershell
cd frontend
npm install
npm run dev
```

## 🎮 Pierwsze Uruchomienie

1. **Uruchom środowisko**: `.\activate_live_env.ps1`
2. **Uruchom backend**: `.\start_backend.ps1` (Terminal 1)
3. **Uruchom frontend**: `.\start_frontend.ps1` (Terminal 2)
4. **Otwórz przeglądarkę**: http://localhost:3000

Jeśli wszystko działa poprawnie, zobaczysz interfejs symulacji Live 2.0!

## 📞 Wsparcie

Jeśli masz problemy:
1. Sprawdź czy wszystkie skrypty są wykonywalne w PowerShell
2. Sprawdź logi błędów w terminalu
3. Sprawdź czy porty nie są zajęte
4. Restart terminali i spróbuj ponownie
