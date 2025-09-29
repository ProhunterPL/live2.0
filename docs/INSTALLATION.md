# 🚀 Przewodnik Uruchamiania Live 2.0 Lokalnie

## 📋 Wymagania Systemowe

### ✅ Co już masz zainstalowane:
- ✅ **Docker** (v28.4.0) - Gotowy do użycia!
- ✅ **Node.js** - Gotowy do użycia!
- ✅ **npm** - Gotowy do użycia!

### ❌ Co musisz zainstalować:
- ❌ **Python 3.9+** - Wymagany dla backend
- ❌ **pip** - Menedżer pakietów Pythona

## 🛠️ Instalacja Wymaganych Narzędzi

### 1. Instalacja Python (Windows)

**Opcja A: Z oficjalnej strony**
1. Idź na https://www.python.org/downloads/
2. Pobierz Python 3.9+ (zalecane 3.11)
3. Podczas instalacji **ZAZNACZ** "Add Python to PATH"
4. Zainstaluj

**Opcja B: Przez Chocolatey (jeśli masz)**
```powershell
choco install python
```

**Opcja C: Przez winget (Windows 10+)**
```powershell
winget install Python.Python.3.11
```

### 2. Weryfikacja instalacji
```powershell
python --version
pip --version
```

## 🎯 Sposoby Uruchomienia

### Opcja 1: Docker (NAJŁATWIEJSZA) 🐳

**Zalety:**
- ✅ Nie wymaga instalacji Python
- ✅ Automatyczna konfiguracja
- ✅ Izolowane środowisko
- ✅ Gotowe do produkcji

**Uruchomienie:**
```powershell
# W katalogu projektu
docker-compose up -d

# Sprawdź status
docker-compose ps

# Zatrzymaj
docker-compose down
```

**Dostęp:**
- Frontend: http://localhost:3000
- Backend: http://localhost:8000
- API Docs: http://localhost:8000/docs

### Opcja 2: Lokalne Środowisko 🏠

**Backend:**
```powershell
# Przejdź do katalogu backend
cd backend

# Zainstaluj zależności
pip install -r requirements.txt

# Uruchom serwer
python -m api.server
```

**Frontend (nowe okno terminala):**
```powershell
# Przejdź do katalogu frontend
cd frontend

# Zainstaluj zależności
npm install

# Uruchom serwer dev
npm run dev
```

### Opcja 3: Makefile (jeśli masz make) 🔧

```powershell
# Instaluj wszystko
make install

# Uruchom w trybie dev
make dev

# Uruchom testy
make test
```

## 🔧 Rozwiązywanie Problemów

### Problem: "Python nie znaleziony"
**Rozwiązanie:**
1. Zainstaluj Python z https://python.org
2. Podczas instalacji zaznacz "Add Python to PATH"
3. Restart terminala

### Problem: "pip nie znaleziony"
**Rozwiązanie:**
```powershell
# Jeśli masz Python ale nie pip
python -m ensurepip --upgrade
```

### Problem: "Taichi nie może znaleźć GPU"
**Rozwiązanie:**
- Taichi automatycznie przełączy się na CPU
- Symulacja będzie działać, ale wolniej
- Dla GPU: zainstaluj CUDA Toolkit

### Problem: "Port już zajęty"
**Rozwiązanie:**
```powershell
# Sprawdź co używa portu
netstat -ano | findstr :8000
netstat -ano | findstr :3000

# Zabij proces (zastąp PID)
taskkill /PID <PID> /F
```

## 📊 Sprawdzenie Instalacji

### Test Backend:
```powershell
curl http://localhost:8000/
# Powinno zwrócić: {"message": "Live 2.0 Simulation API", "version": "1.0.0"}
```

### Test Frontend:
```powershell
curl http://localhost:3000/
# Powinno zwrócić HTML strony
```

### Test WebSocket:
```powershell
# Użyj narzędzia jak Postman lub przeglądarki
# Połącz się z: ws://localhost:8000/simulation/{id}/stream
```

## 🎮 Pierwsze Uruchomienie

1. **Uruchom aplikację** (Docker lub lokalnie)
2. **Otwórz przeglądarkę** na http://localhost:3000
3. **Aplikacja automatycznie:**
   - Utworzy symulację
   - Połączy się z backend
   - Rozpocznie symulację
4. **Zobaczysz:**
   - Heatmapę cząstek
   - Metryki w czasie rzeczywistym
   - Kontrolki symulacji

## 🚀 Zalecane Uruchomienie

**Dla początkujących:** Docker
```powershell
docker-compose up -d
```

**Dla deweloperów:** Lokalne środowisko
```powershell
# Terminal 1 - Backend
cd backend && pip install -r requirements.txt && python -m api.server

# Terminal 2 - Frontend  
cd frontend && npm install && npm run dev
```

## 📞 Wsparcie

Jeśli masz problemy:
1. Sprawdź logi: `docker-compose logs`
2. Sprawdź status: `docker-compose ps`
3. Restart: `docker-compose restart`
4. Pełny restart: `docker-compose down && docker-compose up -d`
