# 🚀 Jak Uruchomić Live 2.0 Lokalnie

## 📋 Co Masz Gotowe ✅
- ✅ **Docker** (v28.4.0) - Zainstalowany
- ✅ **Docker Compose** (v2.39.4) - Zainstalowany  
- ✅ **Node.js** - Zainstalowany
- ✅ **npm** - Zainstalowany

## ❌ Co Musisz Zainstalować
- ❌ **Python 3.9+** - Wymagany dla backend

## 🎯 Najłatwiejszy Sposób: Docker

### Krok 1: Uruchom Docker Desktop
1. Otwórz Docker Desktop
2. Poczekaj aż status będzie "Running"
3. Sprawdź: `docker --version`

### Krok 2: Uruchom Aplikację
```powershell
# W katalogu projektu
docker compose up -d --build
```

### Krok 3: Sprawdź Status
```powershell
docker compose ps
```

### Krok 4: Otwórz Aplikację
- **Frontend**: http://localhost:3000
- **Backend**: http://localhost:8000
- **API Docs**: http://localhost:8000/docs

## 🏠 Alternatywa: Lokalne Środowisko

### Krok 1: Zainstaluj Python
1. Idź na https://www.python.org/downloads/
2. Pobierz Python 3.11
3. **WAŻNE**: Zaznacz "Add Python to PATH"
4. Zainstaluj i restart terminala

### Krok 2: Sprawdź Instalację
```powershell
python --version
pip --version
```

### Krok 3: Uruchom Backend
```powershell
cd backend
pip install -r requirements.txt
python -m api.server
```

### Krok 4: Uruchom Frontend (nowe okno)
```powershell
cd frontend
npm install
npm run dev
```

## 🔧 Rozwiązywanie Problemów

### Problem: Docker nie odpowiada
**Rozwiązanie:**
1. Uruchom Docker Desktop
2. Poczekaj na "Running" status
3. Restart terminala

### Problem: Python nie znaleziony
**Rozwiązanie:**
1. Zainstaluj Python z python.org
2. Zaznacz "Add Python to PATH"
3. Restart terminala

### Problem: Port zajęty
**Rozwiązanie:**
```powershell
# Sprawdź co używa portu
netstat -ano | findstr :8000
netstat -ano | findstr :3000

# Zabij proces
taskkill /PID <PID> /F
```

## 🚀 Szybki Start (Zalecany)

**Opcja 1: Docker (Najłatwiejsza)**
```powershell
# 1. Uruchom Docker Desktop
# 2. W katalogu projektu:
docker compose up -d --build
# 3. Otwórz: http://localhost:3000
```

**Opcja 2: Lokalnie**
```powershell
# 1. Zainstaluj Python
# 2. Backend:
cd backend && pip install -r requirements.txt && python -m api.server
# 3. Frontend (nowe okno):
cd frontend && npm install && npm run dev
# 4. Otwórz: http://localhost:3000
```

## 📊 Sprawdzenie Czy Działa

### Test Backend:
```powershell
curl http://localhost:8000/
# Powinno zwrócić: {"message": "Live 2.0 Simulation API"}
```

### Test Frontend:
- Otwórz http://localhost:3000
- Powinieneś zobaczyć interfejs symulacji
- Symulacja powinna się automatycznie uruchomić

## 🎮 Co Zobaczysz

Po uruchomieniu zobaczysz:
- **Heatmapę cząstek** w czasie rzeczywistym
- **Metryki symulacji** (liczba cząstek, wiązań, novelty rate)
- **Kontrolki** (start, pause, stop, reset)
- **Wykrywanie nowych substancji**
- **Podgląd struktur molekularnych**

## 🆘 Jeśli Nic Nie Działa

1. **Sprawdź logi Docker:**
   ```powershell
   docker compose logs
   ```

2. **Restart wszystkiego:**
   ```powershell
   docker compose down
   docker compose up -d --build
   ```

3. **Sprawdź porty:**
   ```powershell
   netstat -ano | findstr :3000
   netstat -ano | findstr :8000
   ```

4. **Użyj lokalnego środowiska** jako backup

## 🎯 Zalecenie

**Dla pierwszego uruchomienia:** Użyj Docker - jest najłatwiejszy i nie wymaga instalacji Python.

**Dla rozwoju:** Użyj lokalnego środowiska - szybsze debugowanie.
