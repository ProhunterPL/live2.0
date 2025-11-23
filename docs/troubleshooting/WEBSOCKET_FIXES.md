# WebSocket Reconnection Fixes

## Problem
WebSocket klient próbował się nieprzerwanie reconnectować do nieistniejących symulacji, co powodowało:
- Ciągłe błędy 403/1008 w logach backendu
- Niepotrzebne zużycie zasobów CPU i pamięci
- Wielokrotne próby połączenia co 30-40 sekund

## Rozwiązanie

### 1. Inteligentne zarządzanie reconnection (frontend/src/lib/ws.ts)

**Zmiany:**
- Zmniejszono maksymalną liczbę prób reconnection z 5 do 3
- Dodano flagę `shouldReconnect` zapobiegającą reconnection po permanentnych błędach
- Dodano wykrywanie kodów błędów oznaczających brak symulacji:
  - **1008** - Policy Violation (backend zwraca to dla nieistniejącej symulacji)
  - **1003** - Unsupported Data
  - **1002** - Protocol Error

**Nowe eventy:**
- `simulation_not_found` - emitowany gdy backend zwróci kod 1008 (symulacja nie istnieje)
- `reconnect_failed` - emitowany gdy przekroczono limit prób reconnection

### 2. Automatyczna rekreacja symulacji (frontend/src/App.tsx)

**Zmiany:**
- Dodano nasłuchiwanie eventu `simulation_not_found`
- Przy wykryciu braku symulacji: automatycznie tworzona jest nowa symulacja
- Wyczyszczenie stanu gdy symulacja nie istnieje w `updateStatus()`

### 3. Inteligentne polling metryk (frontend/src/components/Controls.tsx)

**Zmiany:**
- Dodano detekcję błędów 404 w `updateMetrics()` i `updateNovelSubstances()`
- Polling przestaje działać gdy wykryje że symulacja nie istnieje
- Logowanie ostrzeżeń zamiast ciągłych błędów

## Rezultaty

### Przed zmianami:
```
2025-10-01 22:21:48,148 - Simulation sim_1759349761253 not found
INFO: connection rejected (403 Forbidden)
// ... powtarzane co ~30 sekund w nieskończoność
```

### Po zmianach:
```
WebSocket closed with code 1008 (simulation not found or permanent error). 
Stopping reconnection attempts.
Simulation not found on backend. Creating new simulation...
```

## Parametry reconnection

- **Max retry attempts**: 3 (zmniejszone z 5)
- **Retry delay**: 1000ms * numer próby (1s, 2s, 3s)
- **Permanent failure codes**: 1008, 1003, 1002
- **Normal closure code**: 1000 (nie powoduje reconnection)

## Użycie

### Normalne użycie
Wszystkie zmiany działają automatycznie:
1. Frontend łączy się z backendem
2. Jeśli symulacja nie istnieje → automatycznie tworzy nową
3. Jeśli połączenie się zerwie → próbuje maksymalnie 3 razy
4. Jeśli backend zwróci "nie ma symulacji" → przestaje próbować i tworzy nową

### Manualne zatrzymanie reconnection
```typescript
wsClient.stopReconnecting() // Zatrzymuje próby reconnection bez zamykania połączenia
wsClient.disconnect()        // Zamyka połączenie i zatrzymuje reconnection
```

## Testing

1. **Test 1: Restart backendu**
   ```bash
   # Uruchom frontend i backend
   .\start_frontend.ps1
   .\start_backend.ps1
   
   # Otwórz frontend w przeglądarce
   # Następnie zabij backend
   .\kill_backend.ps1
   
   # Uruchom backend ponownie
   .\start_backend.ps1
   
   # Odśwież przeglądarkę
   # Rezultat: Frontend automatycznie utworzy nową symulację
   ```

2. **Test 2: Sprawdzanie logów**
   ```bash
   # Po restarcie backendu, nie powinno być wielokrotnych błędów:
   # ❌ BAD: "Simulation sim_XXX not found" (powtarzane w kółko)
   # ✅ GOOD: "Stopping reconnection attempts" (raz, następnie nowa symulacja)
   ```

## Dodatkowe usprawnienia

Jeśli problem nadal występuje, możliwe dalsze usprawnienia:
1. Dodać mechanizm rate limiting w backendzie
2. Dodać cache sprawdzania aktywnych symulacji
3. Dodać mechanizm "circuit breaker" w frontendzie
4. Persystencja symulacji w backendzie (pliki/baza danych)

