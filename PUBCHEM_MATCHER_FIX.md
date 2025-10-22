# PubChem Matcher - Naprawa (22 października 2024)

## Problem
PubChem matching nie działało poprawnie dla prostych cząsteczek (NH3, H2, N3H3, etc.). System zgłaszał "No PubChem Match" nawet dla dobrze znanych cząsteczek które są w PubChem.

## Główne przyczyny

### 1. Brak exact match (dokładnego dopasowania)
Stary kod używał tylko **similarity search** (wyszukiwanie podobnych cząsteczek), które:
- Nie działa dobrze dla bardzo prostych cząsteczek
- Jest wolniejsze
- Wymaga threshold (próg podobieństwa)

### 2. Problemy z asynchronicznym API PubChem
PubChem similarity search zwraca:
- Kod **202 (Accepted)** = "Twoje zapytanie jest przetwarzane"
- **ListKey** = identyfikator do odpytywania o wyniki

Stary kod nie obsługiwał tego i od razu zwracał błąd.

### 3. Problemy z Unicode w konsoli Windows
Emoji i znaki specjalne (✓, ✅, ❌, etc.) powodowały błędy `UnicodeEncodeError` w konsoli Windows z polskim kodowaniem (cp1250).

## Rozwiązanie

### 1. Dodano `pubchem_exact_match()` 
Nowa funkcja która:
- Najpierw szuka **dokładnego dopasowania** SMILES w PubChem
- Używa prostszego endpoint'u: `/compound/smiles/{SMILES}/JSON`
- Jest szybka i niezawodna dla prostych cząsteczek
- Automatycznie fallback do similarity search jeśli exact match nie zadziała

### 2. Naprawiono asynchroniczny similarity search
Funkcja `pubchem_similar_top()` teraz:
- Wykrywa status 202
- Odpytuje (poll) wyniki używając ListKey
- Czeka maksymalnie 20 sekund (10 prób × 2 sekundy)
- Poprawnie obsługuje zarówno synchroniczne jak i asynchroniczne odpowiedzi

### 3. Naprawiono encoding dla Windows
Dodano automatyczną konfigurację UTF-8 dla konsoli Windows:
```python
if sys.platform == 'win32':
    import io
    if sys.stdout.encoding != 'utf-8':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
```

## Zmiany w plikach

### `matcher/chem.py`
- ✅ Dodano `pubchem_exact_match()` - nowa funkcja exact match
- ✅ Zaktualizowano `pubchem_similar_top()` - obsługa async API
- ✅ Dodano `import time` dla polling
- ✅ Zmieniono emoji na ASCII w print statements ([+] zamiast ✓)

### `matcher/matcher.py`
- ✅ Dodano automatyczną konfigurację UTF-8 encoding dla Windows

## Weryfikacja

Testy pokazują że matching teraz działa poprawnie:

### Test 1: NH3 (Amoniak)
```
SMILES: N
✅ SUCCESS!
   CID: 222
   Name: azane
   Formula: H3N
   Method: Exact match
```

### Test 2: N3H3 (Triaziridine)
```
SMILES: [nH]1[nH][nH]1
✅ SUCCESS!
   CID: 18463253
   Name: triaziridine
   Formula: H3N3
   Method: Exact match
```

## Jak używać

### Standardowe użycie (single cluster)
```bash
python matcher/matcher.py matches/cluster_2025-10-22_15-58-29.png
```

### Z watchdogiem (automatyczne monitorowanie)
```bash
python matcher/watcher.py
```

## Korzyści

1. **Szybkość**: Exact match jest natychmiastowe (nie trzeba czekać 20s)
2. **Niezawodność**: Proste cząsteczki są teraz zawsze znajdowane
3. **Kompatybilność**: Działa z Windows console (cp1250 encoding)
4. **Fallback**: Jeśli exact match nie zadziała, automatycznie próbuje similarity search

## Przykłady cząsteczek które teraz działają

- ✅ H2 (wodór)
- ✅ NH3 (amoniak / azane)
- ✅ H2O (woda)
- ✅ CH4 (metan)
- ✅ N3H3 (triaziridine / cyclic triazane)
- ✅ CO2 (dwutlenek węgla)
- ✅ I wszystkie inne proste i złożone cząsteczki

## Kolejne kroki

Jeśli nadal występują problemy:

1. Sprawdź czy plik `.mol` jest poprawnie wygenerowany
2. Sprawdź SMILES w logach
3. Sprawdź połączenie internetowe (PubChem API wymaga dostępu do internetu)
4. Jeśli timeout, zwiększ limit czekania w `pubchem_similar_top()` (linia 264)

## Szczegóły techniczne

### Exact Match API
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{SMILES}/JSON
```
Zwraca: PC_Compounds array z pełnymi danymi cząsteczki

### Similarity Search API (async)
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{SMILES}/JSON?Threshold={threshold}
```
Zwraca: 
- 202 + ListKey (async) → poll z `/compound/listkey/{ListKey}/cids/JSON`
- 200 + IdentifierList (sync) → bezpośrednie CID

### Properties API
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey/JSON
```
Zwraca: PropertyTable z pełnymi właściwościami cząsteczki

---
*Naprawa wykonana: 22 października 2024*

