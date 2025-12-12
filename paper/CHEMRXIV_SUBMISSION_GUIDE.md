---
date: 2025-12-04
label: guide
---

# chemRxiv Submission Guide - Live 2.0 Manuscript

**Status**: ✅ Manuscript prepared for chemRxiv submission  
**Journal Submission**: Origins of Life and Evolution of Biospheres (2025-12-04)  
**chemRxiv Policy**: ✅ Allowed (preprint before peer review)  
**Advantage**: ✅ No endorsement required (unlike arXiv)

---

## 📋 Podsumowanie

Ten dokument zawiera instrukcje krok po kroku dla submissionu manuskryptu na chemRxiv. chemRxiv to preprint server specjalizujący się w chemii, który:
- ✅ **Nie wymaga endorsementu** (w przeciwieństwie do arXiv)
- ✅ Przyjmuje PDF (nie wymaga LaTeX source)
- ✅ Szybka publikacja (zwykle <24h)
- ✅ Idealnie pasuje do tematyki (prebiotic chemistry, computational chemistry)

---

## ✅ Przygotowanie Plików

### 1. Wersja PDF dla chemRxiv

**Rekomendacja**: Użyj tej samej wersji co dla czasopisma, ale **bez line numbers**.

**Opcje**:
- **Opcja A**: Skompiluj `manuscript_arxiv.tex` (już bez line numbers) → PDF
- **Opcja B**: Skompiluj `manuscript_draft.tex` z wyłączonymi line numbers → PDF

**Sprawdź przed submissionem**:
- [ ] PDF ma wszystkie figury (9 figur)
- [ ] PDF ma wszystkie tabele (2 tabele)
- [ ] Wszystkie referencje są poprawnie wyświetlone
- [ ] Brak line numbers (linie numerowane nie są potrzebne)
- [ ] Rozmiar PDF < 10 MB (chemRxiv limit)

### 2. Wymagane Pliki

**Minimalne wymagania**:
- ✅ **Main PDF**: Manuscript jako pojedynczy plik PDF
- ✅ **Title**: Tytuł manuskryptu
- ✅ **Abstract**: Abstract (250 słów)
- ✅ **Keywords**: Słowa kluczowe
- ✅ **Author information**: Imię, nazwisko, afiliacja

**Opcjonalne**:
- Supplementary materials (można dodać później)
- Figure files (jeśli chcesz osobno, ale nie jest wymagane)

---

## 🔧 Krok 1: Rejestracja na chemRxiv

### 1.1 Utworzenie Konta

1. **Przejdź do**: https://chemrxiv.org/engage/chemrxiv/account/register
2. **Wypełnij formularz**:
   - Email (używaj profesjonalnego adresu)
   - Imię i nazwisko: Michał Klawikowski
   - Afiliacja: Independent Researcher, Pruszcz Gdański, Poland
   - Hasło (silne, min. 8 znaków)
   - ORCID (opcjonalnie, ale zalecane)
3. **Potwierdź email** (sprawdź skrzynkę i kliknij link)

### 1.2 Weryfikacja Konta

- chemRxiv może wymagać weryfikacji email
- Sprawdź spam folder, jeśli nie otrzymasz emaila
- Po weryfikacji możesz od razu submitować

**Uwaga**: chemRxiv **NIE wymaga endorsementu** - to główna zaleta w stosunku do arXiv!

---

## 📤 Krok 2: Submission na chemRxiv

### 2.1 Przygotowanie PDF

1. **Skompiluj PDF** (jeśli jeszcze nie masz):
   ```bash
   cd paper
   # Użyj manuscript_arxiv.tex (bez line numbers)
   pdflatex manuscript_arxiv.tex
   bibtex manuscript_arxiv
   pdflatex manuscript_arxiv.tex
   pdflatex manuscript_arxiv.tex
   ```

2. **Sprawdź PDF**:
   - Otwórz `manuscript_arxiv.pdf`
   - Upewnij się, że wszystkie figury są widoczne
   - Upewnij się, że wszystkie tabele są widoczne
   - Sprawdź, że nie ma line numbers

3. **Sprawdź rozmiar**:
   - PDF powinien być < 10 MB
   - Jeśli większy, zoptymalizuj figury (zmniejsz rozdzielczość do 300 DPI)

### 2.2 Submission przez Web Interface

1. **Zaloguj się**: https://chemrxiv.org/engage/chemrxiv/login
2. **Kliknij**: "Submit a Preprint" lub "New Submission"
3. **Wypełnij formularz**:

   **Title**:
   ```
   Emergent Molecular Complexity in Prebiotic Chemistry Simulations Using a Physics-Based Approach
   ```
   *(Użyj zaktualizowanego tytułu z response_to_springer_title_data.txt)*

   **Authors**:
   ```
   Michał Klawikowski
   Independent Researcher, Pruszcz Gdański, Poland
   Email: klawikowski@klawikowski.pl
   ```

   **Abstract**:
   ```
   [Wklej abstract z manuskryptu - 250 słów]
   ```

   **Keywords**:
   ```
   prebiotic chemistry, origin of life, molecular dynamics, autocatalysis, emergent complexity
   ```

   **Subject Area** (wybierz najbardziej pasujące):
   - ✅ **Computational Chemistry** (najlepsze dopasowanie)
   - ✅ **Physical Chemistry**
   - ✅ **Biochemistry**
   - ✅ **Theoretical Chemistry**

   **Comments** (Opcjonalnie):
   ```
   This manuscript has been submitted to Origins of Life and Evolution of Biospheres for peer review.
   ```

4. **Upload PDF**:
   - Kliknij "Upload Manuscript"
   - Wybierz plik PDF (`manuscript_arxiv.pdf`)
   - Poczekaj na upload (może zająć kilka minut dla dużego pliku)

5. **Sprawdź preview**:
   - chemRxiv wyświetli preview PDF
   - Sprawdź, czy wszystko wygląda dobrze
   - Jeśli są problemy, popraw PDF i upload ponownie

6. **Submit**:
   - Kliknij "Submit" lub "Publish"
   - Otrzymasz email z potwierdzeniem
   - Submission ID będzie w emailu

---

## ⏱️ Krok 3: Proces Review na chemRxiv

### 3.1 Timeline

- **Automatyczna weryfikacja**: 1-4 godziny
- **Publikacja**: Zwykle w ciągu 24 godzin (często <12h)

**Uwaga**: chemRxiv jest szybszy niż arXiv w publikacji!

### 3.2 Status Submissionu

Sprawdzaj status na: https://chemrxiv.org/engage/chemrxiv/user-dashboard

**Możliwe statusy**:
- `Submitted` - Oczekuje na przetworzenie
- `Under Review` - Weryfikacja techniczna
- `Published` - Opublikowane! ✅
- `Revision Required` - Wymaga poprawek (sprawdź email)

### 3.3 Jeśli Submission Wymaga Poprawek

**Możliwe przyczyny**:
- Problemy z formatowaniem PDF
- Brakujące informacje (abstract, keywords)
- Problemy z metadanymi

**Co zrobić**:
1. Sprawdź email od chemRxiv (będzie szczegółowy opis problemu)
2. Popraw błędy
3. Resubmit (użyj tego samego submission ID)

---

## 📝 Krok 4: Po Publikacji

### 4.1 Aktualizacja Dokumentacji

**Zaktualizuj** `paper/SUBMISSION_LOG.md`:
```markdown
## chemRxiv Submission

- **Date**: [DATA]
- **chemRxiv DOI**: [np. 10.26434/chemrxiv-XXXXX]
- **URL**: https://chemrxiv.org/engage/chemrxiv/article-details/[ID]
- **Status**: ✅ Published
```

### 4.2 Aktualizacja Manuskryptu (Opcjonalnie)

Jeśli chcesz dodać link do chemRxiv w manuskrypcie (po publikacji):
- Dodaj w sekcji "Data and Code Availability"
- Format: `Preprint available at: https://chemrxiv.org/engage/chemrxiv/article-details/[ID]`

### 4.3 Integracja z Journal Submission

**Ważne**: chemRxiv pozwala na preprint przed recenzją, więc:
- ✅ Możesz submitować do chemRxiv nawet po submission do czasopisma
- ✅ Po publikacji w czasopiśmie, możesz zaktualizować chemRxiv z informacją o publikacji
- ✅ chemRxiv automatycznie linkuje do finalnej publikacji (jeśli podasz DOI)

---

## ✅ Checklist Przed Submissionem

### Przygotowanie Plików
- [ ] PDF skompilowany bez line numbers
- [ ] Wszystkie figury są widoczne w PDF
- [ ] Wszystkie tabele są widoczne w PDF
- [ ] Wszystkie referencje są poprawnie wyświetlone
- [ ] Rozmiar PDF < 10 MB
- [ ] Tytuł zaktualizowany (bez dwukropka)

### Rejestracja
- [ ] Konto na chemRxiv utworzone
- [ ] Email potwierdzony
- [ ] ORCID dodany (opcjonalnie, ale zalecane)

### Submission
- [ ] Formularz submissionu wypełniony
- [ ] PDF uploadowany
- [ ] Preview sprawdzony
- [ ] Submission wysłany

### Po Submissionie
- [ ] Email z potwierdzeniem otrzymany
- [ ] Status submissionu monitorowany
- [ ] Dokumentacja zaktualizowana po publikacji

---

## 🎯 Zalety chemRxiv vs arXiv

### chemRxiv ✅
- ✅ **Nie wymaga endorsementu** (główna zaleta!)
- ✅ Szybsza publikacja (<24h)
- ✅ Przyjmuje PDF (nie wymaga LaTeX source)
- ✅ Specjalizuje się w chemii (lepsze dopasowanie)
- ✅ Integracja z czasopismami chemicznymi
- ✅ Automatyczne linkowanie do finalnej publikacji

### arXiv
- ⚠️ Wymaga endorsementu dla niektórych kategorii
- ⚠️ Wymaga LaTeX source (kompilacja na serwerze)
- ⚠️ Mniej specjalistyczny dla chemii

**Rekomendacja**: chemRxiv jest lepszym wyborem dla tego manuskryptu!

---

## 📞 Pomoc i Zasoby

### chemRxiv Support
- **Email**: support@chemrxiv.org
- **FAQ**: https://chemrxiv.org/engage/chemrxiv/help
- **Submission Guide**: https://chemrxiv.org/engage/chemrxiv/help/submission-guide

### ORCID
- **Utworzenie konta**: https://orcid.org/register
- **Dlaczego warto**: Ułatwia śledzenie publikacji i cytowań

---

## 🔄 Aktualizacje i Wersje

**Wersja 1.0** (2025-12-04):
- ✅ Instrukcje dla submissionu na chemRxiv
- ✅ Checklist i timeline
- ✅ Porównanie z arXiv
- ✅ Informacje o integracji z czasopismami

**Uwagi**:
- chemRxiv pozwala na update preprintu po publikacji w czasopiśmie
- Można dodać DOI finalnej publikacji po akceptacji

---

**Ostatnia aktualizacja**: 2025-12-04  
**Status**: ✅ Gotowe do użycia  
**Zmiana z arXiv**: ✅ Ze względu na brak wymogu endorsementu

