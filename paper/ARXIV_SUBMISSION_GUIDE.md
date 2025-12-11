---
date: 2025-12-04
label: guide
---

# arXiv Submission Guide - Live 2.0 Manuscript

**Status**: ✅ Manuscript prepared for arXiv submission  
**Journal Submission**: Origins of Life and Evolution of Biospheres (2025-12-04)  
**arXiv Policy**: ✅ Allowed (preprint before peer review)

---

## 📋 Podsumowanie

Ten dokument zawiera instrukcje krok po kroku dla submissionu manuskryptu na arXiv. Manuskrypt został już złożony do czasopisma (Origins of Life and Evolution of Biospheres), więc preprint na arXiv jest zgodny z polityką wydawcy.

---

## ✅ Przygotowanie Plików

### 1. Wersja Manuskryptu dla arXiv

**Plik**: `paper/manuscript_arxiv.tex`

**Zmiany w stosunku do wersji dla czasopisma**:
- ✅ Usunięto `\usepackage{lineno}` i `\linenumbers` (line numbers nie są potrzebne na arXiv)
- ✅ Dodano notkę o submission do czasopisma (w komentarzu na początku)
- ✅ Wszystkie pakiety LaTeX są kompatybilne z arXiv

**Sprawdź przed submissionem**:
- [ ] Wszystkie figury są dostępne w folderze `figures/`
- [ ] Wszystkie tabele są dostępne w folderze `tables/`
- [ ] Bibliografia (`references.bib`) jest kompletna
- [ ] Kompilacja LaTeX działa lokalnie

### 2. Struktura Plików dla arXiv

```
arxiv_submission/
├── manuscript_arxiv.tex          # Główny plik LaTeX
├── references.bib                 # Bibliografia
├── figures/                       # Wszystkie figury
│   ├── figure1_thermodynamic_validation.png
│   ├── figure2_benchmark_validation.png
│   ├── figure2_formose_validation.png
│   ├── figure3_molecular_diversity.png
│   ├── figure4_reaction_networks.png
│   ├── figure5_autocatalytic_cycles.png
│   ├── figure6_novel_molecules.png
│   ├── figure6b_novel_structures.png
│   └── molecular_structures_panel.png
└── tables/                        # Wszystkie tabele
    ├── table5_hub_molecules.tex
    └── table6_novel_molecules.tex
```

**Uwaga**: arXiv automatycznie kompiluje LaTeX, więc wszystkie pliki muszą być w jednym katalogu lub używać względnych ścieżek.

---

## 🔧 Krok 1: Rejestracja na arXiv

### 1.1 Utworzenie Konta

1. **Przejdź do**: https://arxiv.org/register
2. **Wypełnij formularz**:
   - Email (używaj profesjonalnego adresu)
   - Imię i nazwisko
   - Afiliacja: "Independent Researcher, Pruszcz Gdański, Poland"
   - Hasło (silne, min. 8 znaków)
3. **Potwierdź email** (sprawdź skrzynkę i kliknij link)

### 1.2 Endorsement (Wymagane dla pierwszego submissionu)

**Dla kategorii `q-bio.BM` (Quantitative Biology - Biomolecules)**:
- Potrzebujesz endorsement od aktywnego użytkownika arXiv w tej kategorii
- **Alternatywnie**: Użyj kategorii `physics.chem-ph` (Physics - Chemical Physics), która nie wymaga endorsementu dla pierwszego submissionu

**Rekomendacja**: Użyj `physics.chem-ph` dla pierwszego submissionu (łatwiejsze), możesz później zmienić kategorię.

---

## 📤 Krok 2: Submission na arXiv

### 2.1 Przygotowanie Plików

1. **Utwórz folder** `arxiv_submission/` w `paper/`
2. **Skopiuj pliki**:
   ```bash
   cd paper
   mkdir -p arxiv_submission/figures arxiv_submission/tables
   cp manuscript_arxiv.tex arxiv_submission/
   cp references.bib arxiv_submission/
   cp figures/*.png arxiv_submission/figures/
   cp tables/*.tex arxiv_submission/tables/
   ```

3. **Sprawdź kompilację lokalnie**:
   ```bash
   cd arxiv_submission
   pdflatex manuscript_arxiv.tex
   bibtex manuscript_arxiv
   pdflatex manuscript_arxiv.tex
   pdflatex manuscript_arxiv.tex
   ```
   
   **Upewnij się, że**:
   - PDF kompiluje się bez błędów
   - Wszystkie figury są widoczne
   - Wszystkie tabele są widoczne
   - Wszystkie referencje są poprawnie wyświetlone

### 2.2 Utworzenie Archiwum

**Opcja A: ZIP (Rekomendowane)**
```bash
cd paper/arxiv_submission
zip -r ../arxiv_submission.zip .
```

**Opcja B: TAR.GZ**
```bash
cd paper/arxiv_submission
tar -czf ../arxiv_submission.tar.gz .
```

### 2.3 Submission przez Web Interface

1. **Zaloguj się**: https://arxiv.org/login
2. **Kliknij**: "Submit to arXiv"
3. **Wybierz**: "New submission"
4. **Wypełnij formularz**:

   **Primary Classification**:
   - `physics.chem-ph` (Physics - Chemical Physics) - **REKOMENDOWANE** (nie wymaga endorsementu)
   - LUB `q-bio.BM` (Quantitative Biology - Biomolecules) - wymaga endorsementu

   **Title**:
   ```
   Emergent Molecular Complexity in Prebiotic Chemistry Simulations: A Physics-Based Approach
   ```

   **Authors**:
   ```
   Michał Klawikowski (Independent Researcher, Pruszcz Gdański, Poland)
   ```

   **Abstract**:
   ```
   [Wklej abstract z manuskryptu - 250 słów]
   ```

   **Comments** (Opcjonalnie):
   ```
   This manuscript has been submitted to Origins of Life and Evolution of Biospheres for peer review.
   ```

   **Keywords**:
   ```
   prebiotic chemistry, origin of life, molecular dynamics, autocatalysis, emergent complexity
   ```

5. **Upload plików**:
   - Wybierz archiwum ZIP lub TAR.GZ
   - LUB uploaduj pliki pojedynczo (mniej wygodne)

6. **Sprawdź preview**:
   - arXiv automatycznie skompiluje LaTeX
   - Sprawdź PDF preview
   - Jeśli są błędy, popraw i upload ponownie

7. **Submit**:
   - Kliknij "Submit"
   - Otrzymasz email z potwierdzeniem
   - Submission ID będzie w emailu

---

## ⏱️ Krok 3: Proces Review na arXiv

### 3.1 Timeline

- **Automatyczna kompilacja**: 5-30 minut
- **Moderation** (jeśli potrzebne): 1-2 dni
- **Publikacja**: Zwykle w ciągu 24-48 godzin

### 3.2 Status Submissionu

Sprawdzaj status na: https://arxiv.org/user

**Możliwe statusy**:
- `Submitted` - Oczekuje na przetworzenie
- `Processing` - Kompilacja LaTeX
- `On hold` - Wymaga moderacji (sprawdź email)
- `Announced` - Opublikowane! ✅

### 3.3 Jeśli Submission jest "On Hold"

**Możliwe przyczyny**:
- Błędy kompilacji LaTeX
- Problemy z formatowaniem
- Potrzeba endorsementu (dla niektórych kategorii)

**Co zrobić**:
1. Sprawdź email od arXiv (będzie szczegółowy opis problemu)
2. Popraw błędy
3. Resubmit (użyj tego samego submission ID)

---

## 📝 Krok 4: Po Publikacji

### 4.1 Aktualizacja Dokumentacji

**Zaktualizuj** `paper/SUBMISSION_LOG.md`:
```markdown
## arXiv Submission

- **Date**: [DATA]
- **arXiv ID**: [np. 2412.XXXXX]
- **URL**: https://arxiv.org/abs/[ID]
- **Status**: ✅ Published
```

### 4.2 Aktualizacja Manuskryptu (Opcjonalnie)

Jeśli chcesz dodać link do arXiv w manuskrypcie (po publikacji):
- Dodaj w sekcji "Data and Code Availability"
- Format: `Preprint available at: https://arxiv.org/abs/[ID]`

---

## ✅ Checklist Przed Submissionem

### Przygotowanie Plików
- [ ] `manuscript_arxiv.tex` jest gotowy (bez line numbers)
- [ ] Wszystkie figury są w folderze `figures/`
- [ ] Wszystkie tabele są w folderze `tables/`
- [ ] `references.bib` jest kompletny
- [ ] Kompilacja lokalna działa bez błędów

### Rejestracja
- [ ] Konto na arXiv utworzone
- [ ] Email potwierdzony
- [ ] Endorsement uzyskany (jeśli potrzebny dla wybranej kategorii)

### Submission
- [ ] Archiwum ZIP/TAR.GZ utworzone
- [ ] Formularz submissionu wypełniony
- [ ] PDF preview sprawdzony
- [ ] Submission wysłany

### Po Submissionie
- [ ] Email z potwierdzeniem otrzymany
- [ ] Status submissionu monitorowany
- [ ] Dokumentacja zaktualizowana po publikacji

---

## 🎯 Kategorie arXiv

### Rekomendowane

**`physics.chem-ph` (Physics - Chemical Physics)**
- ✅ Nie wymaga endorsementu dla pierwszego submissionu
- ✅ Dobrze pasuje do tematyki (computational chemistry, molecular dynamics)
- ✅ Szybka publikacja (zwykle <24h)

**`q-bio.BM` (Quantitative Biology - Biomolecules)**
- ⚠️ Wymaga endorsementu dla pierwszego submissionu
- ✅ Bardziej pasuje do tematyki (prebiotic chemistry, origin of life)
- ✅ Można użyć po uzyskaniu endorsementu

### Alternatywne

- `q-bio.QM` (Quantitative Biology - Quantitative Methods)
- `physics.bio-ph` (Physics - Biological Physics)

---

## 📞 Pomoc i Zasoby

### arXiv Support
- **Email**: help@arxiv.org
- **FAQ**: https://arxiv.org/help
- **Submission FAQ**: https://arxiv.org/help/submit

### Endorsement
- **Jak uzyskać endorsement**: https://arxiv.org/help/endorsement
- **Lista endorserów**: https://arxiv.org/help/endorsement#list

### LaTeX Help
- **arXiv LaTeX Guide**: https://arxiv.org/help/submit_tex
- **Common Issues**: https://arxiv.org/help/faq/mistakes

---

## 🔄 Aktualizacje i Wersje

**Wersja 1.0** (2025-12-04):
- ✅ Instrukcje dla pierwszego submissionu
- ✅ Checklist i timeline
- ✅ Informacje o kategoriach

**Uwagi**:
- Po pierwszym submissionie, kolejne są szybsze (nie wymagają endorsementu)
- Można updateować preprint po publikacji w czasopiśmie (z odpowiednią notką)

---

**Ostatnia aktualizacja**: 2025-12-04  
**Status**: ✅ Gotowe do użycia

