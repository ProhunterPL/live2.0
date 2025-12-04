truth-filter2.md

04.12.2025 

🎯 Cel TRUTH-FILTER 2.0

Celem nie jest „idealna chemia”, tylko:

Wyciąć rzeczy ewidentnie nierealistyczne przy Twoim modelu (Morse + LJ, brak QM).

Ostrożnie traktować wszystkie aromaty i mocno napięte pierścienie.

Jawnie oznaczać poziom zaufania do każdej cząsteczki: REJECT / FLAG / ACCEPT.

Dawać jasny ślad w logach, dlaczego coś zostało wycięte lub oflagowane.

🧱 Ogólny pipeline TRUTH-FILTER 2.0

Wejście:

lista wykrytych cząsteczek (SMILES / graf wiązań / skład),

metadane z symulacji (scenariusz, czas, energia, liczba wystąpień, match do PubChem).

Wyjście:
dla każdej cząsteczki:

{
  "smiles": "...",
  "formula": "C8H12N2O3",
  "mass": 184.0,
  "validity": "ACCEPT | FLAG | REJECT",
  "reasons": ["VALENCE_OK", "AROMATIC_UNSUPPORTED", "HIGH_STRAIN", ...],
  "confidence": 0.0–1.0
}


Kroki:

Valence & basic sanity check

Charge & connectivity

Ring strain & geometry heuristics

Aromaticity policy (NEW)

Model-compatibility (czy to w ogóle może wyjść z Morse + LJ)

Cross-check z bazami (PubChem/ChEBI)

Statystyka wystąpień / stabilność w czasie

Scoring + decyzja (ACCEPT/FLAG/REJECT)

1️⃣ Valence check (twardy filtr)

Zadanie dla agenta (implementer):

Użyj RDKit (jeśli dostępne) albo własne reguły walencyjne.

Dla każdego atomu policz liczbę wiązań (z uwzględnieniem podwójnych/potrójnych).

Reguły „hard fail” (→ REJECT):

C: valence not in {4} (po uwzględnieniu wiązań wielokrotnych).

N: >4 lub <2 (bez sensownych wyjątków typu N+
O).

O: >2.

H: >1.

Jakiekolwiek „floating” atomy bez połączeń (poza monoatomowym H2, O2 itd.).

Instrukcja dla agenta:

Jeśli valence check nie przechodzi → ustaw validity="REJECT", confidence=0.0, dodaj reason "VALENCE_ERROR" i nie analizuj dalej tej cząsteczki.

2️⃣ Charge & connectivity

Sprawdź, czy suma ładunków jest rozsądna (preferuj 0, ±1).

Sprawdź, czy graf jest spójny (jedna cząsteczka, nie „dwupak”).

Hard fail:

|total_charge| > 2

cząsteczka składa się z >1 niepołączonych komponentów.

→ REJECT z reason "UNPHYSICAL_CHARGE" / "DISCONNECTED_COMPONENTS".

3️⃣ Ring strain & geometry heuristics

Model jest 2D + klasyczne potencjały → nie chcemy ekstremalnie napiętych pierścieni.

Agent:

Znajdź wszystkie pierścienie (DFS / RDKit GetRingInfo()).

Dla każdego pierścienia oblicz:

liczbę atomów n,

liczbę heteroatomów (N, O),

obecność wiązań podwójnych.

Soft / hard rules:

Hard REJECT:

pierścień 3-członowy z >1 heteroatomem,

pierścień 4-członowy z >1 heteroatomem + wiązania podwójne.

FLAG (podejrzane, ale możliwe):

bicykle/cage (więcej niż jeden pierścień współdzielący atomy),

pierścienie 5–6-członowe z 2+ heteroatomami.

Dodaj reason "HIGH_STRAIN_RING".

4️⃣ Aromaticity policy (KLUCZOWE)

Założenie modelu:
Morse + LJ nie opisuje delokalizacji π, więc aromatyczność jako efekt kwantowy nie jest modelowana.
Dlatego wszystkie pierścienie aromatyczne traktujemy jako „chemically plausible, model-incompatible”.

Agent:

Włącz detekcję aromatyczności (RDKit, lub własne heurystyki: naprzemienne wiązania podwójne, planarne 5/6-członowe itp.).

Jeśli wykryto pierścień aromatyczny:

Ustaw:

validity = "FLAG"
reasons += ["AROMATIC_UNSUPPORTED_BY_MODEL"]
confidence *= 0.5


Dodaj pole:

"model_compatibility": "low"


UWAGA: Nie musisz tego wyrzucać → w paperze możesz uczciwie napisać:

„Aromatyczne cząsteczki są przedstawione jako topologicznie możliwe, ale nasz model nie zawiera pełnej stabilizacji aromatycznej.”

5️⃣ Model-compatibility score

Chcemy „zbić” score dla rzeczy, które wymagają chemii spoza Twojego modelu.

Heurystyki, które agent może zakodować:

+0.1 jeśli cząsteczka jest:

≤10 heavy atoms,

bez aromatów,

bez pierścieni o dużym naprężeniu.

–0.2 jeśli:

ma >1 pierścień,

ma 3+ heteroatomy w jednym pierścieniu,

ma wiązanie typu C≡C/C≡N w mocno rozgałęzionej strukturze.

–0.3 jeśli:

ma aromat (AROMATIC_UNSUPPORTED_BY_MODEL).

Startowy confidence: 0.8 po przejściu walencji.
Następnie modyfikujesz go heurystykami do zakresu [0,1].

6️⃣ Cross-check z bazami (PubChem / internal db)

Jeśli już masz matcher:

Jeśli jest match z PubChem / znaną bazą:

reasons += ["MATCH_KNOWN_CHEMISTRY"]

confidence += 0.1 (max 1.0)

Jeśli brak matchu:

reasons += ["NO_DB_MATCH"]

nie obniżaj mocno — to przecież „novel”.

To jest dodatkowa warstwa, nie hard-fail.

7️⃣ Statystyka wystąpień / stabilność w czasie

Agent powinien:

policzyć, w ilu różnych krokach i snapshotach dana cząsteczka się pojawia,

policzyć maksymalną liczebność w czasie.

Reguły:

Jeśli cząsteczka pojawia się tylko raz w jednym kroku i znika →
reasons += ["TRANSIENT_SINGLETON"], confidence *= 0.5.

Jeśli cząsteczka utrzymuje się przez >N kroków i/lub > X kopii:

reasons += ["PERSISTENT_SPECIES"], confidence += 0.1.

To jest ważne, żeby nie brać za „novel compound” pojedynczego glitcha w grafie wiązań.

8️⃣ Final decision logic

Na końcu agent robi mapowanie:

if "VALENCE_ERROR" in reasons or "UNPHYSICAL_CHARGE" in reasons or "DISCONNECTED_COMPONENTS" in reasons:
    validity = "REJECT"
elif confidence < 0.4:
    validity = "REJECT"
elif "AROMATIC_UNSUPPORTED_BY_MODEL" in reasons or "HIGH_STRAIN_RING" in reasons or "TRANSIENT_SINGLETON" in reasons:
    validity = "FLAG"
else:
    validity = "ACCEPT"


Interpretacja:

ACCEPT → można używać w analizie głównej (species counts, networks, etc.).

FLAG → można pokazywać w wynikach, ale:

w paperze opisujesz je jako „putative / tentative”,

nie opierasz głównych tez tylko na flagged molecules.

REJECT → ignorowane przy budowaniu sieci, cykli, statystyk.