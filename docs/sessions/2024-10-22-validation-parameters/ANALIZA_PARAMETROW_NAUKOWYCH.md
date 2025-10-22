# Analiza Parametrów Symulacji vs Literatura Naukowa

## Data: 2025-10-22

---

## 🔬 PARAMETRY OBECNE W KODZIE

### 1. Zasięg Wiązania (Binding Range)

**Kod: `backend/sim/core/binding.py:310`**
```python
if r <= PARTICLE_RADIUS_COMPILE * 2.0:  # Binding range
```

**Obecne:**
- `particle_radius = 0.5` Å (z frontendu)
- Binding range = 2 × 0.5 = **1.0 Å**

---

### 2. Parametry Lennard-Jones

**Z `data/physics_parameters_example.json`:**
```json
"C": {
  "epsilon": 0.105,  // kJ/mol
  "sigma": 3.431     // Å (UFF force field)
}
```

**W kodzie `backend/sim/core/spatial_hash.py:159-160`:**
```python
sigma = 1.0       # Å - ZANIŻONE!
epsilon = 0.5     # kJ/mol
```

---

### 3. Siła Wiązań Kowalencyjnych

**Kod: `backend/sim/core/binding.py:519`**
```python
1: {'k_spring': 10.0, 'rest_len': 0.8, 'strength': 20.0}  # covalent
```

---

### 4. Gęstość Cząsteczek

**Obecna konfiguracja frontend:**
```typescript
grid_height: 128 Å
grid_width: 128 Å  
max_particles: 5000
```

**Objętość:**
```
V = 128 × 128 × (zakładam głębokość 10 Å dla 2D) = 163,840 Å³
```

**Gęstość:**
```
ρ = 5000 / 163,840 = 0.0305 cząsteczek/Å³ = 30.5 cząsteczek/nm³
```

---

## 📚 LITERATURA NAUKOWA - STANDARDY

### 1. Zasięg Wiązań (Bond Formation Distance)

**Literatura:**
- **Van der Waals radii:**
  - C: 1.70 Å (Bondi, 1964)
  - N: 1.55 Å
  - O: 1.52 Å
  - H: 1.20 Å

- **Suma vdW radii** (maksymalny zasięg kontaktu):
  - C-C: 3.40 Å
  - C-N: 3.25 Å
  - C-O: 3.22 Å
  - N-H: 2.75 Å
  - O-H: 2.72 Å

- **Długości wiązań kowalencyjnych** (z physics_db):
  - C-C (single): 1.54 Å
  - C-O: 1.43 Å
  - O-H: 0.96 Å
  - N-H: 1.01 Å

**WNIOSEK:**
- ✅ Kowalencyjne bonds: 0.96-1.54 Å
- ✅ Van der Waals kontakt: 2.7-3.4 Å
- ❌ **Nasze 1.0 Å to ZA MAŁO!** (tylko ścisłe bonds, brak vdW)

**Standardy MD Simulations:**
- **GROMACS/NAMD:** cutoff dla niewiązanych = 2.5-3.0 × sigma
- Dla sigma = 3.431 Å → cutoff = 8.6-10.3 Å
- **Dla wiązań:** typowo 1.2-1.5 × długość równowagi
- Dla r_e = 1.54 Å → cutoff = 1.8-2.3 Å

---

### 2. Parametry Lennard-Jones

**Literatura (UFF Force Field, Rappé et al. 1992):**
```
C:  ε = 0.105 kJ/mol,  σ = 3.431 Å  ✅ W physics_db
N:  ε = 0.069 kJ/mol,  σ = 3.261 Å
O:  ε = 0.060 kJ/mol,  σ = 3.118 Å
H:  ε = 0.044 kJ/mol,  σ = 2.571 Å
```

**OPLS-AA Force Field (Jorgensen et al.):**
```
C(sp3): ε = 0.276 kJ/mol, σ = 3.50 Å
```

**AMBER Force Field:**
```
C:  ε = 0.36 kJ/mol,   σ = 3.40 Å
```

**NASZE w spatial_hash.py:**
```python
sigma = 1.0     # ❌ ZANIŻONE 3.4× (powinno 3.4 Å)
epsilon = 0.5   # ⚠️ Zawyżone 2-5× (powinno 0.1-0.3 kJ/mol)
```

---

### 3. Siła Wiązań Kowalencyjnych

**Literatura (z physics_db - Luo 2007):**
```
C-C single bond:
  D_e = 348 kJ/mol (dissociation energy)
  r_e = 1.54 Å
  k_spring = 2255 kJ/(mol·Å²)  [harmonic approximation]
```

**Harmonic spring constant:**
```
k = 2 × D_e × a²
gdzie a = 1.8 Å⁻¹ (Morse width)
k = 2 × 348 × 1.8² = 2254 kJ/(mol·Å²)
```

**NASZE:**
```python
k_spring = 10.0   # ❌ ZANIŻONE 225× (powinno ~2250!)
strength = 20.0   # ❌ ZANIŻONE 17× (powinno ~348)
```

---

### 4. Gęstość Cząsteczek

**Literatura - typowe gęstości w chemii prebiotycznej:**

**Miller-Urey (1953) - redukcyjna atmosfera:**
- CH₄: ~1% atmosfery = 0.01 bar
- NH₃: ~0.1% = 0.001 bar
- H₂O(para): zmienna
- Temperatura: ~298 K (25°C)

**Ideal Gas Law:**
```
n/V = P/(RT)
P = 0.01 bar = 1000 Pa
R = 8.314 J/(mol·K)
T = 298 K

n/V = 1000 / (8.314 × 298) = 0.403 mol/m³
```

**Przeliczenie na cząsteczki/Å³:**
```
0.403 mol/m³ × 6.022×10²³ / 10³⁰ Å³ = 2.43×10⁻⁴ cząsteczek/Å³
```

**Woda ciekła (dla porównania):**
```
ρ = 1 g/cm³
M = 18 g/mol
n = ρ/M = 0.0556 mol/cm³ = 55.6 mol/L

Cząsteczki/Å³ = 55.6 × 6.022×10²³ / 10²⁷ = 0.0335 cząsteczek/Å³
```

**NASZE:**
```
0.0305 cząsteczek/Å³ ≈ 91% gęstości wody ciekłej ✅ WYSOKIE!
```

**Typowe MD simulations (gas phase):**
- Gas phase: 10⁻⁴ - 10⁻³ cząsteczek/Å³
- Liquid phase: 0.03-0.05 cząsteczek/Å³
- **NASZA: 0.0305 = LIQUID PHASE DENSITY** ✅

---

## 🎯 WNIOSKI

### ✅ CO JEST POPRAWNE:

1. **Gęstość cząsteczek: 0.0305/Å³**
   - Odpowiada gęstości cieczy (jak woda)
   - WYSTARCZAJĄCO GĘSTE dla tworzenia wiązań!
   - 📊 Literatura: 0.03-0.05/Å³ dla cieczy

2. **Energia pulsów: 8.0**
   - Miller-Urey: 50-100 kJ/mol wyładowania
   - Nasze ~8 (skala bezwymiarowa) wydaje się OK

3. **Physics Database:**
   - Parametry z literatury (UFF, Luo 2007)
   - ✅ Poprawne wartości C-C bond: 348 kJ/mol, 1.54 Å

### ❌ CO JEST ZŁE:

1. **Zasięg wiązania: 1.0 Å → ZA MAŁO!**
   - Literatura: 2.7-3.4 Å dla vdW
   - Literatura: 8.6-10.3 Å dla LJ cutoff
   - **POWINNO BYĆ: 3.0-4.0 Å minimum**

2. **LJ Sigma w spatial_hash: 1.0 Å → ZANIŻONE 3.4×**
   - Literatura: 3.431 Å (UFF)
   - **POWINNO BYĆ: 3.4 Å**

3. **Siła wiązań: k=10 → ZANIŻONE 225×**
   - Literatura: 2255 kJ/(mol·Å²)
   - **POWINNO BYĆ: 500-2000** (nawet 1/4 byłoby lepsze)

4. **Próg probability: 0.6 → ZA WYSOKI**
   - Przy exp(-r/2) dla r=1.0 Å → exp(-0.5) = 0.606
   - Tylko cząsteczki w BEZPOŚREDNIM kontakcie tworzą bonds
   - **POWINNO BYĆ: 0.2-0.3**

---

## 💡 REKOMENDACJE

### OPCJA A: POPRAW PARAMETRY (ZALECANE!)

**Gęstość jest OK! Problem to parametry binding!**

**Zmiany w `backend/sim/core/binding.py`:**

1. **Linia 310 - zwiększ zasięg:**
```python
if r <= PARTICLE_RADIUS_COMPILE * 6.5:  # było: 2.0 → teraz 3.25 Å
```
**Uzasadnienie:** 3.25 Å = średnia suma vdW radii (C-N, C-O)

2. **Linia 315 - zmniejsz próg:**
```python
if binding_probability > 0.25:  # było: 0.6
```
**Uzasadnienie:** Pozwala na wiązania przy r=2-3 Å (realistyczne)

3. **Linia 329 - zmniejsz mass_ratio:**
```python
if mass_ratio > 0.4:  # było: 0.7
```
**Uzasadnienie:** 
- O-H: ratio = 1/16 = 0.063 → potrzebujemy niższego progu
- Ale ratio > 0.4 to sensowny kompromis (C-O = 0.75 ✅)

4. **Linia 517-522 - zwiększ siłę wiązań:**
```python
self.bond_type_params = {
    0: {'k_spring': 2.0, 'rest_len': 1.0, 'damping': 0.1, 'strength': 5.0},     # vdW
    1: {'k_spring': 500.0, 'rest_len': 0.8, 'damping': 0.2, 'strength': 100.0}, # covalent - ZWIĘKSZONE
    2: {'k_spring': 50.0, 'rest_len': 1.2, 'damping': 0.15, 'strength': 30.0},  # H-bond
    3: {'k_spring': 100.0, 'rest_len': 0.9, 'damping': 0.25, 'strength': 50.0}  # metallic
}
```
**Uzasadnienie:**
- k=500 to 1/4 literaturowego 2255 (kompromis dla stabilności)
- strength=100 to ~1/3 literaturowego 348 kJ/mol
- Nadal naukowe, ale bardziej stabilne numerycznie

5. **`backend/sim/core/spatial_hash.py:159` - popraw sigma:**
```python
sigma = 3.4  # było: 1.0 → teraz zgodne z UFF
```

---

### OPCJA B: ZWIĘKSZ LICZBĘ CZĄSTECZEK (NIE ZALECANE)

**Dlaczego NIE:**
- ✅ Mamy już 0.0305/Å³ = 91% gęstości wody!
- ✅ To LIQUID PHASE density - wystarczająco gęste!
- ❌ Więcej cząsteczek = wolniejsza symulacja
- ❌ Więcej cząsteczek = większy memory footprint
- ❌ NIE ROZWIĄŻE problemu za małego zasięgu wiązania

**Jeśli chcesz zwiększyć:**
```typescript
max_particles: 8000  // było: 5000 → +60%
```
- Nowa gęstość: 0.049/Å³ = 146% wody (BARDZO GĘSTE!)
- ⚠️ Może być za gęste i destabilizować symulację

---

## 🎯 OSTATECZNA REKOMENDACJA

### ✅ WYKONAJ: OPCJA A - Popraw Parametry

**Dlaczego:**
1. ✅ Gęstość jest już OPTYMALNA (0.0305/Å³)
2. ✅ Parametry są NIENAUKOWE (binding range 3× za mały!)
3. ✅ Literaturowe parametry: sigma=3.4 Å, cutoff=3.0 Å
4. ✅ Proste zmiany w 1 pliku (binding.py)
5. ✅ Zachowuje naukowy realizm
6. ✅ Nie spowalnia symulacji

**Oczekiwane rezultaty:**
- ✅ Klastry 8-20 atomów (glikol, formamid, urea)
- ✅ Stabilne wiązania C-O, C-N, O-H
- ✅ Zgodność z Miller-Urey (1953) produktami
- ✅ Symulacja stabilna >200k steps

**NIE zwiększaj liczby cząsteczek** - gęstość jest już bardzo wysoka!

---

## 📋 PODSUMOWANIE ZMIAN

```python
# backend/sim/core/binding.py

# Linia 310:
if r <= PARTICLE_RADIUS_COMPILE * 6.5:  # 3.25 Å zasięg

# Linia 315:
if binding_probability > 0.25:  # Niższy próg

# Linia 329:
if mass_ratio > 0.4:  # Pozwala C-O bonds

# Linia 519:
1: {'k_spring': 500.0, 'rest_len': 0.8, 'strength': 100.0},  # Silniejsze bonds
```

```python
# backend/sim/core/spatial_hash.py

# Linia 159:
sigma = 3.4  # UFF standard
```

---

## 📚 CYTOWANA LITERATURA

1. **Bondi, A. (1964)** - Van der Waals Volumes and Radii
   - J. Phys. Chem., 68(3), 441-451

2. **Rappé et al. (1992)** - UFF Force Field
   - J. Am. Chem. Soc., 114(25), 10024-10035

3. **Luo, Y.-R. (2007)** - Comprehensive Handbook of Chemical Bond Energies
   - CRC Press (już w physics_db!)

4. **Miller & Urey (1953)** - Organic Compound Synthesis on Primitive Earth
   - Science, 117(3046), 528-529

5. **GROMACS Manual (2023)** - MD Simulation Best Practices
   - www.gromacs.org

---

## ❓ DECYZJA

**Czy wprowadzić OPCJĘ A (popraw parametry)?**
- TAK → wprowadzę zmiany w backend/sim/core/binding.py
- NIE → wyjaśnij dlaczego, zaproponuję alternatywę

