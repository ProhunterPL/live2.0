---
date: 2025-12-04
label: plan
---

# Paper 2: "Autocatalytic Networks in Prebiotic Chemistry" - Outline

**Status**: 📋 Planning  
**Target Journal**: Origins of Life and Evolution of Biospheres (lub JCTC)  
**Timeline**: Luty-kwiecień 2025 (analysis + writing)  
**Target Submission**: Maj 2025

---

## 📋 Podsumowanie

Paper 2 skupia się na szczegółowej analizie autocatalytic networks wykrytych w Phase 2B. Podczas gdy Paper 1 pokazuje ogólną złożoność molekularną, Paper 2 zagłębia się w mechanizmy autocatalysis i topologię sieci reakcji.

**Key Message**: Autocatalytic networks są kluczowym mechanizmem emergencji złożoności w prebiotic chemistry, a ich topologia różni się znacząco między scenariuszami.

---

## 📊 Available Data

### Phase 2B Results
- ✅ **43 runs completed**: 18 Miller-Urey, 17 Hydrothermal, 8 Formamide
- ✅ **769,315 autocatalytic cycles** detected across all scenarios
- ✅ **Reaction networks** analyzed (molecules as nodes, reactions as edges)
- ✅ **All snapshots** saved (every 50K steps)
- ✅ **All checkpoints** saved (every 100K steps)

### Analysis Needed
- ⏳ **Detailed cycle classification**: Direct, indirect, hypercycles
- ⏳ **Amplification factors**: Quantitative analysis
- ⏳ **Network topology**: Degree distribution, hubs, centrality, motifs
- ⏳ **Scenario comparison**: Statistical analysis of differences
- ⏳ **Temporal evolution**: How networks evolve over time

---

## 📝 Manuscript Structure

### Title
"Autocatalytic Networks in Prebiotic Chemistry: Topology, Mechanisms, and Scenario-Dependent Emergence"

**Alternative**: "Emergent Autocatalytic Networks in Simulated Prebiotic Environments: A Network Topology Analysis"

---

### Abstract (~250 words)

**Background** (2-3 sentences):
- Autocatalysis jako kluczowy mechanizm self-organization w prebiotic chemistry
- Znaczenie network topology dla zrozumienia emergencji złożoności
- Brak systematycznych badań topologii autocatalytic networks w różnych scenariuszach

**Methods** (3-4 sentences):
- Analiza 43 symulacji Phase 2B (Miller-Urey, Hydrothermal, Formamide)
- Detekcja autocatalytic cycles (direct, indirect, hypercycles)
- Network topology analysis (degree distribution, hubs, centrality, motifs)
- Statistical comparison między scenariuszami

**Results** (4-5 sentences):
- [X] autocatalytic cycles detected (769,315 total)
- [X]% direct autocatalysis, [X]% indirect cycles, [X]% hypercycles
- Amplification factors: [X] (Miller-Urey), [X] (Hydrothermal), [X] (Formamide)
- Network topology różni się znacząco między scenariuszami (p < 0.01)
- Hub molecules identified: [top 5 per scenario]

**Significance** (2-3 sentences):
- Autocatalytic networks są powszechne i różnorodne w prebiotic chemistry
- Topologia sieci determinuje dynamikę emergencji złożoności
- Scenariusze różnią się mechanizmami autocatalysis (testable predictions)

---

### 1. Introduction (~1500 words)

#### 1.1 Autocatalysis in Prebiotic Chemistry (400 words)
- **Autocatalysis definition**: A + B → 2A (direct) vs A → B → C → A (indirect)
- **Historical context**: Breslow (1959) formose reaction, Eigen (1971) hypercycles
- **Modern relevance**: RAF sets, autocatalytic sets, chemical evolution
- **Key question**: Jak często i w jakich formach występuje autocatalysis w prebiotic chemistry?

#### 1.2 Network Topology and Chemical Evolution (400 words)
- **Network science**: Graph theory applied to chemical reaction networks
- **Topology metrics**: Degree distribution, hubs, centrality, motifs
- **Biological networks**: Scale-free, small-world properties
- **Prebiotic networks**: Czy mają podobne właściwości?
- **Key question**: Czy topologia sieci determinuje dynamikę emergencji?

#### 1.3 Scenario-Dependent Chemistry (400 words)
- **Miller-Urey**: Reducing atmosphere, electrical discharge
- **Hydrothermal vents**: Alkaline, high temperature, mineral surfaces
- **Formamide**: UV radiation, versatile precursor
- **Hypothesis**: Różne scenariusze → różne network topologies → różne mechanizmy autocatalysis

#### 1.4 Study Overview (300 words)
- **Our approach**: Physics-based simulation, emergent reactions
- **Data**: 43 runs, 769,315 cycles, comprehensive network analysis
- **Key questions**:
  1. Jakie typy autocatalytic cycles występują?
  2. Jak topologia sieci różni się między scenariuszami?
  3. Jakie są amplification factors?
  4. Czy hub molecules są konserwowane między scenariuszami?

---

### 2. Methods (~1800 words)

#### 2.1 Simulation Framework (400 words)
- **Brief summary**: Physics-based particle simulation (refer to Paper 1)
- **Scenarios**: Miller-Urey, Hydrothermal, Formamide (compositions, conditions)
- **Runs**: 18 Miller-Urey, 17 Hydrothermal, 8 Formamide (500K steps each)
- **Reference**: Paper 1 for detailed methods

#### 2.2 Autocatalytic Cycle Detection (500 words)
- **Algorithm**: DFS-based cycle detection in reaction graph
- **Direct autocatalysis**: A + B → 2A (detection method)
- **Indirect cycles**: A → B → C → A (DFS algorithm)
- **Hypercycles**: Mutual catalysis (A catalyzes B, B catalyzes A)
- **Validation**: Manual verification of top 100 cycles
- **Limitations**: Cycle detection in large networks (computational complexity)

#### 2.3 Network Topology Analysis (500 words)
- **Graph construction**: Molecules as nodes, reactions as edges
- **Degree distribution**: In-degree, out-degree, total degree
- **Hub identification**: Top 10% by degree centrality
- **Centrality metrics**: Betweenness, closeness, eigenvector centrality
- **Network motifs**: 3-node and 4-node motifs (Mfinder algorithm)
- **Statistical analysis**: Comparison between scenarios (ANOVA, t-tests)

#### 2.4 Amplification Factor Calculation (400 words)
- **Definition**: Rate increase per cycle iteration
- **Method**: Temporal analysis of cycle participants
- **Calculation**: d[A]/dt before vs after cycle activation
- **Validation**: Comparison with theoretical predictions

---

### 3. Results (~2000 words)

#### 3.1 Autocatalytic Cycle Classification (500 words)
- **Total cycles**: 769,315 across all scenarios
- **Direct autocatalysis**: [X]% ([X] cycles)
  - Examples: [top 5 direct cycles]
  - Amplification factors: [range, mean, std]
- **Indirect cycles**: [X]% ([X] cycles)
  - Examples: [top 5 indirect cycles]
  - Cycle lengths: [distribution]
- **Hypercycles**: [X]% ([X] cycles)
  - Examples: [top 3 hypercycles]
  - Mutual catalysis patterns
- **Figure 1**: Cycle type distribution (pie chart + examples)

#### 3.2 Amplification Factors (400 words)
- **Miller-Urey**: Mean [X], range [Y-Z], std [W]
- **Hydrothermal**: Mean [X], range [Y-Z], std [W]
- **Formamide**: Mean [X], range [Y-Z], std [W]
- **Statistical comparison**: ANOVA results (p < 0.01)
- **Figure 2**: Amplification factors by scenario (violin plots)

#### 3.3 Network Topology Analysis (600 words)
- **Degree distribution**: Power-law? Scale-free? (fit to power-law)
- **Hub molecules**: Top 10 per scenario
  - Miller-Urey: [list]
  - Hydrothermal: [list]
  - Formamide: [list]
  - Overlap analysis: [X]% shared hubs
- **Centrality metrics**: 
  - Betweenness: [top 5 per scenario]
  - Closeness: [top 5 per scenario]
  - Eigenvector: [top 5 per scenario]
- **Network motifs**: 
  - 3-node motifs: [frequency per scenario]
  - 4-node motifs: [frequency per scenario]
  - Statistical significance: [Z-scores]
- **Figure 3**: Network topology comparison (degree distribution, hubs, motifs)

#### 3.4 Scenario Comparison (500 words)
- **Cycle frequency**: Cycles per 1000 reactions (by scenario)
- **Cycle types**: Distribution differences (chi-square test)
- **Amplification factors**: Statistical comparison (ANOVA)
- **Network topology**: Topology metrics comparison (t-tests)
- **Hub molecules**: Overlap and uniqueness analysis
- **Figure 4**: Scenario comparison (multi-panel figure)
- **Table 1**: Summary statistics by scenario

---

### 4. Discussion (~1500 words)

#### 4.1 Autocatalysis as Emergent Mechanism (400 words)
- **Prevalence**: Autocatalysis jest powszechne w prebiotic chemistry
- **Diversity**: Różne typy cykli (direct, indirect, hypercycles)
- **Amplification**: Mechanizmy wzmacniania reakcji
- **Comparison with literature**: Eigen hypercycles, RAF sets, autocatalytic sets

#### 4.2 Network Topology and Chemical Evolution (400 words)
- **Scale-free properties**: Czy sieci są scale-free?
- **Hub molecules**: Rola hub molecules w network stability
- **Motifs**: Znaczenie network motifs dla funkcjonalności
- **Biological relevance**: Porównanie z biological networks

#### 4.3 Scenario-Dependent Mechanisms (400 words)
- **Miller-Urey**: [charakterystyka network topology]
- **Hydrothermal**: [charakterystyka network topology]
- **Formamide**: [charakterystyka network topology]
- **Implications**: Różne scenariusze → różne mechanizmy emergencji
- **Testable predictions**: Eksperymentalne weryfikacje

#### 4.4 Limitations and Future Work (300 words)
- **Computational limitations**: Cycle detection w bardzo dużych sieciach
- **Temporal resolution**: Snapshots co 50K steps (może przegapić krótkie cykle)
- **Validation**: Potrzeba eksperymentalnej weryfikacji hub molecules
- **Future work**: 
  - Temporal evolution analysis
  - Larger networks (10M+ steps)
  - Experimental validation of top cycles

---

### 5. Conclusions (~250 words)

- **Summary**: Autocatalytic networks są powszechne i różnorodne
- **Key findings**: 
  1. [X]% cycles są autocatalytic
  2. Topologia różni się między scenariuszami
  3. Hub molecules są konserwowane
- **Significance**: Mechanizmy emergencji złożoności w prebiotic chemistry
- **Future directions**: Eksperymentalna weryfikacja, większe sieci

---

## 📊 Figures Plan

### Figure 1: Cycle Type Distribution
- **Panel A**: Pie chart (direct, indirect, hypercycles)
- **Panel B**: Examples of each type (reaction diagrams)
- **Panel C**: Cycle length distribution (histogram)

### Figure 2: Amplification Factors
- **Panel A**: Violin plots by scenario
- **Panel B**: Amplification vs cycle length (scatter)
- **Panel C**: Temporal evolution of amplification (time series)

### Figure 3: Network Topology
- **Panel A**: Degree distribution (log-log plot, all scenarios)
- **Panel B**: Hub molecules (network visualization, top 20)
- **Panel C**: Network motifs (frequency comparison)

### Figure 4: Scenario Comparison
- **Panel A**: Cycle frequency comparison (bar chart)
- **Panel B**: Amplification factors comparison (box plots)
- **Panel C**: Hub molecule overlap (Venn diagram)
- **Panel D**: Topology metrics comparison (heatmap)

### Figure 5: Temporal Evolution (Optional)
- **Panel A**: Network growth over time (all scenarios)
- **Panel B**: Cycle emergence timeline
- **Panel C**: Hub molecule stability (persistence over time)

---

## 📋 Tables Plan

### Table 1: Summary Statistics by Scenario
- Cycle counts (total, direct, indirect, hypercycles)
- Amplification factors (mean, std, range)
- Network metrics (nodes, edges, density, clustering)
- Hub molecules (top 5 per scenario)

### Table 2: Top 10 Hub Molecules (All Scenarios)
- Molecule name, formula, SMILES
- Degree centrality
- Betweenness centrality
- Scenario(s) where found
- PubChem match (if available)

### Table S1: All Autocatalytic Cycles (Supplementary)
- Cycle ID, type, length
- Participants (molecules)
- Amplification factor
- Scenario, run ID
- First appearance (step)

### Table S2: Network Motifs (Supplementary)
- Motif type (3-node, 4-node)
- Frequency per scenario
- Z-score (statistical significance)
- Examples

---

## ✅ Analysis Checklist

### Data Analysis (Luty 2025)
- [ ] Classify all 769,315 cycles (direct, indirect, hypercycles)
- [ ] Calculate amplification factors for all cycles
- [ ] Build reaction networks (all scenarios)
- [ ] Calculate network topology metrics (degree, centrality, motifs)
- [ ] Identify hub molecules (top 10 per scenario)
- [ ] Statistical comparison between scenarios

### Figure Generation (Marzec 2025)
- [ ] Figure 1: Cycle type distribution
- [ ] Figure 2: Amplification factors
- [ ] Figure 3: Network topology
- [ ] Figure 4: Scenario comparison
- [ ] Figure 5: Temporal evolution (optional)

### Writing (Marzec-Kwiecień 2025)
- [ ] Abstract (250 words)
- [ ] Introduction (1500 words)
- [ ] Methods (1800 words)
- [ ] Results (2000 words)
- [ ] Discussion (1500 words)
- [ ] Conclusions (250 words)

### Finalization (Kwiecień-Maj 2025)
- [ ] Internal review
- [ ] Statistical validation
- [ ] Figure quality check
- [ ] References complete
- [ ] Supplementary materials
- [ ] Submission preparation

---

## 📅 Timeline

### Luty 2025 (Weeks 1-4)
- **Week 1-2**: Detailed cycle classification
- **Week 3-4**: Network topology analysis

### Marzec 2025 (Weeks 5-8)
- **Week 5-6**: Figure generation
- **Week 7-8**: Methods & Results writing

### Kwiecień 2025 (Weeks 9-12)
- **Week 9-10**: Introduction & Discussion writing
- **Week 11-12**: Finalization & internal review

### Maj 2025 (Weeks 13-16)
- **Week 13-14**: Final formatting & checks
- **Week 15-16**: Submission

---

**Last Updated**: 2025-12-04  
**Status**: Outline ready, ready to start analysis (styczeń 2026)

