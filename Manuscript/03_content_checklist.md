# Content Checklist: What to Include in Each Section

**Use this as a checklist while writing—tick off items as you complete them.**

---

## Abstract (250 words)

### Must Include:
- [ ] Problem statement (why DNA geometry matters)
- [ ] Method name: "Shape-Height Transform"
- [ ] Dataset size: "202 PDB structures"
- [ ] Three geometric units: Vij, Shambhian, Mamton
- [ ] **Key finding 1:** r=0.83 correlation (with p-value)
- [ ] **Key finding 2:** 91% classification accuracy
- [ ] **Key finding 3:** Blind nucleosome detection (1F66)
- [ ] Significance: "first systematic basepair-resolution framework"

### Avoid:
- ❌ Vague claims without numbers
- ❌ Excessive background (save for Intro)
- ❌ Citations (abstracts typically don't cite)

---

## 1. Introduction (1500-2000 words)

### 1.1 DNA Structure and Function (~400 words)

**Must Include:**
- [ ] Opening: DNA 3D structure regulates function
- [ ] Example 1: Nucleosome positioning affects accessibility
- [ ] Example 2: Protein-DNA recognition depends on shape
- [ ] Example 3: A-tracts create rigid regions
- [ ] Citations: Watson & Crick, Olson, Luger, Richmond

**Why This Matters:** Establishes that geometry is biologically important.

### 1.2 Current Methods and Limitations (~500 words)

**Must Include:**
- [ ] X-ray crystallography (atomic resolution, static)
- [ ] NMR (solution structure, ensemble)
- [ ] MD simulations (dynamics, computationally expensive)
- [ ] **Computational tools:**
  - [ ] Curves+ (helical parameters)
  - [ ] 3DNA (base-pair steps)
  - [ ] Canal (curvature prediction)
- [ ] **The gap:** Existing tools miss cross-sectional shape
- [ ] Citations: Drew 1981, Lavery & Sklenar 1989, Lu & Olson 2003

**Why This Matters:** Shows what's missing in the field (your opportunity).

### 1.3 Shape-Height Transform (~300 words)

**Must Include:**
- [ ] Concept: DNA as curve + tangent frames + cross-sections
- [ ] **Vij definition:** Local curvature ||C''(h)||
- [ ] **Shambhian definition:** Cross-sectional eccentricity
- [ ] **Mamton definition:** Skewness magnitude
- [ ] Advantages: basepair-scale, orthogonal, predictive
- [ ] Optional: Schematic diagram (Figure 1)

**Why This Matters:** Introduces your method before Results.

### 1.4 Research Questions (~300 words)

**Must Include:**
- [ ] **Question 1:** Are geometric units independent or coupled?
- [ ] **Question 2:** Do functional classes have distinct signatures?
- [ ] **Question 3:** Can geometry predict function without annotation?
- [ ] Brief mention of approach: "We analyzed 202 structures..."
- [ ] Roadmap: "Section 3.2 reports coupling, 3.5 reports classification..."

**Why This Matters:** Gives readers expectations for Results.

---

## 2. Methods (2000-2500 words)

### 2.1 Dataset Construction (~400 words)

**Must Include:**
- [ ] Data source: RCSB PDB (www.rcsb.org)
- [ ] **Selection criteria:**
  - [ ] Resolution <3.0 Å
  - [ ] DNA chain length >20 bp
  - [ ] No missing backbone atoms
- [ ] Final dataset: 202 structures
- [ ] Preprocessing: extracted DNA chains, assigned bases
- [ ] **Diversity statistics:**
  - [ ] Free duplex count
  - [ ] Protein-DNA complex count
  - [ ] Nucleosome count
  - [ ] Non-canonical count
- [ ] Reference: "Supplementary Table S1 lists all structures"

**Why This Matters:** Readers need to know what data you analyzed.

### 2.2 SHT Calculation (~600 words)

**Must Include:**
- [ ] **Height function h(s):**
  - [ ] Arc length parameterization
  - [ ] Cubic spline smoothing
- [ ] **Tangent recovery:**
  - [ ] T(h) = dC/dh
  - [ ] C''(h) = d²C/dh² (curvature vector)
- [ ] **Cross-sectional frame:**
  - [ ] {T, N₁, N₂} orthonormal basis
  - [ ] Atom projection onto (N₁, N₂) plane
- [ ] **Moment calculation:**
  - [ ] 2nd moment tensor: I_jk = Σ mᵢ (rᵢ)_j (rᵢ)_k
  - [ ] 3rd moment tensor: S_jkl
- [ ] Reference: "Full derivation in Supplementary Methods"

**Why This Matters:** Core of your method—must be reproducible.

### 2.3 Geometric Unit Definitions (~400 words)

**Must Include:**
- [ ] **Vij formula:** v₀(h) = ||C''(h)||
  - [ ] Units: [Length⁻¹]
  - [ ] Physical meaning: inverse radius of curvature
- [ ] **Shambhian formula:** η_r = 1 - λ₂/λ₁
  - [ ] Dimensionless, range [0,1]
  - [ ] Physical meaning: cross-sectional flatness
- [ ] **Mamton formula:** Ω = ||S||_F
  - [ ] Units: [Length³]
  - [ ] Physical meaning: skewness magnitude (roughness)
- [ ] **Aggregation:** Computed mean, max, std per structure

**Why This Matters:** Defines variables used in Results.

### 2.4 Statistical Analysis (~500 words)

**Must Include:**
- [ ] **Correlation analysis:**
  - [ ] Pearson correlation (SciPy 1.15.1)
  - [ ] Two-tailed significance testing
  - [ ] Bonferroni correction (α = 0.01/3)
- [ ] **Causal arrow testing:**
  - [ ] Random Forest regression (100 trees)
  - [ ] 5-fold cross-validation
  - [ ] R² scoring
- [ ] **Clustering:**
  - [ ] K-means in 3D [Vij, η, Ω] space
  - [ ] Elbow method + silhouette score for optimal K
  - [ ] Z-score normalization
- [ ] **Outlier detection:**
  - [ ] Z-score >2.5σ threshold
  - [ ] Per-unit independent testing

**Why This Matters:** Establishes rigor of analysis.

### 2.5 Functional Classification (~600 words)

**Must Include:**
- [ ] **Class definitions:**
  - [ ] Nucleosome: histone-wrapped DNA
  - [ ] Protein-bound: TF/polymerase complexes
  - [ ] Free duplex: canonical B-form/A-form
  - [ ] Non-canonical: quadruplex, junctions, damage
- [ ] **Training set:**
  - [ ] 33 manually annotated structures
  - [ ] Geometry-based heuristics for augmentation
  - [ ] Class distribution: 16/11/4/2
- [ ] **Model architecture:**
  - [ ] Random Forest (200 trees, depth 10)
  - [ ] Balanced class weights
  - [ ] Features: [Vij_mean, η_mean, Ω_mean]
  - [ ] Standardization via StandardScaler
- [ ] **Evaluation:**
  - [ ] Stratified 5-fold cross-validation
  - [ ] Metrics: accuracy, precision, recall, F1
  - [ ] Feature importance via Gini impurity

**Why This Matters:** Model must be reproducible by others.

---

## 3. Results (2500-3000 words, 6-8 figures)

### 3.1 Distribution (~400 words, Figure 1)

**Must Include:**
- [ ] Vij statistics: mean, std, range
- [ ] Shambhian statistics: mean, std, range
- [ ] Mamton statistics: mean, std, range
- [ ] Distribution shapes: "Vij right-skewed, Shambhian left-skewed..."
- [ ] **Figure 1:** Histograms + box plots for all three units
- [ ] Interpretation: "Most DNA has low curvature (Vij<0.15)..."

### 3.2 Vij-Mamton Coupling ⭐ (~600 words, Figure 2)

**Must Include:**
- [ ] **Primary result:** "r = 0.83, p = 1.6×10⁻¹⁷"
- [ ] 95% confidence interval: [0.72, 0.90]
- [ ] Sample size: N = 65 structures
- [ ] R² = 0.69 (69% variance explained)
- [ ] Linear fit equation: Mamton = 580·Vij - 45
- [ ] **Contrast:** Vij-Shambhian r=-0.09, p=0.50 (no correlation)
- [ ] **Causal arrows:**
  - [ ] Vij → Mamton: R² = 0.19
  - [ ] Mamton → Vij: R² = 0.41
- [ ] **Figure 2:**
  - [ ] Panel A: Vij vs Mamton scatter with fit line
  - [ ] Panel B: Vij vs Shambhian (null result)
  - [ ] Panel C: Correlation matrix heatmap
- [ ] Mechanism hypothesis: "Bending disrupts base stacking..."

**Why This Matters:** Your primary discovery—strongest result.

### 3.3 Geometric States (~500 words, Figure 3)

**Must Include:**
- [ ] K-means clustering: K=4 optimal
- [ ] Silhouette score: 0.52
- [ ] **Cluster statistics table:**
  - [ ] State 0: N=27, Vij=0.12, η=0.99, Ω=168 (normal B-DNA)
  - [ ] State 1: N=2, Vij=0.30, η=0.91, Ω=26366 (nucleosomes)
  - [ ] State 2: N=11, Vij=0.10, η=0.85, Ω=248 (protein-deformed)
  - [ ] State 3: N=25, Vij=0.08, η=0.98, Ω=150 (stiff duplex)
- [ ] **Key finding:** State 1 is 5σ separated (discrete cluster)
- [ ] **Figure 3:**
  - [ ] Panel A: 3D scatter colored by cluster
  - [ ] Panel B: Silhouette scores K=2-8
  - [ ] Panel C: PCA projection (PC1 vs PC2)
- [ ] Biological interpretation per state

### 3.4 Outliers (~600 words, Figure 4)

**Must Include:**
- [ ] 5 structures exceed 2.5σ threshold
- [ ] **Outlier table:**
  - [ ] 1AOI: Vij +4.7σ, Mamton +5.5σ → nucleosome
  - [ ] 1F66: Vij +4.4σ, Mamton +5.6σ → nucleosome
  - [ ] 1AHD: Shambhian -4.2σ → homeodomain
- [ ] **1F66 validation** ⭐:
  - [ ] "Predicted as nucleosome from geometry alone"
  - [ ] "PDB lookup confirmed: H2A.Z variant nucleosome"
  - [ ] "Suto et al. (2000), Nat Struct Biol"
- [ ] **1AHD validation:**
  - [ ] "Antennapedia homeodomain-DNA complex"
  - [ ] "Billeter et al. (1993)"
  - [ ] "Protein binding reduces cross-sectional circularity"
- [ ] **Figure 4:**
  - [ ] Panel A: Z-score plots
  - [ ] Panel B: 3D scatter with outliers labeled
  - [ ] Panel C: Structure visualizations (1AOI, 1F66, 1AHD)
  - [ ] Panel D: PDB metadata table

**Why This Matters:** Validates your method via blind prediction.

### 3.5 Classification Model (~600 words, Figures 5-6)

**Must Include:**
- [ ] **Cross-validation performance:**
  - [ ] Accuracy: 91% ± 12%
  - [ ] Per-fold: [100%, 71%, 83%, 100%, 100%]
- [ ] **Per-class metrics table:**
  - [ ] Free duplex: precision/recall/F1
  - [ ] Protein-bound: precision/recall/F1
  - [ ] Non-canonical: precision/recall/F1
  - [ ] Nucleosome: precision/recall/F1
- [ ] **Feature importance:**
  - [ ] Mamton: 41%
  - [ ] Shambhian: 33%
  - [ ] Vij: 26%
- [ ] **Figure 5:**
  - [ ] Panel A: CV accuracy per fold
  - [ ] Panel B: Confusion matrix 4x4
  - [ ] Panel C: Feature importance bars
  - [ ] Panel D: ROC curves (if multi-class AUC computed)
- [ ] Interpretation: "Mamton dominates predictive power..."

### 3.6 Predictions on Unlabeled (~400 words, Figure 6)

**Must Include:**
- [ ] 32 unlabeled structures predicted
- [ ] **Class distribution:**
  - [ ] 25 free duplex (78%)
  - [ ] 4 non-canonical (13%)
  - [ ] 3 protein-bound (9%)
  - [ ] 0 nucleosome
- [ ] **Top 10 high-confidence predictions:**
  - [ ] 1CF7: 94% free duplex
  - [ ] 1BF4: 91% free duplex
  - [ ] 1B8I: 89% free duplex
  - [ ] ... (list top 10)
- [ ] **Figure 6:**
  - [ ] Panel A: Predicted class counts (bar chart)
  - [ ] Panel B: Confidence histogram
  - [ ] Panel C: Top 10 table
  - [ ] Panel D: 3D scatter colored by prediction
- [ ] Validation where possible: "1EMH predicted protein-bound, confirmed TF complex"

---

## 4. Discussion (1800-2200 words)

### 4.1 Principal Findings (~300 words)

**Must Include:**
- [ ] Restate three key results:
  - [ ] Vij-Mamton coupling r=0.83
  - [ ] Four geometric states
  - [ ] 91% classification accuracy
- [ ] Novelty: "First basepair-resolution cross-sectional framework"
- [ ] Significance: "Enables functional prediction from geometry alone"

### 4.2 Mechanism of Coupling (~500 words)

**Must Include:**
- [ ] **Hypothesis 1:** Bending disrupts base stacking
  - [ ] Differential strain (inner vs outer edge)
  - [ ] Sugar pucker transitions (C2'-endo ↔ C3'-endo)
  - [ ] Creates asymmetric cross-section
- [ ] **Supporting evidence:**
  - [ ] A-tracts: low Vij + low Mamton (both rigid)
  - [ ] TA steps: high Vij + high Mamton (both flexible)
  - [ ] Nucleosomes: extreme both
- [ ] **Literature comparison:**
  - [ ] Cite MD studies showing bending-induced unwinding
  - [ ] Cite experimental work on sequence-dependent mechanics
- [ ] **Predictive application:** "Given Vij, estimate Mamton via Ω ≈ 580·v₀"

### 4.3 Geometric States (~500 words)

**Must Include:**
- [ ] Not continuous spectrum—discrete states
- [ ] **State 1 (nucleosomes):**
  - [ ] Discrete outlier, 5σ separated
  - [ ] Geometric signature of DNA wrapping
  - [ ] Could predict nucleosome occupancy from geometry
- [ ] **State 2 (protein-deformed):**
  - [ ] Intermediate geometry
  - [ ] Transcription factors, polymerases
  - [ ] May correlate with binding affinity
- [ ] **States 0 & 3 (canonical):**
  - [ ] Most genomic DNA
  - [ ] State 3 likely A-form or A-tracts
- [ ] **Biological significance:** State transitions require energy

### 4.4 Comparison to Existing Methods (~400 words)

**Must Include:**
- [ ] **Curves+:** Helical parameters (twist, roll)
  - [ ] Strengths: Comprehensive, well-established
  - [ ] Gap: Doesn't capture cross-sectional shape
- [ ] **3DNA:** Base-pair step analysis
  - [ ] Strengths: Atomic-level detail
  - [ ] Gap: Misses surface roughness
- [ ] **SHT advantages:**
  - [ ] Orthogonal geometric information
  - [ ] Basepair-scale resolution
  - [ ] Predictive capability (91% classifier)
- [ ] **Complementarity:** "SHT complements, not replaces"
- [ ] **Example:** "Curves+ cannot distinguish nucleosomes from bent DNA; SHT identifies via Mamton"

### 4.5 Applications (~600 words)

**Must Include:**
- [ ] **Application 1: Chromatin biology**
  - [ ] Predict nucleosome positioning from sequence→geometry
  - [ ] Compare to NuPoP, NuCLE (sequence-only tools)
  - [ ] SHT adds geometric constraints
- [ ] **Application 2: Protein-DNA docking**
  - [ ] Filter poses by geometric compatibility
  - [ ] Proteins prefer State 2 (moderate deformation)
  - [ ] Could improve docking scores
- [ ] **Application 3: DNA nanotechnology**
  - [ ] Design structures with controlled curvature
  - [ ] Avoid excessive roughness (use Vij-Mamton relationship)
  - [ ] Predict mechanical stability
- [ ] **Application 4: Drug design**
  - [ ] Small molecules targeting DNA geometry
  - [ ] Screen compounds for geometric effects in silico
- [ ] **Future experiments:**
  - [ ] FRET: Measure Vij via end-to-end distance
  - [ ] AFM: Directly measure roughness (compare to Mamton)
  - [ ] CD: Correlate with Shambhian (B vs A form)

### 4.6 Limitations (~400 words)

**Must Include:**
- [ ] **Limitation 1: Dataset size**
  - [ ] 202 structures (vs ~1500 total DNA in PDB)
  - [ ] Merged dataset only 65 structures
  - [ ] Mitigation: Ongoing expansion
- [ ] **Limitation 2: Resolution effects**
  - [ ] Requires high-resolution (<3 Å)
  - [ ] Coordinate errors affect geometric measurements
  - [ ] B-factors not accounted for
- [ ] **Limitation 3: Static snapshots**
  - [ ] Crystal structures are frozen
  - [ ] DNA breathes, fluctuates in solution
  - [ ] Future: Apply to MD trajectories
- [ ] **Limitation 4: Sequence-geometry gap**
  - [ ] Found correlations, not mechanisms
  - [ ] Need sequence models
  - [ ] Future: Deep learning to predict geometry from sequence

**Be honest but constructive:** State limitation + mitigation/future work.

---

## 5. Conclusions (~400 words)

**Must Include:**
- [ ] **Paragraph 1:** What you did
  - [ ] Developed SHT framework
  - [ ] Analyzed 202 structures
  - [ ] Extracted Vij, η, Ω
- [ ] **Paragraph 2:** Key findings (with numbers)
  - [ ] r=0.83 coupling
  - [ ] 91% accuracy
  - [ ] Blind nucleosome detection
- [ ] **Paragraph 3:** Validation
  - [ ] All outliers have biological explanations
  - [ ] 1F66 confirmed via PDB lookup
  - [ ] p<10⁻¹⁶ statistical significance
- [ ] **Paragraph 4:** Broader impact
  - [ ] First basepair-scale cross-sectional framework
  - [ ] Enables functional prediction
  - [ ] Opens avenues: chromatin, drug design, nanotech
- [ ] **Paragraph 5:** Future outlook
  - [ ] Expand to full PDB (~1500 structures)
  - [ ] Integrate sequence models
  - [ ] Experimental validation (FRET, AFM)
- [ ] **Final sentence:** Memorable closing
  - [ ] "The Shape-Height Transform provides a novel geometric lens for understanding DNA structure-function relationships at unprecedented basepair-scale resolution."

---

## Figures (6-8 main, 8+ supplementary)

### Main Figures (Must Have)

**Figure 1: Method Overview & Distributions**
- [ ] Panel A: SHT pipeline schematic
- [ ] Panel B: Vij histogram
- [ ] Panel C: Shambhian histogram
- [ ] Panel D: Mamton histogram

**Figure 2: Vij-Mamton Coupling** ⭐ MOST IMPORTANT
- [ ] Panel A: Vij vs Mamton scatter (r=0.83)
- [ ] Panel B: Vij vs Shambhian (null)
- [ ] Panel C: Shambhian vs Mamton (null)
- [ ] Panel D: Correlation matrix heatmap

**Figure 3: Geometric States**
- [ ] Panel A: 3D K-means clusters
- [ ] Panel B: Silhouette scores
- [ ] Panel C: PCA projection
- [ ] Panel D: Per-cluster statistics bars

**Figure 4: Outliers & Validation**
- [ ] Panel A: Z-score plots
- [ ] Panel B: 3D scatter with outliers labeled
- [ ] Panel C: Structure visualizations (1AOI, 1F66, 1AHD)
- [ ] Panel D: PDB metadata table

**Figure 5: Classification Performance**
- [ ] Panel A: CV accuracy per fold
- [ ] Panel B: Confusion matrix 4x4
- [ ] Panel C: Feature importance
- [ ] Panel D: ROC curves

**Figure 6: Predictions on Unlabeled**
- [ ] Panel A: Predicted class distribution
- [ ] Panel B: Confidence histogram
- [ ] Panel C: Top 10 high-confidence table
- [ ] Panel D: 3D scatter colored by prediction

**Optional Figure 7: Sequence Analysis** (if you do this)
- [ ] Panel A: AT% vs Vij
- [ ] Panel B: Dinucleotide step geometry
- [ ] Panel C: GC-rich vs AT-rich profiles

### Supplementary Figures

- [ ] **S1:** Per-structure geometric profiles (all 202)
- [ ] **S2:** Correlation analysis (all pairs, all 202 structures)
- [ ] **S3:** Clustering K=2 through K=8
- [ ] **S4:** Permutation feature importance
- [ ] **S5:** Model calibration curves
- [ ] **S6:** Learning curves (training set size vs accuracy)
- [ ] **S7:** Prediction confidence distributions
- [ ] **S8:** Additional structure visualizations

---

## Supplementary Tables

**Must Include:**
- [ ] **S1:** Complete dataset (202 PDB IDs, resolution, year, function, chains)
- [ ] **S2:** Geometric measurements (Vij, η, Ω mean/max/std for all 202)
- [ ] **S3:** Clustering assignments (PDB ID, cluster, distances to centroids)
- [ ] **S4:** Classification predictions (32 unlabeled, predicted class, confidence)
- [ ] **S5:** Software versions (Python 3.12.3, NumPy 2.2.4, etc.)

---

## References (~80-100 citations)

### Must Cite (High Priority)

**DNA Structure Classics:**
- [ ] Watson & Crick (1953) - Double helix
- [ ] Drew et al. (1981) - B-DNA dodecamer
- [ ] Dickerson (1989/1992) - DNA bending and curvature
- [ ] Olson et al. (1998) - DNA deformability

**Nucleosomes:**
- [ ] Luger et al. (1997) - 1AOI nucleosome structure ⭐
- [ ] Suto et al. (2000) - 1F66 H2A.Z variant ⭐
- [ ] Richmond & Davey (2003) - Nucleosome structure review
- [ ] Segal et al. (2006) - Nucleosome positioning from sequence

**Protein-DNA:**
- [ ] Billeter et al. (1993) - 1AHD homeodomain ⭐
- [ ] von Hippel (2007) - Protein-DNA interactions

**DNA Mechanics:**
- [ ] Calladine & Drew (1984) - Understanding DNA
- [ ] Lankas et al. (2003) - DNA flexibility from MD
- [ ] Lavery & Hartmann (2010) - Nucleic acid conformational analysis

**Computational Tools:**
- [ ] Lavery & Sklenar (1989) - Curves
- [ ] Lu & Olson (2003) - 3DNA
- [ ] Vlahovicek et al. (2003) - Canal

**Machine Learning:**
- [ ] Breiman (2001) - Random Forests
- [ ] Pedregosa et al. (2011) - scikit-learn

**Software:**
- [ ] Harris et al. (2020) - NumPy
- [ ] McKinney (2010) - Pandas
- [ ] Virtanen et al. (2020) - SciPy
- [ ] Hunter (2007) - Matplotlib

### Citation Checklist by Section

**Introduction:**
- [ ] 30-40 citations covering DNA structure, function, and existing tools

**Methods:**
- [ ] 10-15 citations for software, algorithms, statistical tests

**Results:**
- [ ] 5-10 citations comparing to specific prior findings

**Discussion:**
- [ ] 30-40 citations interpreting mechanisms, comparing to literature, proposing applications

---

## Writing Progress Tracker

### Week 1: Methods + Results
- [ ] Day 1: Methods 2.1-2.2 written
- [ ] Day 2: Methods 2.3-2.4 written
- [ ] Day 3: Methods 2.5 written
- [ ] Day 4: Results 3.1-3.3 written
- [ ] Day 5: Results 3.4-3.6 written
- [ ] Days 6-7: Revise Methods + Results

### Week 2: Figures
- [ ] Day 1: Generate Figures 1-3
- [ ] Day 2: Generate Figures 4-6
- [ ] Day 3: Polish all figures (fonts, colors, resolution)
- [ ] Day 4: Write all figure captions
- [ ] Days 5-7: Revise Results based on figures

### Week 3: Introduction + Discussion
- [ ] Day 1: Read 15 papers for Intro
- [ ] Day 2: Read 15 papers for Discussion
- [ ] Day 3: Write Introduction draft
- [ ] Day 4: Write Discussion sections 4.1-4.3
- [ ] Day 5: Write Discussion sections 4.4-4.6
- [ ] Days 6-7: Revise Intro + Discussion

### Week 4: Polish
- [ ] Day 1: Write Abstract + Conclusions
- [ ] Day 2: Format references (BibTeX)
- [ ] Day 3: Create supplementary materials
- [ ] Day 4: Format per NAR guidelines
- [ ] Day 5: Read through, check flow
- [ ] Days 6-7: Final polishing

### Week 5-6: Review
- [ ] Send to collaborators/advisor
- [ ] Incorporate feedback
- [ ] Address all comments

### Week 7: Submit
- [ ] Write cover letter
- [ ] Prepare submission files
- [ ] Upload to journal portal
- [ ] SUBMIT! 🎉

---

## Final Pre-Submission Checklist

### Content Complete
- [ ] Abstract: 250 words, includes key findings
- [ ] Introduction: 1500-2000 words, 30+ citations
- [ ] Methods: 2000-2500 words, fully reproducible
- [ ] Results: 2500-3000 words, 6-8 figures
- [ ] Discussion: 1800-2200 words, interprets findings
- [ ] Conclusions: 400 words, memorable closing
- [ ] References: 60-100 citations, properly formatted
- [ ] Supplementary: Tables S1-S5, Figures S1-S8

### Quality Checks
- [ ] No typos (spell-checked)
- [ ] Consistent notation (v₀/Vij, η/Shambhian, Ω/Mamton)
- [ ] All abbreviations defined on first use
- [ ] All figures referenced in text
- [ ] All citations formatted correctly
- [ ] Word count within NAR limits (<9000 words main text)
- [ ] Line numbers added
- [ ] Page numbers continuous

### Statistical Rigor
- [ ] All correlations include p-values
- [ ] All comparisons include effect sizes
- [ ] All sample sizes reported (N=...)
- [ ] All confidence intervals reported (95% CI)
- [ ] Multiple testing correction applied where needed

### Figures Quality
- [ ] All figures 300 DPI minimum
- [ ] Font sizes 12pt minimum
- [ ] Color-blind friendly palettes
- [ ] Clear axis labels
- [ ] Legends included
- [ ] Captions self-contained (explain without reading text)

### NAR-Specific
- [ ] Formatted per NAR author guidelines
- [ ] Data availability statement
- [ ] Author contributions section
- [ ] Competing interests declared
- [ ] Acknowledgments section
- [ ] Cover letter drafted

**You're ready to submit when ALL boxes are checked!**

---

## After Submission

### What to Expect
- **Initial review:** 2-4 weeks (editor decision: desk reject, review, or accept)
- **Peer review:** 4-8 weeks (2-3 reviewers provide feedback)
- **Revision:** 2-4 weeks (you address reviewer comments)
- **Final decision:** 1-2 weeks
- **Total timeline:** 3-6 months from submission to publication

### Possible Outcomes

**Accept (rare on first submission):** 🎉 Publish as is!

**Minor revisions (best realistic outcome):**
- [ ] Address specific comments
- [ ] Add requested figures/analyses
- [ ] Revise within 2-4 weeks
- [ ] Resubmit → Accept

**Major revisions (common):**
- [ ] Add significant new data/analysis
- [ ] Expand discussion
- [ ] Address conceptual issues
- [ ] Revise within 4-8 weeks
- [ ] Resubmit → possible second round

**Reject with resubmission option:**
- [ ] Substantial changes needed
- [ ] Essentially rewrite
- [ ] Consider as new submission

**Reject (if happens at NAR, try PLOS ONE or BMC Bioinformatics):**
- [ ] Don't give up!
- [ ] Revise based on feedback
- [ ] Submit to another journal

---

**YOU'VE GOT THIS! Follow the checklist, write section by section, and you'll have a paper ready to submit in 6-8 weeks.**
