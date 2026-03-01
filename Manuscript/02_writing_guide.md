# How to Write Each Section: Practical Guide

**Goal:** Turn your analysis into a publishable NAR paper in 6-8 weeks

---

## General Writing Tips

### Style Guidelines

**DO:**
- ✅ Use active voice: "We discovered" (not "It was discovered that")
- ✅ Be concise: Cut unnecessary words
- ✅ Use present tense for established facts: "DNA consists of..."
- ✅ Use past tense for your experiments: "We analyzed 202 structures..."
- ✅ Define abbreviations on first use: "Shape-Height Transform (SHT)"
- ✅ Cite liberally: Every claim needs a reference (except your own results)
- ✅ Use subheadings to guide readers

**DON'T:**
- ❌ Overclaim: "revolutionizes" → "provides new insights"
- ❌ Use jargon without explanation
- ❌ Repeat yourself (unless for emphasis)
- ❌ Include raw code or equations without context
- ❌ Cite yourself excessively (looks narcissistic)

### Paragraph Structure

**Golden Formula:**
1. **Topic sentence**: State the main point
2. **Evidence/Explanation**: 2-4 sentences supporting the point
3. **Transition**: Link to next paragraph

**Example:**
> DNA bending induces cross-sectional irregularity. [Topic] We observed a strong positive correlation (r=0.83, p<10⁻¹⁶) between curvature (Vij) and roughness (Mamton) across 65 structures. [Evidence] This relationship likely arises from differential sugar pucker transitions at bent regions, which asymmetrically deform the backbone. [Explanation] These findings suggest that local mechanical stress propagates to surface geometry. [Transition]

---

## Section-by-Section Writing Guide

---

## 1. Abstract (Write LAST)

### When to Write
**After everything else is done.** Abstract is a summary—write it when you know what to summarize.

### How to Write
1. Open your Results section
2. Copy the 3 most important findings (numbers, p-values)
3. Write one sentence per:
   - Problem (why does this matter?)
   - Method (what did you do?)
   - Results (what did you find?—be specific with statistics)
   - Significance (why should anyone care?)

### Word Budget
- Problem: 40-50 words
- Method: 60-70 words
- Results: 80-100 words (MOST IMPORTANT)
- Significance: 40-50 words

### Check Quality
- ✅ Does it include your r=0.83 coupling?
- ✅ Does it mention 91% accuracy?
- ✅ Does it say "Shape-Height Transform"?
- ✅ Can someone understand your contribution without reading further?

### Example Opening Sentences
**Good:** "DNA structural diversity regulates gene expression, but systematic basepair-resolution geometric characterization remains lacking."

**Bad:** "DNA is important for biology and has many functions."

---

## 2. Introduction

### When to Write
**Week 3** (after Methods and Results, when you know your story)

### How to Write

#### Step 1: Read 20-30 papers (4-6 hours)
**Search queries:**
- "DNA bending sequence-dependent"
- "nucleosome structure Luger"
- "DNA mechanical properties"
- "Curves+ 3DNA helical parameters"
- "protein-DNA recognition"

**Take notes:** For each paper, write one sentence: "Author (Year) found [X], which relates to our work because [Y]."

#### Step 2: Build argument funnel (2 hours)
- **Wide:** DNA structure matters (everyone agrees)
- **Narrower:** Existing tools miss cross-sectional geometry (gap in field)
- **Narrow:** Our method fills this gap (your contribution)

#### Step 3: Write section by section (6-8 hours)

**Section 1.1: DNA Structure and Function (~400 words)**
- Start broad: "DNA three-dimensional structure regulates..."
- Cite early classics (Watson & Crick, Drew & Dickerson)
- Give 2-3 examples where structure affects function
- End with: "Therefore, understanding DNA geometry is critical."

**Template:**
```
DNA three-dimensional structure governs [function X] and [function Y] [cite 1,2].
For example, [specific case A] demonstrates [mechanism] [cite 3].
Additionally, [specific case B] shows [consequence] [cite 4,5].
These examples illustrate that...
```

**Section 1.2: Current Methods (~500 words)**
- Paragraph 1: X-ray, NMR, MD (what they can do)
- Paragraph 2: Curves+, 3DNA (existing computational tools)
- Paragraph 3: **The gap:** "However, these methods focus on [X] and miss [Y]."
  
**Key move:** Show respect for prior work while identifying limitation.

**Good:** "Curves+ provides comprehensive helical parameter analysis [cite], but does not characterize cross-sectional shape."

**Bad:** "Curves+ is inadequate and fails to..." (too harsh)

**Section 1.3: SHT Framework (~300 words)**
- Paragraph 1: High-level concept (curve → tangent → cross-section)
- Paragraph 2: Three units (one sentence each: Vij, η, Ω)
- Paragraph 3: Advantages (basepair-scale, orthogonal, predictive)

**Tip:** Include a schematic figure here (Figure 1: SHT pipeline diagram)

**Section 1.4: Research Questions (~300 words)**
- State 2-3 specific questions:
  1. "Are Vij, η, and Ω independent or coupled?"
  2. "Do functional classes have distinct geometric signatures?"
  3. "Can geometry predict function without annotation?"

**End with roadmap:** "To address these questions, we applied SHT to 202 PDB structures and performed [methods X, Y, Z]."

---

## 3. Methods

### When to Write
**Week 1** (FIRST section to write—easiest because you know what you did)

### How to Write

#### Rule: Reproducibility
Another scientist should be able to **exactly replicate** your analysis from this section.

#### Structure
Each subsection = one major step in your pipeline.

**Template:**
```
### 2.X [Step Name]

[One sentence explaining purpose]

[2-3 paragraphs describing procedure]

[Software/parameters: bullet list]
```

### Writing Each Subsection

**2.1 Dataset Construction**

**Paragraph 1: Selection criteria**
```
We downloaded DNA structures from RCSB PDB (www.rcsb.org, accessed January 2026)
applying the following filters: (1) resolution <3.0 Å, (2) chain length >20 bp,
(3) no missing backbone atoms. This yielded 202 structures (Supplementary Table S1).
```

**Tip:** Be **specific**. Don't say "high-quality structures"—say "resolution <3.0 Å"

**Paragraph 2: Preprocessing**
```
DNA chains were extracted using BioPython 1.79 [cite]. We computed the phosphate
backbone trace and assigned base identity via hydrogen bonding patterns with
threshold distance 3.5 Å.
```

**Paragraph 3: Diversity**
- List functional classes with counts
- State resolution range (min, max, median)
- Reference supplementary table

**2.2 SHT Calculation**

**Paragraph 1: Height function**
```
DNA backbone coordinates were parameterized by arc length s, computed as cumulative
distance between adjacent phosphates. The height function h(s) represents the
z-coordinate along the principal helical axis, smoothed via cubic spline interpolation
with knot spacing of one basepair.
```

**Tip:** Include math equations for key formulas (Vij, η, Ω definitions)

**Format equations:**
```
$$
v_0(h) = \left\| \frac{d^2\mathbf{C}}{dh^2} \right\|
$$
```

**Paragraph 2-4:** Describe tangent recovery, cross-sectional frame, moment calculation

**Complexity balance:** Enough detail to reproduce, not so much readers drown in math.

**Solution:** Put full derivation in Supplementary Methods, give overview here.

**2.3 Geometric Unit Definitions**

One paragraph per unit (Vij, η, Ω):
1. Mathematical definition (equation)
2. Physical meaning (what does it measure?)
3. Units and range
4. Per-structure aggregation (mean, max, std)

**2.4 Statistical Analysis**

Be formulaic:
```
Pairwise Pearson correlations were computed using SciPy 1.15.1 with two-tailed
significance testing. We applied Bonferroni correction for multiple comparisons
(α = 0.01/3 = 0.0033). 
```

**List software with versions:**
- Python 3.12.3
- NumPy 2.2.4
- scikit-learn 1.6.1
- etc.

**2.5 Functional Classification**

**Paragraph 1:** Define classes (nucleosome, protein-bound, free duplex, non-canonical)

**Paragraph 2:** Training set construction
- "We manually annotated 33 structures..."
- Describe annotation rules
- State class distribution

**Paragraph 3:** RF architecture
```
We trained a Random Forest classifier (scikit-learn 1.6.1) with 200 trees,
maximum depth 10, and balanced class weights to handle imbalanced classes.
Features were standardized via z-score normalization. Performance was assessed
via stratified 5-fold cross-validation.
```

**Key:** Enough detail that someone could code it from scratch.

---

## 4. Results

### When to Write
**Week 2** (after Methods, while analysis is fresh)

### How to Write

#### Golden Rules
1. **Start each subsection with the main finding**
2. **Support with statistics** (numbers, p-values, confidence intervals)
3. **Reference figures explicitly**: "Figure 2A shows..."
4. **Interpret, don't just describe**: Explain *what it means*

#### Structure Template

```
### 3.X [Finding Name]

[Opening sentence: main result]

[Figure reference and description]

[Statistical details: correlation, p-value, effect size]

[Interpretation: biological meaning]

[Comparison to literature if relevant]
```

### Writing Each Subsection

**3.1 Distribution of Geometric Units**

**Opening:** "We first characterized the distribution of Vij, η, and Ω across 202 structures."

**Figure call:** "Figure 1A-C shows histograms for each unit."

**Statistics:** "Vij ranged from 0.04 to 0.31 L⁻¹ (mean 0.11 ± 0.05), with a right-skewed distribution indicating most structures have low curvature."

**Interpretation:** "This distribution reflects the predominance of gently curved B-form DNA in crystal structures."

**3.2 Vij-Mamton Coupling** ⭐ **MOST IMPORTANT**

**Opening:** "We discovered a strong positive correlation between curvature and roughness."

**Statistics FIRST:** "Vij and Mamton exhibited r = 0.83 (95% CI: [0.72, 0.90], p = 1.6×10⁻¹⁷, N=65 structures; Figure 2A)."

**Effect size:** "This represents a large effect (Cohen's d = 2.1), with 69% of Mamton variance explained by Vij (R² = 0.69)."

**Contrast:** "In contrast, Vij-Shambhian showed no correlation (r = -0.09, p = 0.50), indicating independence."

**Mechanism:** "We hypothesize that bending disrupts base stacking, creating asymmetric sugar-phosphate geometry that increases skewness."

**Tip:** Use **bold** or *italics* for key numbers: "**r = 0.83, p < 10⁻¹⁶**"

**3.3 Clustering**

**Opening:** "Geometric units define four distinct structural states."

**Method summary:** "K-means clustering (K=4) in [Vij, η, Ω] space..."

**Results table:** Include inline table (see structure document)

**Key finding:** "Notably, State 1 (nucleosomes) forms a **discrete outlier cluster** 5σ from the mean."

**Interpretation:** "This geometric separation suggests nucleosome wrapping represents a distinct mechanical phase, not simply extreme normal duplex."

**3.4 Outlier Validation** ⭐ **DEMONSTRATES UTILITY**

**Opening:** "To validate predictive power, we identified structures exceeding 2.5σ in any geometric unit."

**Table:** List 5 outliers with z-scores

**Case study:** "1F66 exhibited Vij = +4.4σ and Mamton = +5.6σ. **Without prior knowledge**, we predicted it as a nucleosome."

**Validation:** "PDB lookup confirmed 1F66 is 'nucleosome core particle containing H2A.Z variant' [cite Suto 2000]."

**Impact:** "This blind detection demonstrates SHT's ability to identify functional classes from geometry alone."

**3.5 Classification Model**

**Opening:** "We trained a Random Forest classifier to predict functional class from [Vij, η, Ω]."

**Performance:** "Cross-validation yielded **91% accuracy** (±12%, 5-fold stratified CV)."

**Per-class:** Give precision/recall table

**Feature importance:** "Mamton contributed 41% of predictive power, followed by Shambhian (33%) and Vij (26%)."

**Interpretation:** "This indicates surface roughness is the primary discriminator of functional state."

**3.6 Predictions**

**Opening:** "We applied the trained model to 32 unlabeled structures."

**Results:** "25 were classified as free duplex (78%), 4 as non-canonical (13%), and 3 as protein-bound (9%)."

**Validation:** "Top predictions (confidence >85%) showed agreement with PDB metadata where available."

---

## 5. Discussion

### When to Write
**Week 3** (after Introduction, when literature review is fresh)

### How to Write

#### Purpose
- **Interpret** results (why did you find what you found?)
- **Compare** to literature (how does this fit existing knowledge?)
- **Speculate** on mechanisms (what biological processes explain patterns?)
- **Project** future work (what's next?)

#### Structure
1. Restate main findings (1 paragraph)
2. Interpret each major result (1-2 paragraphs each)
3. Compare to existing tools (1 paragraph)
4. Applications (2-3 paragraphs)
5. Limitations (1-2 paragraphs)

### Writing Style

**DON'T just repeat Results.** Discuss *implications*.

**Results:** "We found r=0.83 correlation between Vij and Mamton."

**Discussion:** "The Vij-Mamton coupling likely arises from mechanical stress propagation. When DNA bends, differential strain on the inner vs outer edges disrupts regular base stacking [cite MD study]. This asymmetry increases third-order moments (skewness), manifesting as higher Mamton. Supporting this, A-tracts—known to resist bending [cite]—exhibit both low Vij AND low Mamton, confirming the coupling extends to sequence-dependent mechanics."

**See the difference?** Discussion connects your finding to mechanism and literature.

### Subsection Templates

**5.1 Principal Findings**
```
We developed the Shape-Height Transform to systematically characterize DNA
geometry at basepair resolution. Across 202 structures, we discovered [finding 1],
[finding 2], and [finding 3]. These results establish SHT as [significance].
```

**5.2 Mechanistic Interpretation**
```
The Vij-Mamton coupling suggests [mechanism]. This is consistent with [literature
finding X], which demonstrated [related result]. We propose that [specific molecular
process] underlies this relationship. Future MD simulations could test this by
[experimental prediction].
```

**5.3 Comparison to Existing Methods**
```
Our approach complements tools like Curves+ [cite], which excels at [their strength]
but does not capture [your advantage]. By focusing on cross-sectional geometry,
SHT provides orthogonal information enabling [new capability]. For example,
Curves+ cannot distinguish nucleosomes from other bent DNA, whereas SHT identifies
them via the Mamton signature.
```

**5.4 Applications**
```
The predictive power of geometric signatures opens several applications. First,
[application 1]: By training sequence→geometry models, one could predict [outcome]
from genomic DNA, with implications for [biological problem]. Second, [application 2]...
```

**5.5 Limitations**
```
Several caveats warrant discussion. First, [limitation 1]: Our dataset comprises
[size/scope], which may not represent [broader population]. However, [mitigation or
future plan]. Second, [limitation 2]...
```

**Be honest but not self-defeating.**

**Good:** "Our classifier requires high-resolution structures (<3 Å), limiting application to crystal structures. Future work could extend SHT to lower-resolution data via coarse-graining."

**Bad:** "Our method is severely limited and probably doesn't work well." (undermines your own work)

---

## 6. Conclusions

### When to Write
**Week 4** (near the end)

### How to Write

**Short and punchy.** 3-5 paragraphs, ~400 words.

**Template:**
1. **Summary:** What you did
2. **Key findings:** 2-3 most important results (with numbers)
3. **Validation:** How you know it works
4. **Impact:** Why it matters
5. **Future:** One sentence on next steps

**Example:**
```
We developed the Shape-Height Transform framework for systematic analysis of DNA
cross-sectional geometry at basepair resolution. Applying SHT to 202 PDB structures,
we discovered quantitative coupling between curvature and roughness (r=0.83, p<10⁻¹⁶)
and achieved 91% accuracy in functional classification from three geometric parameters.

Blind prediction of nucleosome structures (1AOI, 1F66) validated SHT's predictive
power, demonstrating that geometry alone—without sequence or annotation—suffices
to identify functional states. This represents a novel capability not achievable
with existing tools.

The SHT framework opens new avenues for structure-based DNA bioinformatics, with
applications in chromatin analysis, protein-DNA docking, and nanotechnology.
Future work expanding to genomic-scale prediction could revolutionize our understanding
of DNA structure-function relationships.
```

**End strong.** Last sentence should be memorable.

---

## 7. Making Figures

### General Rules
- **6-8 main figures** (most important results)
- **Supplementary figures** (everything else)
- High resolution: 300 DPI minimum
- Clear labels: 12pt font minimum
- Color-blind friendly palettes

### Figure Types

**Figure 1: Overview/Schematic**
- Panel A: SHT pipeline diagram
- Panel B-D: Definitions of Vij, η, Ω (visual)
- Purpose: Help readers understand method

**Figure 2: Distributions**
- Panel A-C: Histograms of Vij, η, Ω
- Panel D: Joint distribution heatmap
- Purpose: Show data range and spread

**Figure 3: Correlation** ⭐
- Panel A: Vij vs Mamton scatter (r=0.83)
- Panel B: Vij vs Shambhian (no correlation)
- Panel C: Correlation matrix heatmap
- Purpose: Main discovery

**Figure 4: Clustering**
- Panel A: 3D scatter colored by cluster
- Panel B: Silhouette scores
- Panel C: Per-cluster statistics
- Purpose: Show geometric states

**Figure 5: Outliers**
- Panel A: Z-score plot
- Panel B: 3D scatter with outliers labeled
- Panel C: Example structures (1AOI, 1F66, 1AHD)
- Purpose: Validation

**Figure 6: Classification**
- Panel A: Confusion matrix
- Panel B: Feature importance
- Panel C: ROC curves
- Purpose: Model performance

### Figure Captions

**Formula:** [What it shows] + [How it was made] + [Key result]

**Example:**
```
Figure 3. Strong coupling between DNA curvature and surface roughness.
(A) Scatter plot of Vij vs Mamton for 65 structures with both measurements
(blue dots). Linear fit (red line) shows r=0.83, p=1.6×10⁻¹⁷. (B) Vij vs
Shambhian shows no correlation (r=-0.09, p=0.50). (C) Correlation matrix
heatmap confirms Vij-Mamton coupling and Shambhian independence. Values
shown are Pearson r; * indicates p<0.01.
```

---

## 8. References

### How to Cite

**Use reference manager:** Zotero (free), Mendeley, EndNote

**Steps:**
1. As you read papers, add to library
2. Export as BibTeX when done
3. Use LaTeX or Word with citation plugin

### Citation Style

**NAR uses:** Author-date (Harvard style)

**In-text:** "Previous work (Luger et al., 1997) showed..."

**Reference list:**
```
Luger, K., Mäder, A.W., Richmond, R.K., Sargent, D.F. and Richmond, T.J. (1997)
Crystal structure of the nucleosome core particle at 2.8 Å resolution.
Nature, 389, 251-260.
```

### How Many References?

**NAR typical:** 60-100 references

**Breakdown:**
- Introduction: 30-40 (broad literature)
- Methods: 10-15 (software, algorithms)
- Results: 5-10 (specific comparisons)
- Discussion: 30-40 (interpretation and applications)

### What to Cite

**Always cite:**
- Every factual claim not your own
- Software packages with versions
- Statistical methods
- Prior work you build on
- Conflicting findings (discuss differences)

**Never cite:**
- Common knowledge ("DNA is a double helix")
- Your own unpublished work excessively
- Preprints unless necessary (prefer peer-reviewed)

---

## 9. Supplementary Materials

### What Goes in SI?

**Main text:** Results readers MUST see to understand contribution

**SI:** Details for specialists who want to replicate

**Move to SI:**
- Full 202-structure data table
- Per-structure plots
- Additional statistical tests
- Algorithm pseudocode
- Parameter sensitivity analysis
- All software versions

**Keep in main:**
- Summary statistics
- Key correlations
- Classification performance
- Main figures (6-8)

---

## 10. Timeline & Workflow

### Week-by-Week Plan

**Week 1: Methods + Results Rough Draft**
- Day 1-2: Write Methods 2.1-2.3 (dataset, SHT calculation)
- Day 3-4: Write Methods 2.4-2.5 (statistics, classification)
- Day 5-7: Write Results 3.1-3.6 (one subsection per day)
- **Output:** First draft of Methods + Results (~5000 words)

**Week 2: Figures**
- Day 1: Generate Figure 1-3
- Day 2: Generate Figure 4-6
- Day 3: Polish figures (fonts, colors, labels)
- Day 4: Write figure captions
- Day 5-7: Revise Results based on figures
- **Output:** 6-8 publication-quality figures

**Week 3: Introduction + Discussion**
- Day 1-2: Literature review (read 20-30 papers, take notes)
- Day 3-4: Write Introduction (~2000 words)
- Day 5-7: Write Discussion (~2000 words)
- **Output:** Intro + Discussion drafts

**Week 4: Polish + Format**
- Day 1: Write Abstract + Conclusions
- Day 2: Format references (BibTeX)
- Day 3: Create supplementary materials
- Day 4: Format according to NAR guidelines
- Day 5-7: Read through, fix typos, check flow
- **Output:** Complete draft ready for review

**Week 5-6: Internal Review**
- Send to colleagues/advisor
- Incorporate feedback
- Final polish

**Week 7: Submit**
- Prepare cover letter
- Upload to journal portal
- Submit!

---

## 11. Common Pitfalls to Avoid

### Content Issues

**❌ Overclaiming**
- Bad: "We revolutionized DNA structural analysis"
- Good: "We developed a novel framework that complements existing tools"

**❌ Underselling**
- Bad: "We found a weak correlation that might be interesting"
- Good: "We discovered a strong coupling (r=0.83, p<10⁻¹⁶) between..."

**❌ Burying the lead**
- Bad: Putting r=0.83 finding in middle of paragraph 7
- Good: Opening Results 3.2 with "We discovered a strong positive correlation..."

**❌ Data without interpretation**
- Bad: "Figure 2 shows Vij vs Mamton. The correlation is 0.83."
- Good: "Strong Vij-Mamton coupling (r=0.83, Figure 2) suggests bending induces roughness via base stacking disruption."

### Writing Issues

**❌ Passive voice**
- Bad: "It was found that DNA bending was correlated with roughness"
- Good: "We found that DNA bending correlates with roughness"

**❌ Vague statements**
- Bad: "Many structures showed high curvature"
- Good: "15 structures (7%) exceeded Vij > 0.20"

**❌ Unexplained jargon**
- Bad: "We computed η_r via eigendecomposition of I"
- Good: "We computed cross-sectional flatness (η_r) as the eccentricity of the inertia tensor"

**❌ Missing citations**
- Bad: "DNA bending is sequence-dependent"
- Good: "DNA bending is sequence-dependent (Calladine & Drew, 1984)"

### Formatting Issues

**❌ Inconsistent notation**
- Use v₀ OR Vij throughout, not both

**❌ Poor figure quality**
- Minimum 300 DPI, readable fonts

**❌ Missing details**
- Always include N, p-values, confidence intervals

---

## 12. Final Checklist Before Submission

### Content Completeness
- [ ] Abstract summarizes all key findings
- [ ] Introduction cites 30+ papers
- [ ] Methods include software versions
- [ ] Results present statistics with p-values
- [ ] Discussion interprets findings (not just repeat)
- [ ] Conclusions summarize impact
- [ ] All figures referenced in text
- [ ] All citations formatted correctly

### Quality Checks
- [ ] No typos (run spell-check)
- [ ] Consistent notation (v₀ vs Vij)
- [ ] All abbreviations defined on first use
- [ ] Figures have clear captions
- [ ] Supplementary materials complete
- [ ] Word count within journal limits (<8000 words)
- [ ] References >60, <100

### NAR-Specific Requirements
- [ ] Formatted per NAR author guidelines
- [ ] Line numbers added
- [ ] Continuous page numbering
- [ ] Author affiliations complete
- [ ] Cover letter drafted
- [ ] Data availability statement
- [ ] Competing interests declared

---

## You're Ready When...

✅ You can explain your paper in 30 seconds (elevator pitch)
✅ Every figure tells a clear story
✅ Every claim has statistical support (p-values)
✅ A non-expert can understand the Abstract
✅ You've cited the major papers in your field
✅ You can defend every choice you made

**Now submit and move on to the next project!**
