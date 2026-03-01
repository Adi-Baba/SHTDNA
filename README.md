# SHTDNA — Shape-Height Transform for DNA

**Author:** Aditya (Adi-Baba)  
**Status:** Active Research  

---

## Overview

The **Shape-Height Transform (SHT)** is a novel geometric framework for analysing 3-D DNA structure using three orthogonal quantitative units derived from atomic coordinates in PDB files.

| Unit | Symbol | Measures |
|------|--------|---------|
| **Vij** (Curvature) | $v_0$ | Local bending of the helical axis |
| **Shambhian** (Eccentricity) | $\eta_r$ | Cross-sectional ellipticity |
| **Mamton** (Roughness) | $\Omega$ | Backbone height variance |

Key results (N = 65 merged dataset):
- Vij–Mamton coupling: **r = 0.83** (p < 10⁻¹⁵)  
- Random-forest classifier: **91% accuracy** distinguishing structural classes  
- All three units are statistically orthogonal (each captures a distinct geometric mode)

---

## Repository Structure

```
SHTDNA/
├── PDB/code/               # Core SHT computation scripts
├── CurvatureMode/          # Vij (curvature) extraction & validation
├── EccentricityMode/       # Shambhian (eccentricity) extraction & validation
├── RoughnessMode/          # Mamton (roughness) extraction & validation
├── BreathingMode/          # B-factor breathing analysis
├── DiscoveryMode/          # State-space reduction & trajectory analysis
├── NovelMode/              # Hydration prediction & novel analyses
├── PredictionMode/         # Nucleosome energy prediction
├── ThermodynamicsMode/     # Thermodynamic analysis
├── Spectral_Allostery_Discovery/  # Spectral allosteric signal analysis
├── Matrix_DNA_Analysis/    # Julia-based matrix DNA analyser
├── SHT_DNA/                # SHT literature extraction
├── Manuscript/             # Technical Monograph (LaTeX source → PDF)
└── Readme/                 # Extended documentation
```

> **Dataset not included.** Raw PDB structure files (8.9 GB) are not redistributed here. Download structures from [rcsb.org](https://www.rcsb.org) using the IDs listed in `PDB/dataset/analysis_summary.txt` or in Appendix A of the monograph.

---

## Quick Start

```bash
git clone git@github.com:Adi-Baba/SHTDNA.git
cd SHTDNA
python3 -m venv .venv && source .venv/bin/activate
pip install numpy scipy scikit-learn matplotlib seaborn biopython
```

### Compute all three SHT units for a PDB file

```bash
python3 PDB/code/sht_dna_analysis.py path/to/structure.pdb
```

### Reproduce the curvature (Vij) results

```bash
python3 CurvatureMode/extract_curvature.py
```

### Reproduce the full N = 65 merged analysis

```bash
python3 DiscoveryMode/extract_full_signature.py
```

---

## Technical Monograph

The full mathematical derivation, validation, and biological interpretation is in [`Manuscript/Technical_Monograph/`](Manuscript/Technical_Monograph/). Compile the PDF with:

```bash
cd Manuscript/Technical_Monograph
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

---

## License

MIT — see [LICENSE](LICENSE)
