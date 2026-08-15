# CAS-X: A Probabilistic Multi-Tissue Functional Genomics Framework for Prioritizing Type 2 Diabetes Susceptibility Genes

## Overview

CAS-X (Candidate Actionability Scoring Framework) is an advanced computational gene-prioritization engine developed to identify and rank candidate Type 2 Diabetes (T2D) susceptibility genes. Moving beyond standard single-tissue models, CAS-X integrates genome-wide association study (GWAS) signals with multi-tissue eQTL evidence, utilizing expected regulatory effect mathematics to isolate systemic targets.

Unlike conventional approaches that focus heavily on the pancreas, CAS-X incorporates evidence from five biologically relevant metabolic tissues: pancreas, subcutaneous adipose, visceral adipose, skeletal muscle, and whole blood. This framework mathematically prioritizes highly pleiotropic, systemic drivers of metabolic disease over localized signals.

---

## Core Methodology

### 1. Expected Regulatory Effect (Phase 2.7)
To accurately quantify target causality, CAS-X calculates the Expected Regulatory Effect for each gene by integrating the causal probability of an eQTL with its magnitude of effect:
`Expected Effect = Posterior Inclusion Probability (PIP) × |aFC|`

### 2. Multi-Tissue Integration
Cross-tissue integration is achieved by identifying the maximum regulatory effect across the 5 target metabolic compartments, ensuring that robust systemic signals are preserved.

### 3. Empirical Specificity Benchmarking
To validate biological specificity, prioritized metabolic candidates were rigorously benchmarked against empirical negative control loci (e.g., Brain Cortex eQTLs) using the Mann-Whitney U rank test (p = 1.04e-299), proving the framework does not merely capture global hyper-expressed genes.

---

## Study Workflow

GWAS Loci Extraction → GTEx eQTL Mapping → Expected Effect Scoring (PIP × \|aFC|) → Multi-Tissue Aggregation → Top 15 Target Prioritization → LOTO Sensitivity Analysis & Negative Control Validation

---

## Key Findings

| Metric | Result |
|----------|----------|
| GWAS loci analyzed | 50 |
| Top-tier prioritized candidates | 15 |
| Biological tissues integrated | 5 |
| Coverage improvement (vs. pancreas-only) | **+133.3%** |
| Highly Systemic Targets Identified | JAZF1, FTO, GRB14 |

---

## Repository Structure

```text
CAS_X_T2D/
├── data/
│   ├── raw/                       # Raw inputs (GTEx v11 continuous matrices)
│   └── processed/                 # Intermediate processed outputs & rankings
├── src/                           # Complete automated Python pipeline
│   ├── 01_data_processing.py      # Parses GTEx data & filters metabolic tissues
│   ├── 02_scoring_engine.py       # Phase 2.7 math (Expected Effect)
│   ├── 03_validation_metrics.py   # Phase 4 Negative Control Benchmarking
│   ├── 04_visualization_main.py   # Generates 600 DPI LZW TIFFs (Main Text)
│   └── 05_visualization_si.py     # Generates 600 DPI LZW TIFFs (Supplementary)
├── results/
│   ├── main_figures/              # Final high-res manuscript figures
│   ├── supplementary_figures/     # Final supplementary figures
│   ├── main_tables/               # Tables 1 & 2 CSVs
│   └── supplementary_tables/      # Table S1 - S4 CSVs
├── README.md
└── requirements.txt


---

## Reproducibility

All analyses were performed using Python-based workflows and publicly available datasets. For immediate review, all fully processed outputs, dataset rankings, and high-resolution figures are readily available in the `results/` directory.

For full transparency, the complete modular codebase required to reproduce the prioritization framework from scratch is provided in the `src/` directory.

---

## Author

**Mohd Mehboob Uddin**

Department of Life Sciences
A.V. College of Arts, Science and Commerce
Osmania University
Hyderabad, Telangana, India

---

## Citation

If you use CAS-X in your research, please cite the associated publication when available.

---

## License

This repository is provided for academic and research purposes.
