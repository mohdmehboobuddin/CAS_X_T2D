CAS-X: Multi-Tissue Prioritization Framework for Type 2 Diabetes Susceptibility Genes 
Overview

CAS-X (Candidate Actionability Scoring Framework) is a computational gene-prioritization framework developed to identify and rank candidate Type 2 Diabetes (T2D) susceptibility genes through the integration of GWAS signals, multi-tissue eQTL evidence, expression support, pathway enrichment, and external validation resources.

The framework was designed to move beyond single-tissue analyses by incorporating evidence across biologically relevant tissues involved in T2D pathogenesis, including pancreas, adipose tissue, skeletal muscle, and blood.

Study Workflow

GWAS Loci → Gene Mapping → GTEx eQTL Integration → Multi-Tissue Scoring → CAS-X Prioritization → External Validation → Pathway Enrichment

Data Sources
GWAS
GWAS Catalog Type 2 Diabetes loci
eQTL
GTEx v11
Expression Validation
GSE38642
GSE81608
External Validation
Open Targets Platform
GenomeCRISPR
Key Results
Candidate Gene Set
50 Type 2 Diabetes GWAS loci analyzed
37 candidate genes prioritized
Multi-Tissue Integration

Evidence collected across:

Pancreas
Adipose Subcutaneous Tissue
Adipose Visceral Tissue
Skeletal Muscle
Whole Blood
Coverage Improvement
+175% increase in candidate gene coverage relative to pancreas-only prioritization approaches
Top Ranked Genes (CAS-X V5)
FTO
JAZF1
VEGFA
NOTCH2
TSPAN8
Validation
Expression Support
GSE38642 validation
GSE81608 human pancreatic islet validation
Open Targets Validation

Validated genes:

GCK
HNF1B
HNF4A
KCNJ11
KCNQ1
PPARG
SLC30A8
TCF7L2
CRISPR Validation

Validated genes:

FTO
IRS1
KCNJ11
PPARG
TCF7L2
Pathway Enrichment

Significant pathways identified include:

Type II Diabetes Mellitus
Insulin Secretion
IRS Activation
Adipogenesis
AMPK Signaling
Energy Metabolism
Repository Structure
CAS_X_T2D/

├── data/
├── scripts/
├── results/
│   ├── figures/
│   └── tables/
├── docs/
├── README.md
├── requirements.txt
├── LICENSE
└── .gitignore
Author

Mohd Mehboob Uddin

Department of Life Sciences
A.V. College of Arts, Science and Commerce
Osmania University
Hyderabad, Telangana, India

Citation

If you use CAS-X in your research, please cite the associated publication when available.

License

This project is distributed under the MIT License.
