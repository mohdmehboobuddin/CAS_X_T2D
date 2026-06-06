import pandas as pd

df = pd.read_csv("data/processed/evidence_matrix.csv")

print("\n=== CAS-X Evidence Matrix ===\n")

# Show loci that already have some evidence
filled = df.dropna(how="all", subset=[
    "GWAS_EVIDENCE",
    "REGULATORY_EVIDENCE",
    "EQTL_EVIDENCE",
    "CHROMATIN_EVIDENCE",
    "EXPRESSION_EVIDENCE",
    "FUNCTIONAL_EVIDENCE",
    "CRISPR_EVIDENCE"
])

print(filled)

print(f"\nTotal loci: {len(df)}")
print(f"Loci with evidence: {len(filled)}")
