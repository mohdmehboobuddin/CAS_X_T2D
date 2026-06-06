import pandas as pd

df = pd.read_csv("data/processed/evidence_matrix.csv")

# Initialize scores
df["GWAS_SCORE"] = 2

# Regulatory / Coding
df["REGULATORY_SCORE"] = df["REGULATORY_EVIDENCE"].map({
    "Coding variant with strong biological evidence": 3,
    "Coding variant with established metabolic role": 3,
    "Strong regulatory evidence": 3,
    "Enhancer-associated regulatory locus": 2,
    "Proximal regulatory locus": 1
})

# Expression
df["EXPRESSION_SCORE"] = df["EXPRESSION_EVIDENCE"].map({
    "Highly expressed in pancreatic islets": 2,
    "Expressed in pancreatic islets": 2,
    "Moderately expressed in pancreatic islets": 1,
    "Moderately expressed in metabolic tissues": 1,
    "Low expression in pancreatic islets": 0
})

# Functional
df["FUNCTIONAL_SCORE"] = df["FUNCTIONAL_EVIDENCE"].map({
    "Strong functional evidence": 3,
    "Functional evidence reported": 2,
    "Moderate functional evidence": 1
})

# CRISPR
df["CRISPR_SCORE"] = 2

# Placeholder until real datasets are added
df["EQTL_SCORE"] = 0
df["CHROMATIN_SCORE"] = 0

# CAS-X Total
score_cols = [
    "GWAS_SCORE",
    "REGULATORY_SCORE",
    "EQTL_SCORE",
    "CHROMATIN_SCORE",
    "EXPRESSION_SCORE",
    "FUNCTIONAL_SCORE",
    "CRISPR_SCORE"
]

df["CASX_TOTAL"] = df[score_cols].sum(axis=1)

# Priority Tier
def tier(score):
    if score >= 15:
        return "High"
    elif score >= 10:
        return "Moderate"
    else:
        return "Low"

df["PRIORITY_TIER"] = df["CASX_TOTAL"].apply(tier)

df.to_csv("results/tables/casx_rankings.csv", index=False)

print("\n=== CAS-X Rankings ===\n")
print(df[["SNP","GENE","CASX_TOTAL","PRIORITY_TIER"]])

print("\nSaved:")
print("results/tables/casx_rankings.csv")
