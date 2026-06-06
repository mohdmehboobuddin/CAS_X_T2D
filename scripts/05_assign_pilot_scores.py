import pandas as pd

# Load pilot evidence table
df = pd.read_csv("data/processed/pilot_evidence_table.csv")

# Replace empty cells with 0 for scoring
score_columns = [
    "GWAS_SCORE",
    "REGULATORY_SCORE",
    "EQTL_SCORE",
    "CHROMATIN_SCORE",
    "EXPRESSION_SCORE",
    "FUNCTIONAL_SCORE",
    "CRISPR_SCORE"
]

for col in score_columns:
    df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

# Calculate CAS-X total
df["CASX_TOTAL"] = df[score_columns].sum(axis=1)

# Assign priority tier
def assign_tier(score):
    if score >= 15:
        return "High"
    elif score >= 10:
        return "Moderate"
    else:
        return "Low"

df["PRIORITY_TIER"] = df["CASX_TOTAL"].apply(assign_tier)

# Save results
output_file = "data/processed/casx_master_dataset.csv"
df.to_csv(output_file, index=False)

print("\n=== CAS-X Pilot Ranking ===\n")
print(df[["SNP", "GENE", "CASX_TOTAL", "PRIORITY_TIER"]])

print(f"\nSaved to: {output_file}")
