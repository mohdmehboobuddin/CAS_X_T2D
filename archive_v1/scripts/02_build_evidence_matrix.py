import pandas as pd

df = pd.read_csv("data/processed/master_t2d_loci.csv")

# Initialize evidence columns

df["GWAS_SCORE"] = 0
df["REGULATORY_SCORE"] = 0
df["EQTL_SCORE"] = 0
df["CHROMATIN_SCORE"] = 0
df["EXPRESSION_SCORE"] = 0
df["FUNCTIONAL_SCORE"] = 0
df["CRISPR_SCORE"] = 0

df["CASX_TOTAL"] = 0
df["PRIORITY_TIER"] = "Unassigned"
df["CRISPR_STRATEGY"] = "Unassigned"

output_file = "data/processed/casx_master_dataset.csv"

df.to_csv(output_file, index=False)

print("\nCAS-X evidence matrix created.")
print(f"Rows: {len(df)}")
print(f"Saved to: {output_file}")
