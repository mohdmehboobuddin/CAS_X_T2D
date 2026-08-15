import pandas as pd

print("\n=== CAS-X GSE38642 Expression Extraction ===\n")

# CAS-X candidate genes and probe IDs
gene_map = {
    7909681: "PROX1",
    8112740: "ZBED3",
    8138789: "JAZF1",
    7964927: "TSPAN8",
    7919095: "NOTCH2",
    7946853: "KCNJ11",
    7995655: "FTO",
    8119898: "VEGFA"
}

# Read GSE38642 matrix
df = pd.read_csv(
    "data/raw/expression/GSE38642_series_matrix.txt.gz",
    sep="\t",
    skiprows=65,
    comment="!"
)

# Rename probe column
df = df.rename(columns={"ID_REF": "PROBE"})

# Convert probe IDs to numeric
df["PROBE"] = pd.to_numeric(df["PROBE"], errors="coerce")

# Keep only CAS-X probes
hits = df[df["PROBE"].isin(gene_map.keys())].copy()

# Add gene names
hits["GENE"] = hits["PROBE"].map(gene_map)

# Expression columns
sample_cols = [
    c for c in hits.columns
    if c.startswith("GSM")
]

# Calculate mean expression
hits["MEAN_EXPRESSION"] = hits[sample_cols].mean(axis=1)

# Final table
result = hits[
    ["GENE", "PROBE", "MEAN_EXPRESSION"]
].sort_values(
    by="MEAN_EXPRESSION",
    ascending=False
)

print(result)

# Save results
result.to_csv(
    "results/tables/gse38642_expression_support.csv",
    index=False
)

print("\nSaved:")
print("results/tables/gse38642_expression_support.csv")
