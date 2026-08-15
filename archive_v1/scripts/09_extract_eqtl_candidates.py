import pandas as pd

# Load pancreas GTEx dataset
df = pd.read_parquet(
    "data/raw/eqtl/extracted/Pancreas.v11.eQTLs.SuSiE_summary.parquet"
)

# CAS-X T2D genes
t2d_genes = [
    "SLC30A8","KCNJ11","PPARG","TCF7L2","CDKAL1",
    "IGF2BP2","MTNR1B","IRS1","KLF14","HHEX",
    "HNF1B","NOTCH2","JAZF1","ADAMTS9","GCK",
    "GCKR","HNF4A","CDKN2A","FTO","ADCY5",
    "ZFAND6","KCNQ1","VEGFA","ZBED3","CDKN2B",
    "CRY2","KCNK16","GRB14","ZMIZ1","FAF1",
    "CHCHD9","PROX1","TSPAN8","CENTD2"
]

# Find matching genes
hits = df[df["gene_name"].isin(t2d_genes)]

# Extract strongest PIP per gene
results = (
    hits.groupby("gene_name")["pip"]
    .max()
    .reset_index()
    .sort_values("pip", ascending=False)
)

# Save results
results.to_csv(
    "results/tables/eqtl_candidates.csv",
    index=False
)

print("\n=== CAS-X eQTL Candidates ===\n")
print(results)

print("\nSaved:")
print("results/tables/eqtl_candidates.csv")
