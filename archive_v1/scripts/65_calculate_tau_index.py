import pandas as pd
import numpy as np
import os

# Load the local eQTL matrix
data_path = "data/processed/casx_v6_continuous_gtex.csv"
if not os.path.exists(data_path):
    print(f"Error: Could not find {data_path}")
    exit()

df = pd.read_csv(data_path)

top_15 = ["JAZF1", "FTO", "GRB14", "CDKN2B", "ZBED3", "VEGFA", "ZFAND6", 
          "NOTCH2", "KLF14", "IRS1", "GCK", "ADCY5", "SUGP1", "PROX1", "CDKAL1"]

# Filter to top 15 genes
df_top = df[df['gene_symbol'].isin(top_15)].copy()

# Ensure we are looking at absolute effect sizes
df_top['abs_slope'] = df_top['slope'].abs()

# PIVOT THE DATA: Compress hundreds of SNP rows into 1 row per gene.
# We take the maximum absolute slope for each tissue to represent the gene's strongest effect.
pivot_df = df_top.pivot_table(
    index='gene_symbol', 
    columns='tissue', 
    values='abs_slope', 
    aggfunc='max'
).fillna(0)

# Calculate Tau Index
# Formula: Tau = sum(1 - (x_i / x_max)) / (N - 1)
results = []
for gene, row in pivot_df.iterrows():
    x_max = row.max()
    if x_max == 0:
        tau = np.nan
    else:
        x_hat = row / x_max
        tau = np.sum(1 - x_hat) / (len(pivot_df.columns) - 1)
    
    # A standard Tau threshold: < 0.8 is considered broadly expressed/systemic
    results.append({
        "Gene": gene,
        "Max_Tissue_Effect": row.idxmax(),
        "Tau_Index": round(tau, 3),
        "Profile": "Systemic/Pleiotropic" if tau < 0.8 else "Tissue-Specific"
    })

results_df = pd.DataFrame(results).sort_values(by="Tau_Index")

print("\n=== Tissue Specificity Index (Tau) for Top CAS-X Targets ===\n")
print(results_df.to_string(index=False))

# Save the output
out_dir = "results/tables"
os.makedirs(out_dir, exist_ok=True)
out_path = f"{out_dir}/TableS4_Tau_Index.csv"
results_df.to_csv(out_path, index=False)
print(f"\nSaved Tau Index table to: {out_path}")
