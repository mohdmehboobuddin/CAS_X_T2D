import pandas as pd
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
GTEX_FILE = os.path.join(PROJECT_ROOT, "data", "processed", "casx_v6_continuous_gtex.csv")
RANKINGS_FILE = os.path.join(PROJECT_ROOT, "data", "processed", "casx_v6_probabilistic_rankings.csv")

gtex_df = pd.read_csv(GTEX_FILE)
casx_v6 = pd.read_csv(RANKINGS_FILE)
top_5 = casx_v6.head(5)['GENE'].tolist()

print("--- LEAD EQTL VARIANTS FOR TOP 5 SYSTEMIC TARGETS ---")
print("Extracting lead variants to cross-reference with GWAS peaks for COLOC analysis...\n")

for gene in top_5:
    gene_data = gtex_df[gtex_df['gene_symbol'] == gene]
    if not gene_data.empty:
        # Find the single variant with the strongest P-value across all tissues
        lead_idx = gene_data['pval_nominal'].idxmin()
        lead_variant = gene_data.loc[lead_idx]
        print(f"Target: {gene}")
        print(f"  Lead Variant (GTEx ID): {lead_variant['variant_id']}")
        print(f"  Driving Tissue: {lead_variant['tissue']}")
        print(f"  Effect Size (Slope): {lead_variant['slope']:.3f}")
        print(f"  Nominal P-Value: {lead_variant['pval_nominal']:.2e}\n")
