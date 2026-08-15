import pandas as pd
import os

# Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

RANKINGS_FILE = os.path.join(PROCESSED_DIR, "casx_v6_probabilistic_rankings.csv")
GTEX_FILE = os.path.join(PROCESSED_DIR, "casx_v6_continuous_gtex.csv")

casx_v6 = pd.read_csv(RANKINGS_FILE)
gtex_df = pd.read_csv(GTEX_FILE)

pd.options.display.float_format = '{:.2f}'.format

print("\n" + "="*70)
print("  TABLE 1: CAS-X FINAL TOP 10 TARGETS")
print("="*70)
cols = ['GENE', 'CASX_v6_Final_Score', 'eQTL_Tissue_Count', 'eQTL_Max_Abs_Slope']
print(casx_v6[cols].head(10).to_string(index=False))

print("\n" + "="*70)
print("  TABLE 2: LEAD eQTL VARIANTS FOR TOP 5 TARGETS (COLOC PREP)")
print("="*70)
top_5 = casx_v6.head(5)['GENE'].tolist()
lead_variants = []
for gene in top_5:
    gene_data = gtex_df[gtex_df['gene_symbol'] == gene]
    if not gene_data.empty:
        lead_idx = gene_data['pval_nominal'].idxmin()
        lead = gene_data.loc[lead_idx]
        lead_variants.append({
            'Target Gene': gene,
            'Lead Variant (GTEx ID)': lead['variant_id'],
            'Tissue': lead['tissue'],
            'Slope': lead['slope'],
            'P-Value': f"{lead['pval_nominal']:.2e}"
        })

print(pd.DataFrame(lead_variants).to_string(index=False))
print("="*70 + "\n")
