import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

GTEX_FILE = os.path.join(PROCESSED_DIR, "casx_v6_continuous_gtex.csv")
RANKINGS_FILE = os.path.join(PROCESSED_DIR, "casx_v6_probabilistic_rankings.csv")

gtex_df = pd.read_csv(GTEX_FILE)
casx_v6 = pd.read_csv(RANKINGS_FILE)
gtex_df['abs_slope'] = gtex_df['slope'].abs()

tissues = gtex_df['tissue'].unique()
top_10_original = casx_v6.head(10)['GENE'].tolist()
print(f"Tracking Original Top 10 Genes: {top_10_original}\n")

# Store the rank of each gene across all LOTO iterations
loto_ranks = {gene: [] for gene in top_10_original}

for dropped_tissue in tissues:
    print(f"Running model WITHOUT {dropped_tissue}...")
    loto_df = gtex_df[gtex_df['tissue'] != dropped_tissue].copy()
    
    agg = loto_df.groupby('gene_symbol').agg(
        eQTL_Tissue_Count=('tissue', 'nunique'),
        eQTL_Max_Abs_Slope=('abs_slope', 'max'),
        eQTL_Min_Pval=('pval_nominal', 'min')
    ).reset_index()
    agg['eQTL_NegLog10_Pval'] = -np.log10(agg['eQTL_Min_Pval'] + 1e-300)
    
    merged = pd.merge(casx_v6[['GENE']], agg, left_on='GENE', right_on='gene_symbol', how='left').fillna(0)
    
    X = merged[['eQTL_Tissue_Count', 'eQTL_Max_Abs_Slope', 'eQTL_NegLog10_Pval']].values
    if X.sum() > 0:
        X_scaled = StandardScaler().fit_transform(X)
        merged['Raw_Score'] = PCA(n_components=1).fit_transform(X_scaled)
    else:
        merged['Raw_Score'] = 0
        
    merged['Score'] = ((merged['Raw_Score'] - merged['Raw_Score'].min()) / (merged['Raw_Score'].max() - merged['Raw_Score'].min())) * 100
    merged['Rank'] = merged['Score'].rank(ascending=False, method='min')
    
    for gene in top_10_original:
        rank = merged[merged['GENE'] == gene]['Rank'].values[0]
        loto_ranks[gene].append(rank)

print("\n--- LEAVE-ONE-TISSUE-OUT (LOTO) RESULTS ---")
print(f"{'GENE':<10} | {'Max Rank Drop':<15} | {'Average Rank':<15} | {'Robustness'}")
print("-" * 65)
for gene, ranks in loto_ranks.items():
    max_drop = max(ranks)
    avg_rank = np.mean(ranks)
    status = "Highly Robust" if avg_rank <= 15 else "Sensitive"
    print(f"{gene:<10} | {max_drop:<15.1f} | {avg_rank:<15.1f} | {status}")
