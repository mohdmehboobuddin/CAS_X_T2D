import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "results", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

GTEX_FILE = os.path.join(PROCESSED_DIR, "casx_v6_continuous_gtex.csv")
RANKINGS_FILE = os.path.join(PROCESSED_DIR, "casx_v6_probabilistic_rankings.csv")

gtex_df = pd.read_csv(GTEX_FILE)
casx_v6 = pd.read_csv(RANKINGS_FILE)
gtex_df['abs_slope'] = gtex_df['slope'].abs()

tissues = gtex_df['tissue'].unique()
top_10 = casx_v6.head(10)['GENE'].tolist()

# Calculate LOTO ranks
loto_ranks = {gene: [] for gene in top_10}

for dropped_tissue in tissues:
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
    
    for gene in top_10:
        loto_ranks[gene].append(merged[merged['GENE'] == gene]['Rank'].values[0])

# Prepare data for plotting
plot_data = []
for i, gene in enumerate(top_10):
    original_rank = float(i + 1)
    avg_loto_rank = np.mean(loto_ranks[gene])
    plot_data.append({'Gene': gene, 'Rank Type': 'Original Rank', 'Rank': original_rank})
    plot_data.append({'Gene': gene, 'Rank Type': 'Average LOTO Rank', 'Rank': avg_loto_rank})

plot_df = pd.DataFrame(plot_data)

# Generate Figure
plt.figure(figsize=(10, 6))
sns.set_theme(style="whitegrid")

ax = sns.barplot(
    data=plot_df, 
    x='Gene', 
    y='Rank', 
    hue='Rank Type', 
    palette=['#2B5B84', '#A3C1DA'],
    edgecolor='black'
)

# Invert Y-axis so Rank 1 is at the top
plt.ylim(15, 0)
plt.ylabel('CAS-X Priority Rank (Lower is Better)', fontsize=12, fontweight='bold')
plt.xlabel('Top 10 Prioritized Genes', fontsize=12, fontweight='bold')
plt.title('Leave-One-Tissue-Out (LOTO) Sensitivity Analysis', fontsize=14, fontweight='bold', pad=15)
plt.legend(title='', loc='lower right', frameon=True)

# Clean up layout and save at 500 DPI
plt.tight_layout()
output_file = os.path.join(FIGURES_DIR, "LOTO_Sensitivity_Analysis.png")
plt.savefig(output_file, dpi=500, bbox_inches='tight')
print(f"LOTO figure successfully saved to: {output_file}")
