import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
import os

def run_robustness():
    print("======================================================")
    print(" RUNNING LOTO SENSITIVITY & PCA VARIANCE ANALYSIS")
    print("======================================================")

    # 1. Load Data
    raw_eqtl = pd.read_csv("CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv")
    casx_scores = pd.read_csv("CAS_X_v2/results/tables/Broad_CASX_Prioritization_Scores.csv")
    
    # 2. Rebuild Matrix
    raw_eqtl['eqtl_signal'] = raw_eqtl['pip'] * raw_eqtl['afc'].abs()
    idx_max = raw_eqtl.groupby(['gene_symbol', 'tissue'])['eqtl_signal'].idxmax()
    matrix = raw_eqtl.loc[idx_max].pivot(index='gene_symbol', columns='tissue', values='eqtl_signal').fillna(0.0)
    core_tissues = ['Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Whole_Blood']
    matrix = matrix.reindex(columns=core_tissues).fillna(0.0)

    # 3. PCA Scree Plot Generation (Fig S6 at 600 DPI)
    scaler = StandardScaler()
    Z_matrix = scaler.fit_transform(matrix)
    pca = PCA()
    pca.fit(Z_matrix)
    explained_variance = pca.explained_variance_ratio_ * 100
    cumulative_variance = np.cumsum(explained_variance)

    sns.set_theme(style="whitegrid")
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.bar(range(1, 6), explained_variance, color='#4575b4', alpha=0.7, label='Individual Variance')
    ax1.plot(range(1, 6), cumulative_variance, color='#d73027', marker='o', lw=2, label='Cumulative Variance')
    for i, cum_var in enumerate(cumulative_variance):
        ax1.annotate(f"{cum_var:.1f}%", (i+1, cum_var + 2), ha='center', fontweight='bold', color='#d73027')
    ax1.set_xlabel('Principal Components', fontweight='bold')
    ax1.set_ylabel('Percentage of Variance Explained (%)', fontweight='bold')
    ax1.set_title('PCA Scree Plot: Multi-Tissue eQTL Variance', fontweight='bold')
    ax1.legend(loc='lower right')
    
    os.makedirs("CAS_X_v2/results/figures", exist_ok=True)
    plt.savefig("CAS_X_v2/results/figures/FigS6_PCA_Scree_Plot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Saved PCA Scree Plot (Fig S6) at 600 DPI.")

    # 4. Leave-One-Tissue-Out (LOTO) Analysis
    # BUG FIX: Align the original ranks explicitly to the alphabetical matrix index
    original_ranks = casx_scores.set_index('Target_Gene')['CASX_Score'].reindex(matrix.index).rank(ascending=False)
    
    loto_results = []
    print("\nExecuting LOTO Sensitivity...")
    
    for dropped_tissue in core_tissues:
        loto_tissues = [t for t in core_tissues if t != dropped_tissue]
        loto_matrix = matrix[loto_tissues]
        
        loto_Z = StandardScaler().fit_transform(loto_matrix)
        loto_pca = PCA()
        loto_PC = loto_pca.fit_transform(loto_Z)
        loto_var = loto_pca.explained_variance_ratio_
        
        K = 4
        loto_raw = []
        for g_idx in range(len(loto_matrix)):
            weighted_sq_sum = np.sum(loto_var[:K] * (loto_PC[g_idx, :K] ** 2))
            loto_raw.append(np.sqrt(weighted_sq_sum))
            
        loto_final = MinMaxScaler(feature_range=(0, 100)).fit_transform(np.array(loto_raw).reshape(-1, 1)).flatten()
        loto_ranks = pd.Series(loto_final, index=matrix.index).rank(ascending=False)
        
        # Calculate properly aligned Spearman correlation
        spearman_corr, _ = spearmanr(original_ranks, loto_ranks)
        loto_results.append({'Dropped_Tissue': dropped_tissue.replace('_', ' '), 'Spearman_rho': spearman_corr})
        print(f" -> Dropped {dropped_tissue}: Spearman rho = {spearman_corr:.3f}")

    # 5. Plot LOTO Sensitivity (Fig 4A at 600 DPI)
    loto_df = pd.DataFrame(loto_results).sort_values(by='Spearman_rho', ascending=True)
    
    fig, ax2 = plt.subplots(figsize=(8, 5))
    sns.barplot(data=loto_df, x='Spearman_rho', y='Dropped_Tissue', color='#3182bd', ax=ax2)
    ax2.set_xlim([0.0, 1.0])
    ax2.set_xlabel('Spearman Rank Correlation (ρ)', fontweight='bold')
    ax2.set_ylabel('Excluded Tissue', fontweight='bold')
    ax2.set_title('Leave-One-Tissue-Out (LOTO) Sensitivity Analysis', fontweight='bold', pad=15)
    
    for i, v in enumerate(loto_df['Spearman_rho']):
        ax2.text(v - 0.05, i, f"{v:.3f}", color='white', fontweight='bold', va='center')
        
    plt.savefig("CAS_X_v2/results/figures/Fig4A_LOTO_Sensitivity_Analysis.png", dpi=600, bbox_inches='tight')
    plt.close()
    
    pd.DataFrame(loto_results).to_csv("CAS_X_v2/results/tables/LOTO_Sensitivity_Stats.csv", index=False)
    print("-> Saved LOTO Sensitivity Plot (Fig 4A) at 600 DPI.")
    print("======================================================")

if __name__ == "__main__":
    run_robustness()
