import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, precision_recall_curve, auc
import os

def plot_validation():
    print("Generating final validation ROC/PR curves...")
    
    # Load data
    raw_eqtl = pd.read_csv("CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv")
    casx_scores = pd.read_csv("CAS_X_v2/results/tables/Broad_CASX_Prioritization_Scores.csv")
    
    with open("data/processed/opentargets_t2d_core_genes.txt", "r") as f:
        ot_genes = set(line.strip() for line in f if line.strip())

    # Rebuild baselines
    raw_eqtl['eqtl_signal'] = raw_eqtl['pip'] * raw_eqtl['afc'].abs()
    idx_max = raw_eqtl.groupby(['gene_symbol', 'tissue'])['eqtl_signal'].idxmax()
    matrix = raw_eqtl.loc[idx_max].pivot(index='gene_symbol', columns='tissue', values='eqtl_signal').fillna(0.0)
    
    df = casx_scores[['Target_Gene', 'CASX_Score']].copy()
    df['Label'] = df['Target_Gene'].apply(lambda x: 1 if x in ot_genes else 0)
    df['Pancreas_Signal'] = df['Target_Gene'].map(matrix['Pancreas']).fillna(0.0)
    
    # Calculate curves
    fpr_casx, tpr_casx, _ = roc_curve(df['Label'], df['CASX_Score'])
    roc_auc_casx = auc(fpr_casx, tpr_casx)
    
    fpr_panc, tpr_panc, _ = roc_curve(df['Label'], df['Pancreas_Signal'])
    roc_auc_panc = auc(fpr_panc, tpr_panc)

    # Plotting
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.plot(fpr_casx, tpr_casx, color='#2c7fb8', lw=2.5, label=f'CAS-X Framework (AUC = {roc_auc_casx:.2f})')
    ax.plot(fpr_panc, tpr_panc, color='#f03b20', lw=2, linestyle='--', label=f'Pancreas-Only eQTL (AUC = {roc_auc_panc:.2f})')
    ax.plot([0, 1], [0, 1], color='gray', lw=1.5, linestyle=':', label='Random Expectation (AUC = 0.50)')
    
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontweight='bold')
    ax.set_title('Receiver Operating Characteristic: De Novo Discovery', fontweight='bold', pad=15)
    ax.legend(loc="lower right", frameon=True)
    
    os.makedirs("CAS_X_v2/results/figures", exist_ok=True)
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/FigS7_Quantitative_Validation.png", dpi=300)
    print("SUCCESS: Saved to CAS_X_v2/results/figures/FigS7_Quantitative_Validation.png")

if __name__ == "__main__":
    plot_validation()
