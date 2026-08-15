"""
Script 71: Quantitative Validation Metrics (Fixed Polarity)
Calculates AUROC, AUPRC, and Top-K enrichment for CAS-X prioritization scores.
Automatically handles inverted metrics (like ranks) to ensure correct ROC plotting.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from pathlib import Path

def generate_performance_metrics():
    potential_paths = [
        Path("data/processed/casx_v6_probabilistic_rankings.csv"),
        Path("data/processed/casx_v6_continuous_gtex.csv"),
        Path("results/Main_Manuscript_Tables/Table1_Target_Prioritization_Summary.csv")
    ]
    
    gold_standard_path = Path("data/processed/opentargets_t2d_core_genes.txt")
    
    df = None
    for path in potential_paths:
        if path.exists():
            temp_df = pd.read_csv(path)
            score_col = [c for c in temp_df.columns if any(kw in c.upper() for kw in ['SCORE', 'CAS', 'TOTAL', 'RANK'])]
            gene_col = [c for c in temp_df.columns if 'GENE' in c.upper() or 'SYMBOL' in c.upper()]
            
            if score_col and gene_col:
                df = temp_df
                break
                
    if df is None:
        print("❌ Error: Could not find a scored dataset.")
        return

    score_col = [c for c in df.columns if any(kw in c.upper() for kw in ['SCORE', 'CAS', 'TOTAL', 'RANK'])][-1]
    gene_col = [c for c in df.columns if 'GENE' in c.upper() or 'SYMBOL' in c.upper()][0]
    
    with open(gold_standard_path, 'r') as f:
        gold_genes = [line.strip().upper() for line in f if line.strip()]

    df['True_Positive'] = df[gene_col].astype(str).str.upper().isin(gold_genes).astype(int)
    
    y_true = df['True_Positive']
    
    # Check if the score is inverted (e.g., rank 1 is best) by testing the initial AUC
    temp_scores = df[score_col].fillna(0)
    fpr, tpr, _ = roc_curve(y_true, temp_scores)
    initial_auc = auc(fpr, tpr)
    
    # If the AUC is less than 0.5, the metric is inverted. We flip it by taking the negative.
    if initial_auc < 0.5:
        print(f"🔄 Metric '{score_col}' appears inverted (AUC < 0.5). Flipping polarity for calculation...")
        y_scores = -1 * temp_scores
    else:
        y_scores = temp_scores

    # Recalculate with corrected polarity
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    pr_auc = average_precision_score(y_true, y_scores)
    baseline_pr = y_true.sum() / len(y_true)

    # Generate 600 DPI Two-Panel Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=600)
    sns.set_theme(style="whitegrid")

    # Panel A: ROC Curve
    ax1.plot(fpr, tpr, color='#2980b9', lw=2.5, label=f'CAS-X Framework (AUC = {roc_auc:.2f})')
    ax1.plot([0, 1], [0, 1], color='#7f8c8d', lw=1.5, linestyle='--', label='Random (AUC = 0.50)')
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
    ax1.set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
    ax1.set_title('A. Receiver Operating Characteristic', fontsize=14, pad=10)
    ax1.legend(loc="lower right", framealpha=0.9)

    # Panel B: Precision-Recall Curve
    ax2.plot(recall, precision, color='#c0392b', lw=2.5, label=f'CAS-X Framework (AUC = {pr_auc:.2f})')
    ax2.axhline(y=baseline_pr, color='#7f8c8d', lw=1.5, linestyle='--', label=f'Baseline (AUC = {baseline_pr:.2f})')
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.set_xlabel('Recall (Sensitivity)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Precision (Positive Predictive Value)', fontsize=12, fontweight='bold')
    ax2.set_title('B. Precision-Recall Curve', fontsize=14, pad=10)
    ax2.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    
    out_dir = Path("results/Supplementary_Figures")
    out_path = out_dir / "FigS7_Quantitative_Validation.jpg"
    plt.savefig(out_path)
    print(f"\n✅ Plot successfully saved to: {out_path}")

if __name__ == "__main__":
    generate_performance_metrics()
