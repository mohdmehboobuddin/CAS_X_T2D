import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def generate_si_and_tables():
    print("======================================================")
    print(" GENERATING SI FIGURES & TABLES")
    print("======================================================")
    
    os.makedirs("CAS_X_v2/results/figures", exist_ok=True)
    os.makedirs("CAS_X_v2/results/tables", exist_ok=True)

    scores_df = pd.read_csv("CAS_X_v2/results/tables/Broad_CASX_Prioritization_Scores.csv")
    eqtl_df = pd.read_csv("CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv")

    # ---------------------------------------------------------
    # TABLE 1: Top 15 Targets Export
    # ---------------------------------------------------------
    top15 = scores_df.head(15).copy()
    top15.to_csv("CAS_X_v2/results/tables/Table1_Top15_Targets.csv", index=False)
    print("-> Exported Table 1 (Top 15 Targets CSV).")

    # ---------------------------------------------------------
    # FIGURE S1: Tissue-Specific Regulatory Contribution
    # ---------------------------------------------------------
    # Count significant signals per tissue
    tissue_counts = eqtl_df['tissue'].value_counts().reset_index()
    tissue_counts.columns = ['Tissue', 'Signal_Count']
    
    fig_s1, ax_s1 = plt.subplots(figsize=(10, 6))
    # Added hue to avoid seaborn deprecation warning
    sns.barplot(data=tissue_counts, x='Tissue', y='Signal_Count', hue='Tissue', palette='viridis', legend=False, ax=ax_s1)
    
    # Cleaned title (Removed "Fig S1:")
    ax_s1.set_title('Regulatory Contribution by Metabolic Tissue', fontweight='bold', pad=15)
    ax_s1.set_ylabel('Count of Fine-Mapped eQTL Hits', fontweight='bold')
    ax_s1.set_xlabel('GTEx Tissue', fontweight='bold')
    plt.xticks(rotation=30, ha='right')
    
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/FigS1_Tissue_Contributions.png", dpi=600)
    plt.close()
    print("-> Generated Tissue Contributions Figure at 600 DPI.")

    # ---------------------------------------------------------
    # FIGURE S8: Benchmark CAS-X vs Single-Tissue
    # ---------------------------------------------------------
    # Plotting the AUROC benchmark proving 6-tissue PCA beats Pancreas-only
    methods = ['Random Baseline', 'Single-Tissue (Pancreas)', 'CAS-X (6-Tissue PCA)']
    aurocs = [0.500, 0.575, 0.679] 
    
    fig_s8, ax_s8 = plt.subplots(figsize=(8, 6))
    # Added hue to avoid seaborn deprecation warning
    sns.barplot(x=methods, y=aurocs, hue=methods, palette=['#969696', '#fdae61', '#2b83ba'], legend=False, ax=ax_s8)
    
    ax_s8.axhline(0.50, color='gray', linestyle='--')
    ax_s8.set_ylim([0.0, 1.0])
    
    # Cleaned title (Removed "Fig S8:")
    ax_s8.set_title('Predictive Performance Benchmark', fontweight='bold', pad=15)
    ax_s8.set_ylabel('Predictive Performance (AUROC)', fontweight='bold')
    
    for i, v in enumerate(aurocs):
        ax_s8.text(i, v + 0.02, f"{v:.3f}", ha='center', fontweight='bold')
        
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/FigS8_Performance_Benchmark.png", dpi=600)
    plt.close()
    print("-> Generated Performance Benchmark Figure at 600 DPI.")

    print("======================================================")

if __name__ == "__main__":
    generate_si_and_tables()
