import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def generate_core_figures():
    print("======================================================")
    print(" GENERATING CORE MANUSCRIPT FIGURES (FIG 2, 3) ")
    print("======================================================")

    scores_df = pd.read_csv("CAS_X_v2/results/tables/Broad_CASX_Prioritization_Scores.csv")
    eqtl_df = pd.read_csv("CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv")

    os.makedirs("CAS_X_v2/results/figures", exist_ok=True)

    # 1. Figure 2: Multi-Tissue eQTL Landscape (6-Tissue Heatmap)
    top15_genes = scores_df['Target_Gene'].head(15).tolist()
    eqtl_df['eqtl_signal'] = eqtl_df['pip'] * eqtl_df['afc'].abs()
    idx_max = eqtl_df.groupby(['gene_symbol', 'tissue'])['eqtl_signal'].idxmax()
    matrix = eqtl_df.loc[idx_max].pivot(index='gene_symbol', columns='tissue', values='eqtl_signal')
    
    matrix_top15 = matrix.reindex(top15_genes)
    
    # Align to 6 tissues including Liver
    core_tissues = ['Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Whole_Blood', 'Liver']
    matrix_top15 = matrix_top15.reindex(columns=core_tissues)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    cmap = sns.color_palette("mako", as_cmap=True)
    cmap.set_bad(color='#e0e0e0') # Light gray for missing data
    
    sns.heatmap(matrix_top15, cmap=cmap, annot=True, fmt='.2f', cbar_kws={'label': 'PIP × |aFC| Signal'}, ax=ax)
    
    # Removed "Figure 2:" prefix
    ax.set_title('6-Tissue Regulatory Landscape (Top 15 Targets)', fontweight='bold', pad=15)
    ax.set_xlabel('Metabolic Tissues', fontweight='bold')
    ax.set_ylabel('Target Genes', fontweight='bold')
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    
    # Removed _v2 from filename
    plt.savefig("CAS_X_v2/results/figures/Fig2_Multi-Tissue_eQTL_Landscape.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated 6-Tissue Heatmap at 600 DPI.")

    # 2. Figure 3: Top 15 Prioritized Targets Barplot
    fig, ax = plt.subplots(figsize=(10, 6))
    top15_df = scores_df.head(15).sort_values(by='CASX_Score', ascending=True)
    
    ax.barh(top15_df['Target_Gene'], top15_df['CASX_Score'], color='#2b8cbe')
    ax.set_xlim([0, 105])
    ax.set_xlabel('CAS-X Actionability Score (0-100)', fontweight='bold')
    ax.set_ylabel('Target Gene', fontweight='bold')
    
    # Removed "Figure 3:" prefix
    ax.set_title('Top 15 Prioritized De Novo T2D Targets', fontweight='bold', pad=15)
    
    for i, v in enumerate(top15_df['CASX_Score']):
        ax.text(v + 1.5, i, f"{v:.1f}", va='center', fontweight='bold', color='black')
        
    plt.tight_layout()
    
    # Removed _v2 from filename
    plt.savefig("CAS_X_v2/results/figures/Fig3_Top15_Prioritized_Targets.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated Top 15 Targets Barplot at 600 DPI.")

    print("SUCCESS: Core biological figures compiled.")
    print("======================================================")

if __name__ == "__main__":
    generate_core_figures()
