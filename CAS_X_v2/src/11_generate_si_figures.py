import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def generate_remaining_si_figures():
    print("======================================================")
    print(" GENERATING FINAL SI FIGURES (S2, S3, S4, S5)")
    print("======================================================")
    
    os.makedirs("CAS_X_v2/results/figures", exist_ok=True)
    sns.set_theme(style="whitegrid")

    # ---------------------------------------------------------
    # FIGURE S2: Independent Validation Summary
    # ---------------------------------------------------------
    val_categories = ['GTEx Expression', 'scRNA-seq\nEvidence', 'Open Targets\nEvidence', 'CRISPR Functional\nScreens']
    val_percentages = [73.3, 53.3, 80.0, 46.7] 
    
    fig_s2, ax_s2 = plt.subplots(figsize=(8, 6))
    sns.barplot(x=val_categories, y=val_percentages, hue=val_categories, palette='mako', legend=False, ax=ax_s2)
    
    ax_s2.set_ylim([0, 100])
    ax_s2.set_title('Independent Validation Summary of CAS-X Targets', fontweight='bold', pad=15)
    ax_s2.set_ylabel('Percentage of Top 15 Targets (%)', fontweight='bold')
    
    for i, v in enumerate(val_percentages):
        ax_s2.text(i, v + 2, f"{v}%", ha='center', fontweight='bold')
        
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/FigS2_Validation_Summary.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated Fig S2 (Validation Summary) at 600 DPI.")

    # ---------------------------------------------------------
    # FIGURE S3: External Validation Evidence Landscape
    # ---------------------------------------------------------
    val_counts = [11, 8, 12, 7] # Counts out of 15
    
    fig_s3, ax_s3 = plt.subplots(figsize=(8, 6))
    sns.barplot(x=val_categories, y=val_counts, hue=val_categories, palette='Blues_d', legend=False, edgecolor='black', ax=ax_s3)
    
    ax_s3.set_ylim([0, 15])
    ax_s3.set_title('External Validation Evidence Landscape', fontweight='bold', pad=15)
    ax_s3.set_ylabel('Count of Supported Targets (out of 15)', fontweight='bold')
    
    for i, v in enumerate(val_counts):
        ax_s3.text(i, v + 0.3, f"{v} / 15", ha='center', fontweight='bold')
        
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/FigS3_Validation_Landscape.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated Fig S3 (Validation Landscape) at 600 DPI.")

    # ---------------------------------------------------------
    # FIGURE S4: Pathway Enrichment Profile (Bubble Plot)
    # ---------------------------------------------------------
    pathways = ['Steroid Hormone\nBiosynthesis', 'Type 2 Diabetes /\nInsulin Signaling', 
                'Cholesterol\nHomeostasis', 'Lipoprotein\nMetabolism']
    p_values_neg_log = [3.05, 3.51, 4.35, 5.92] # -log10(P-value)
    overlap_counts = [2, 3, 3, 3]
    
    fig_s4, ax_s4 = plt.subplots(figsize=(9, 5))
    scatter = ax_s4.scatter(p_values_neg_log, pathways, s=[c * 200 for c in overlap_counts], 
                            c=p_values_neg_log, cmap='viridis', alpha=0.8, edgecolors='black')
    
    ax_s4.set_title('Pathway Enrichment Profile of CAS-X Targets', fontweight='bold', pad=15)
    ax_s4.set_xlabel('Enrichment Significance -log10(P-Value)', fontweight='bold')
    ax_s4.set_xlim([2.5, 6.5])
    
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/FigS4_Pathway_Enrichment.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated Fig S4 (Pathway Enrichment) at 600 DPI.")

    # ---------------------------------------------------------
    # FIGURE S5: Epigenomic State of Lead CAS-X Targets
    # ---------------------------------------------------------
    states = ['Super-Enhancer', 'Active Promoter\n(H3K4me3)', 'Active Enhancer\n(H3K27ac)']
    counts = [2, 6, 5]
    
    fig_s5, ax_s5 = plt.subplots(figsize=(7, 4))
    sns.barplot(x=counts, y=states, hue=states, palette='flare', legend=False, ax=ax_s5)
    
    ax_s5.set_title('Epigenomic State of Lead CAS-X Targets', fontweight='bold', pad=15)
    ax_s5.set_xlabel('Number of Lead eQTL Variants', fontweight='bold')
    ax_s5.set_ylabel('Chromatin State (ENCODE)', fontweight='bold')
    ax_s5.set_xlim([0, 8])
    
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/FigS5_Epigenomic_State.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated Fig S5 (Epigenomic State) at 600 DPI.")

    print("======================================================")
    print(" ALL COMPUTATIONAL ASSETS SUCCESSFULLY GENERATED.")
    print("======================================================")

if __name__ == "__main__":
    generate_remaining_si_figures()
