import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure output directory exists
os.makedirs('results/supplementary_figures', exist_ok=True)

# Global Publication Styling (600 DPI, LZW TIFF)
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'Arial', 'font.family': 'sans-serif'})

print("Initializing Complete CAS-X Supplementary 600 DPI TIFF Asset Generation Pipeline (S1 to S8)...")

# ----------------------------------------------------
# Figure S1: Tissue-Specific Regulatory Contribution
# ----------------------------------------------------
print("Generating Supplementary Figure S1...")
fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
tissues = ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Muscle_Skeletal', 'Pancreas', 'Whole_Blood']
counts = [582, 364, 905, 268, 328]
sns.barplot(x=tissues, y=counts, palette='viridis', ax=ax)
ax.set_title("Tissue-Specific Regulatory Contribution to CAS-X Targets", fontsize=12, fontweight='bold', pad=12)
ax.set_ylabel("Count of Significant Signal Hits (|Slope| > 0.1)", fontsize=10, fontweight='bold')
ax.set_xlabel("GTEx Tissue", fontsize=10, fontweight='bold')
plt.xticks(rotation=20, ha='right', fontweight='bold')
plt.tight_layout()
plt.savefig('results/supplementary_figures/Figure_S1_Tissue_Contribution.tiff', format='tiff', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

# ----------------------------------------------------
# Figure S2: Independent Validation Summary (%)
# ----------------------------------------------------
print("Generating Supplementary Figure S2...")
fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
platforms_s2 = ['GTEx / Expression\n(GSE38642)', 'scRNA-seq\n(GSE81608)', 'Open Targets\nT2D Evidence', 'CRISPR Functional\nScreens']
percentages = [73.3, 53.3, 80.0, 46.7]
colors_s2 = ['#36304a', '#445b7f', '#448b94', '#6fc1a1']
bars_s2 = ax.bar(platforms_s2, percentages, color=colors_s2)
ax.set_title("Independent Validation Summary of CAS-X Targets", fontsize=12, fontweight='bold', pad=15)
ax.set_ylabel("Percentage of Top 15 CAS-X Targets (%)", fontsize=10, fontweight='bold')
ax.set_ylim(0, 100)
for bar in bars_s2:
    y = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2.0, y + 2, f"{y}%", ha='center', fontweight='bold', fontsize=10)
plt.tight_layout()
plt.savefig('results/supplementary_figures/Figure_S2_Validation_Summary.tiff', format='tiff', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

# ----------------------------------------------------
# Figure S3: External Validation Evidence Landscape (Count)
# ----------------------------------------------------
print("Generating Supplementary Figure S3...")
fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
platforms_s3 = ['GTEx / Expression\n(GSE38642)', 'scRNA-seq\n(GSE81608)', 'Open Targets\nT2D Evidence', 'CRISPR Functional\nScreens']
counts_val = [11, 8, 12, 7]
colors_s3 = ['#306198', '#5487af', '#92bad1', '#d2e2f0']
bars_s3 = ax.bar(platforms_s3, counts_val, color=colors_s3, edgecolor='black')
ax.set_title("External Validation Evidence Landscape for Top CAS-X Targets", fontsize=11, fontweight='bold', pad=12)
ax.set_ylabel("Count of Supported Target Genes (out of 15)", fontsize=10, fontweight='bold')
ax.set_ylim(0, 15)
for i, v in enumerate(counts_val):
    ax.text(i, v + 0.4, f"{v} / 15", ha='center', fontweight='bold', fontsize=9)
plt.tight_layout()
plt.savefig('results/supplementary_figures/Figure_S3_External_Validation.tiff', format='tiff', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

# ----------------------------------------------------
# Figure S4: Pathway Enrichment Profile (Bubble Plot)
# ----------------------------------------------------
print("Generating Supplementary Figure S4...")
fig, ax = plt.subplots(figsize=(8, 5), dpi=600)
pathways = [
    'FOX0 Signaling Pathway', 
    'Regulation of Lipolysis in Adipocytes', 
    'Adipocytokine Signaling Pathway', 
    'Maturity Onset Diabetes of the Young', 
    'Type 2 Diabetes Mellitus', 
    'Insulin Resistance Pathway'
]
p_values = [2.85, 3.06, 3.49, 4.96, 5.35, 5.92]
sizes = [100, 120, 150, 250, 300, 400]
colors_s4 = sns.color_palette("mako", len(pathways))

ax.scatter(p_values, pathways, s=sizes, c=colors_s4, alpha=0.9)
ax.set_title("Pathway Enrichment Profile of CAS-X Prioritized Targets", fontsize=12, fontweight='bold', pad=15)
ax.set_xlabel("Enrichment Significance -log10(P-Value)", fontsize=10, fontweight='bold')
ax.set_ylabel("Enriched KEGG Biological Pathway", fontsize=10, fontweight='bold')
plt.tight_layout()
plt.savefig('results/supplementary_figures/Figure_S4_Pathway_Enrichment.tiff', format='tiff', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

# ----------------------------------------------------
# Figure S5: Epigenomic State of Lead CAS-X Targets
# ----------------------------------------------------
print("Generating Supplementary Figure S5...")
fig, ax = plt.subplots(figsize=(7, 5), dpi=600)
states = ['Super-Enhancer', 'Active Promoter (H3K4me3)', 'Active Enhancer (H3K27ac)']
state_counts = [2, 4, 9]
sns.barplot(y=states, x=state_counts, palette="flare", ax=ax)
ax.set_title("Epigenomic State of Lead CAS-X Targets", fontsize=12, fontweight='bold', pad=12)
ax.set_xlabel("Number of Lead eQTL Variants", fontsize=11, fontweight='bold')
ax.set_ylabel("Chromatin State (ENCODE)", fontsize=11, fontweight='bold')
ax.set_xlim(0, 10)
plt.tight_layout()
plt.savefig('results/supplementary_figures/Figure_S5_Epigenomic_State.tiff', format='tiff', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

# ----------------------------------------------------
# Figure S6: PCA Scree Plot
# ----------------------------------------------------
print("Generating Supplementary Figure S6...")
fig, ax1 = plt.subplots(figsize=(8, 6), dpi=600)
pcs = ['PC1', 'PC2']
ind_var = [56.3, 43.7]
cum_var = np.cumsum(ind_var)

ax1.bar(pcs, ind_var, color='#5c6b73', alpha=0.9, width=0.6, label='Individual Variance')
ax1.set_ylabel("Percentage of Variance Explained (%)", fontsize=10, color='#2c3e50')
ax1.set_ylim(0, 115)

ax2 = ax1.twinx()
ax2.plot(pcs, cum_var, color='#e0533d', marker='o', linewidth=2.5, markersize=8, label='Cumulative Variance')
ax2.set_ylim(0, 115)
ax2.set_yticks([])

for i, txt in enumerate(cum_var):
    ax2.annotate(f"{txt:.1f}%", (pcs[i], cum_var[i] + 3), ha='center', color='#c0392b', fontweight='bold', fontsize=10)
for i, txt in enumerate(ind_var):
    ax1.annotate(f"{txt:.1f}%", (pcs[i], ind_var[i] - 5), ha='center', color='white', fontweight='bold', fontsize=10)

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines2 + lines, labels2 + labels, loc='lower right')

ax1.set_title("PCA Scree Plot: Multi-Tissue eQTL Variance", fontsize=12, fontweight='bold', pad=12)
ax1.set_xlabel("Principal Components", fontsize=10)
plt.tight_layout()
plt.savefig('results/supplementary_figures/Figure_S6_PCA_Scree.tiff', format='tiff', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

# ----------------------------------------------------
# Figure S7: ROC and PR Curves
# ----------------------------------------------------
print("Generating Supplementary Figure S7...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=600)

fpr = [0.0, 0.17, 0.17, 0.42, 0.42, 0.65, 0.65, 1.0]
tpr = [0.0, 0.25, 0.50, 0.62, 0.75, 0.88, 1.0, 1.0]
ax1.plot(fpr, tpr, color='#2b7bba', linewidth=2.5, label='CAS-X Framework (AUC = 0.72)')
ax1.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random (AUC = 0.50)')
ax1.set_title("A. Receiver Operating Characteristic", fontsize=11, fontweight='bold')
ax1.set_xlabel("False Positive Rate", fontsize=10, fontweight='bold')
ax1.set_ylabel("True Positive Rate", fontsize=10, fontweight='bold')
ax1.set_xlim(0, 1.0)
ax1.set_ylim(0, 1.05)
ax1.legend(loc='lower right')

recall = [0.0, 0.25, 0.50, 0.50, 0.62, 0.62, 0.75, 0.75, 0.88, 0.88, 1.0]
precision = [1.0, 0.28, 0.44, 0.40, 0.45, 0.29, 0.33, 0.30, 0.33, 0.27, 0.30]
ax2.plot(recall, precision, color='#b22222', linewidth=2.5, label='CAS-X Framework (AUC = 0.35)')
ax2.axhline(y=0.22, linestyle='--', color='gray', label='Baseline (AUC = 0.22)')
ax2.set_title("B. Precision-Recall Curve", fontsize=11, fontweight='bold')
ax2.set_xlabel("Recall (Sensitivity)", fontsize=10, fontweight='bold')
ax2.set_ylabel("Precision (Positive Predictive Value)", fontsize=10, fontweight='bold')
ax2.set_xlim(0, 1.0)
ax2.set_ylim(0, 1.05)
ax2.legend(loc='upper right')

plt.suptitle("Quantitative Validation Metrics of the CAS-X Framework", fontsize=13, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig('results/supplementary_figures/Figure_S7_ROC_PR.tiff', format='tiff', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

# ----------------------------------------------------
# Figure S8: Benchmark Comparison
# ----------------------------------------------------
print("Generating Supplementary Figure S8...")
fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
models = ['Random Baseline', 'Standard Single-Tissue\n(Pancreas TWAS Proxy)', 'CAS-X Framework\n(Multi-Tissue PCA)']
auroc_scores = [0.50, 0.82, 0.72]
bars = ax.bar(models, auroc_scores, color=['#7f8c8d', '#f39c12', '#2980b9'], width=0.6, edgecolor='black')
ax.axhline(y=0.50, linestyle='--', color='gray')
ax.set_title("Performance Benchmark: CAS-X vs Single-Tissue Approach", fontsize=12, fontweight='bold', pad=15)
ax.set_ylabel("Predictive Performance (AUROC)", fontsize=10, fontweight='bold')
ax.set_ylim(0, 1.0)
for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.03, f"{yval:.2f}", ha='center', fontweight='bold', fontsize=11)

plt.tight_layout()
plt.savefig('results/supplementary_figures/Figure_S8_Benchmark.tiff', format='tiff', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

print("All 8 Supplementary TIFF figures successfully generated at 600 DPI in results/supplementary_figures/")
