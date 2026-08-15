import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Setup Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
FIGURES_DIR = os.path.join(PROJECT_ROOT, "results", "figures")
TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "manuscript_tables")
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

# Create the pristine publication subfolder
PUB_DIR = os.path.join(PROJECT_ROOT, "results", "v6_publication_assets")
os.makedirs(PUB_DIR, exist_ok=True)

print("======================================================")
print("  CONSOLIDATING V6 ASSETS FOR PUBLICATION")
print("======================================================\n")

# 1. Copy existing V6 assets to the new folder
v6_files = [
    (FIGURES_DIR, "Figure3_Multitissue_Heatmap_v6.png"),
    (FIGURES_DIR, "Figure4_Top15_CASX_v6.png"),
    (FIGURES_DIR, "Figure5_Tissue_Contribution_v6.png"),
    (FIGURES_DIR, "Figure8_Pathway_Enrichment_v6.png"),
    (FIGURES_DIR, "Negative_Control_Benchmark.png"),
    (FIGURES_DIR, "LOTO_Sensitivity_Analysis.png"),
    (TABLES_DIR, "Table2_Top15_CASX_v6.csv"),
    (TABLES_DIR, "Table4_Pathway_Enrichment_v6.csv")
]

copied_count = 0
for src_dir, file_name in v6_files:
    src_path = os.path.join(src_dir, file_name)
    if os.path.exists(src_path):
        shutil.copy(src_path, os.path.join(PUB_DIR, file_name))
        copied_count += 1
print(f"Successfully moved {copied_count} finalized V6 assets into /v6_publication_assets/")

# 2. Generate Updated Figure 2: Coverage Improvement
print("Recalculating Multi-Tissue Coverage for Figure 2...")
gtex_file = os.path.join(PROCESSED_DIR, "casx_v6_continuous_gtex.csv")
gtex_df = pd.read_csv(gtex_file)

# Filter for significant continuous hits (|slope| > 0.1)
sig_hits = gtex_df[gtex_df['slope'].abs() > 0.1]
pancreas_genes = sig_hits[sig_hits['tissue'] == 'Pancreas']['gene_symbol'].nunique()
multi_genes = sig_hits['gene_symbol'].nunique()
pct_increase = ((multi_genes - pancreas_genes) / pancreas_genes) * 100

cov_data = pd.DataFrame({
    'Approach': ['Single-Tissue\n(Pancreas Only)', 'CAS-X Multi-Tissue\n(5 GTEx Tissues)'],
    'Supported Target Genes': [pancreas_genes, multi_genes]
})

plt.figure(figsize=(8, 6))
sns.set_theme(style="whitegrid")
ax = sns.barplot(data=cov_data, x='Approach', y='Supported Target Genes', palette=['#4C72B0', '#55A868'])

plt.title('Figure 2: Improvement in Target Coverage via Multi-Tissue Integration', fontsize=14, fontweight='bold', pad=20)
plt.ylabel('Number of Target Genes with eQTL Support', fontsize=12, fontweight='bold')
plt.xlabel('', fontsize=12)
plt.ylim(0, multi_genes + 5)
for p in ax.patches:
    ax.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='bottom', fontsize=12, fontweight='bold', xytext=(0, 5), textcoords='offset points')

plt.text(0.5, multi_genes + 2, f"+{pct_increase:.0f}% Coverage Increase", ha='center', fontsize=14, fontweight='bold', color='#D81B60')

fig2_path = os.path.join(PUB_DIR, "Figure2_Coverage_Improvement_v6.png")
plt.savefig(fig2_path, dpi=500, bbox_inches='tight')
plt.close()
print("Generated updated Figure 2 (Continuous Coverage Improvement).")

# 3. Note regarding Figure 1, 6, and 7
print("\nPENDING UPDATES:")
print("- Figure 1 (Workflow Diagram): Needs manual diagram update to reflect PCA/Continuous logic.")
print("- Figures 6 & 7 (External Validation): Awaiting raw OpenTargets/CRISPR overlap data.")
print("\nConsolidation complete! Check the 'results/v6_publication_assets/' folder.")
