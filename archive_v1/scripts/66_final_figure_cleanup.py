import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

# Setup Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
PUB_DIR = os.path.join(PROJECT_ROOT, "results", "v6_publication_assets")
os.makedirs(PUB_DIR, exist_ok=True)

RANKINGS_FILE = os.path.join(PROCESSED_DIR, "casx_v6_probabilistic_rankings.csv")
GTEX_FILE = os.path.join(PROCESSED_DIR, "casx_v6_continuous_gtex.csv")

casx_v6 = pd.read_csv(RANKINGS_FILE)
gtex_df = pd.read_csv(GTEX_FILE)
top_15_genes = casx_v6.head(15)['GENE'].tolist()

print("======================================================")
print("  REFORMATTING FIGURES (REMOVING HARDCODED NUMBERS)")
print("======================================================\n")

# ----------------------------------------------------------------------
# 1. FIGURE 1: Minimalist Flowchart (Matching Reference Style)
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 10))
ax.axis('off')
ax.set_xlim(0, 10)
ax.set_ylim(0, 12)

steps = [
    "T2D GWAS Loci Integration",
    "Candidate Gene Mapping",
    "Multi-Tissue GTEx eQTLs",
    "Unsupervised PCA (CAS-X)",
    "Biological & Clinical Validation",
    "Target Prioritization Ranking"
]

y_pos = 10.5
for i, step in enumerate(steps):
    # White box, black border, rounded corners
    box = mpatches.FancyBboxPatch((2.5, y_pos-0.4), 5, 0.8, boxstyle="round,pad=0.3,rounding_size=0.2", 
                                  linewidth=1.5, edgecolor='black', facecolor='white')
    ax.add_patch(box)
    
    # Text centered
    ax.text(5, y_pos, step, fontsize=11, color='black', ha='center', va='center')
    
    # Blue Down Arrow
    if i < len(steps) - 1:
        ax.annotate('', xy=(5, y_pos-0.9), xytext=(5, y_pos-0.4),
                    arrowprops=dict(facecolor='#0072B2', edgecolor='black', width=1.5, headwidth=7))
    
    y_pos -= 1.8

plt.title("CAS-X Framework", fontsize=16, y=0.95)
fig1_path = os.path.join(PUB_DIR, "Figure1_CASX_Framework_Clean.png")
plt.savefig(fig1_path, dpi=500, bbox_inches='tight')
plt.close()
print("[1/5] Regenerated Figure 1 (Minimalist Flowchart).")

# ----------------------------------------------------------------------
# 2. FIGURE 2: Coverage Improvement (Title Fix)
# ----------------------------------------------------------------------
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
ax = sns.barplot(data=cov_data, x='Approach', y='Supported Target Genes', hue='Approach', legend=False, palette=['#4C72B0', '#55A868'])

# CLEAN TITLE
plt.title('Improvement in Target Coverage via Multi-Tissue Integration', fontsize=14, fontweight='bold', pad=20)
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
print("[2/5] Regenerated Figure 2 (Title fixed).")

# ----------------------------------------------------------------------
# 3. FIGURE 9: Clinical Tractability (Title Fix)
# ----------------------------------------------------------------------
tractability_db = {
    'JAZF1':  {'Small Molecule (Pill)': 0.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 1.0},
    'FTO':    {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 1.0},
    'GRB14':  {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 1.0, 'PROTAC (Degradation)': 0.0},
    'CDKN2B': {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 1.0},
    'ZBED3':  {'Small Molecule (Pill)': 0.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 0.5},
    'VEGFA':  {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 1.0, 'PROTAC (Degradation)': 0.0},
    'ZFAND6': {'Small Molecule (Pill)': 0.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 0.5},
    'NOTCH2': {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 1.0, 'PROTAC (Degradation)': 0.0},
    'KLF14':  {'Small Molecule (Pill)': 0.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 0.5},
    'IRS1':   {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 1.0},
    'GCK':    {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 0.0},
    'ADCY5':  {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 0.0},
    'SUGP1':  {'Small Molecule (Pill)': 0.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 0.0},
    'PROX1':  {'Small Molecule (Pill)': 0.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 0.5},
    'CDKAL1': {'Small Molecule (Pill)': 1.0, 'Monoclonal Antibody': 0.0, 'PROTAC (Degradation)': 0.0}
}
tractability_df = pd.DataFrame([{'Target Gene': k, **v} for k, v in tractability_db.items()]).set_index('Target Gene')

plt.figure(figsize=(8, 8))
sns.set_theme(style="white")
cmap = ListedColormap(['#E0E0E0', '#48A9A6', '#2B5B84'])

ax = sns.heatmap(tractability_df, cmap=cmap, linewidths=2, linecolor='white', cbar=False, square=True)

# CLEAN TITLE
plt.title('Clinical Tractability of Top CAS-X Targets', fontsize=14, fontweight='bold', pad=20)
plt.ylabel('Top 15 CAS-X Genes', fontsize=12, fontweight='bold')
plt.xlabel('Therapeutic Modality', fontsize=12, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=11)
plt.yticks(fontsize=11)

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#2B5B84', edgecolor='black', label='Highly Druggable (Known Binding)'),
    Patch(facecolor='#48A9A6', edgecolor='black', label='Theoretically Tractable'),
    Patch(facecolor='#E0E0E0', edgecolor='black', label='Currently Untractable / Unknown')
]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

fig9_path = os.path.join(PUB_DIR, "Figure9_Clinical_Tractability_Matrix.png")
plt.savefig(fig9_path, dpi=500, bbox_inches='tight')
plt.close()
print("[3/5] Regenerated Figure 9 (Title fixed).")

# ----------------------------------------------------------------------
# 4. FIGURE 10: Epigenomic State (Title Fix)
# ----------------------------------------------------------------------
reg_df = pd.read_csv(os.path.join(PUB_DIR, "Table6_Regulatory_Architecture.csv"))

plt.figure(figsize=(7, 5))
sns.set_theme(style="ticks")
ax = sns.countplot(data=reg_df, y='Regulatory Architecture', hue='Regulatory Architecture', legend=False, palette='magma', order=reg_df['Regulatory Architecture'].value_counts().index)

# CLEAN TITLE
plt.title('Epigenomic State of Lead CAS-X Variants', fontsize=12, fontweight='bold', pad=15)
plt.xlabel('Number of Lead eQTL Variants', fontsize=11, fontweight='bold')
plt.ylabel('Chromatin State (ENCODE)', fontsize=11, fontweight='bold')
sns.despine()
plt.tight_layout()

fig10_path = os.path.join(PUB_DIR, "Figure10_Regulatory_Distribution.png")
plt.savefig(fig10_path, dpi=500, bbox_inches='tight')
plt.close()
print("[4/5] Regenerated Figure 10 (Title fixed).")

# ----------------------------------------------------------------------
# 5. FIGURE 11: Pleiotropy Network (Title Fix)
# ----------------------------------------------------------------------
pleiotropy_data = {
    'FTO': ['BMI / Obesity', 'Fasting Insulin', 'Triglycerides'],
    'IRS1': ['Fasting Insulin', 'HDL Cholesterol', 'Coronary Artery Disease'],
    'JAZF1': ['Fasting Glucose', 'HbA1c', 'Prostate Cancer'],
    'CDKN2B': ['Coronary Artery Disease', 'Fasting Glucose'],
    'VEGFA': ['Diabetic Retinopathy', 'BMI / Obesity'],
    'GRB14': ['Fasting Insulin', 'Waist-to-Hip Ratio'],
    'NOTCH2': ['Bone Mineral Density', 'BMI / Obesity'],
    'KLF14': ['HDL Cholesterol', 'Triglycerides', 'Waist-to-Hip Ratio'],
    'GCK': ['Fasting Glucose', 'Birth Weight'],
    'ADCY5': ['Fasting Glucose', 'Birth Weight']
}

G = nx.Graph()
for gene in top_15_genes:
    G.add_node(gene, type='gene')
    traits = pleiotropy_data.get(gene, ['Type 2 Diabetes (Primary)'])
    for trait in traits:
        G.add_node(trait, type='trait')
        G.add_edge(gene, trait)

plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G, seed=42, k=0.6)

genes = [node for node, attr in G.nodes(data=True) if attr['type'] == 'gene']
traits = [node for node, attr in G.nodes(data=True) if attr['type'] == 'trait']

nx.draw_networkx_nodes(G, pos, nodelist=genes, node_color='#2B5B84', node_size=600, alpha=0.9, edgecolors='black')
nx.draw_networkx_nodes(G, pos, nodelist=traits, node_color='#D81B60', node_size=800, node_shape='s', alpha=0.8, edgecolors='black')
nx.draw_networkx_edges(G, pos, width=1.5, alpha=0.5, edge_color='gray')
nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', font_color='black')

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='CAS-X Target Gene', markerfacecolor='#2B5B84', markersize=12, markeredgecolor='black'),
    Line2D([0], [0], marker='s', color='w', label='Systemic Metabolic Trait', markerfacecolor='#D81B60', markersize=12, markeredgecolor='black')
]
plt.legend(handles=legend_elements, loc='upper right', frameon=True)

# CLEAN TITLE
plt.title("Phenome-Wide Pleiotropy (PheWAS) Network of CAS-X Targets", fontsize=14, fontweight='bold', pad=20)
plt.axis('off')
plt.tight_layout()

fig11_path = os.path.join(PUB_DIR, "Figure11_Pleiotropy_Network.png")
plt.savefig(fig11_path, dpi=500, bbox_inches='tight')
plt.close()
print("[5/5] Regenerated Figure 11 (Title fixed).")

print("\n======================================================")
print("SUCCESS: All selected figures are now formatted for publication!")
print("======================================================\n")
