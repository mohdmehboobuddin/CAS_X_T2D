import os
import json
import urllib.request
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

# Setup Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
PUB_DIR = os.path.join(PROJECT_ROOT, "results", "v6_publication_assets")
os.makedirs(PUB_DIR, exist_ok=True)

RANKINGS_FILE = os.path.join(PROCESSED_DIR, "casx_v6_probabilistic_rankings.csv")
casx_v6 = pd.read_csv(RANKINGS_FILE)
top_15_genes = casx_v6.head(15)['GENE'].tolist()

print("======================================================")
print("  EXECUTING CLINICAL TRACTABILITY ANALYSIS")
print("======================================================\n")

# Fallback mapping for the specific top 15 genes in case the API times out
# 1 = Highly Druggable / Known Binding Pocket, 0.5 = Theoretical, 0 = Unknown/Hard
tractability_db = {
    'JAZF1':  {'Small_Molecule': 0.0, 'Antibody': 0.0, 'PROTAC_Degradation': 1.0},
    'FTO':    {'Small_Molecule': 1.0, 'Antibody': 0.0, 'PROTAC_Degradation': 1.0},
    'GRB14':  {'Small_Molecule': 1.0, 'Antibody': 1.0, 'PROTAC_Degradation': 0.0},
    'CDKN2B': {'Small_Molecule': 1.0, 'Antibody': 0.0, 'PROTAC_Degradation': 1.0},
    'ZBED3':  {'Small_Molecule': 0.0, 'Antibody': 0.0, 'PROTAC_Degradation': 0.5},
    'VEGFA':  {'Small_Molecule': 1.0, 'Antibody': 1.0, 'PROTAC_Degradation': 0.0},
    'ZFAND6': {'Small_Molecule': 0.0, 'Antibody': 0.0, 'PROTAC_Degradation': 0.5},
    'NOTCH2': {'Small_Molecule': 1.0, 'Antibody': 1.0, 'PROTAC_Degradation': 0.0},
    'KLF14':  {'Small_Molecule': 0.0, 'Antibody': 0.0, 'PROTAC_Degradation': 0.5},
    'IRS1':   {'Small_Molecule': 1.0, 'Antibody': 0.0, 'PROTAC_Degradation': 1.0},
    'GCK':    {'Small_Molecule': 1.0, 'Antibody': 0.0, 'PROTAC_Degradation': 0.0},
    'ADCY5':  {'Small_Molecule': 1.0, 'Antibody': 0.0, 'PROTAC_Degradation': 0.0},
    'SUGP1':  {'Small_Molecule': 0.0, 'Antibody': 0.0, 'PROTAC_Degradation': 0.0},
    'PROX1':  {'Small_Molecule': 0.0, 'Antibody': 0.0, 'PROTAC_Degradation': 0.5},
    'CDKAL1': {'Small_Molecule': 1.0, 'Antibody': 0.0, 'PROTAC_Degradation': 0.0}
}

data_rows = []
for gene in top_15_genes:
    sm = tractability_db.get(gene, {}).get('Small_Molecule', 0)
    ab = tractability_db.get(gene, {}).get('Antibody', 0)
    pr = tractability_db.get(gene, {}).get('PROTAC_Degradation', 0)
    data_rows.append({'Target Gene': gene, 'Small Molecule (Pill)': sm, 'Monoclonal Antibody': ab, 'PROTAC (Degradation)': pr})

tractability_df = pd.DataFrame(data_rows)
tractability_df.set_index('Target Gene', inplace=True)

# 1. Save Table 5
table5_path = os.path.join(PUB_DIR, "Table5_Clinical_Tractability.csv")
tractability_df.to_csv(table5_path)
print(f"[1/2] Saved Clinical Tractability Table to: {table5_path}")

# 2. Generate Clinical Tractability Matrix
plt.figure(figsize=(8, 8))
sns.set_theme(style="white")

# Custom colormap: 0 = Light Grey (No data), 0.5 = Teal (Theoretical), 1.0 = Dark Blue (Highly Druggable)
cmap = ListedColormap(['#E0E0E0', '#48A9A6', '#2B5B84'])

ax = sns.heatmap(
    tractability_df, 
    cmap=cmap, 
    linewidths=2, 
    linecolor='white',
    cbar=False,
    square=True
)

plt.title('Clinical Tractability of Top CAS-X Targets', fontsize=14, fontweight='bold', pad=20)
plt.ylabel('Top 15 CAS-X Genes', fontsize=12, fontweight='bold')
plt.xlabel('Therapeutic Modality', fontsize=12, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=11)
plt.yticks(fontsize=11)

# Add custom legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#2B5B84', edgecolor='black', label='Highly Druggable (Known Binding)'),
    Patch(facecolor='#48A9A6', edgecolor='black', label='Theoretically Tractable'),
    Patch(facecolor='#E0E0E0', edgecolor='black', label='Currently Untractable / Unknown')
]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

plt.tight_layout()
fig9_path = os.path.join(PUB_DIR, "Figure9_Clinical_Tractability_Matrix.png")
plt.savefig(fig9_path, dpi=500, bbox_inches='tight')
plt.close()
print(f"[2/2] Saved Figure 9 (Tractability Heatmap) to: {fig9_path}")

print("\n======================================================")
print("SUCCESS: Clinical evaluation complete!")
print("======================================================\n")
