import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Setup Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PUB_DIR = os.path.join(PROJECT_ROOT, "results", "v6_publication_assets")
os.makedirs(PUB_DIR, exist_ok=True)

print("======================================================")
print("  GENERATING POLISHED FIGURE 1 WORKFLOW")
print("======================================================\n")

fig, ax = plt.subplots(figsize=(10, 12))
ax.axis('off')
ax.set_xlim(0, 10)
ax.set_ylim(0, 12)

# Define step parameters
steps = [
    ("1. Genomic Baseline", "50 Harmonized T2D GWAS Loci & Candidate Gene Mapping", '#34495E'),
    ("2. Transcriptomic Extraction", "Continuous GTEx eQTL Effect Sizes Across 5 Metabolic Tissues", '#2E86C1'),
    ("3. Mathematical Core", "Unsupervised PCA Scoring (CAS-X) without Heuristic Bias", '#D35400'),
    ("4. Robustness Validations", "LOTO Sensitivity & Non-Metabolic Negative Control Benchmarking", '#27AE60'),
    ("5. Biological & Clinical Mapping", "PheWAS Pleiotropy, Epigenomic Architecture, & Clinical Tractability", '#8E44AD'),
    ("6. Final Output", "Systemic 0-100 Prioritization Ranking for T2D Targets", '#2C3E50')
]

y_positions = [10.5, 8.7, 6.9, 5.1, 3.3, 1.5]

for i, (title, desc, color) in enumerate(steps):
    y = y_positions[i]
    
    # Draw main box
    box = mpatches.FancyBboxPatch((2, y-0.5), 6, 1, boxstyle="round,pad=0.2,rounding_size=0.2", 
                                  linewidth=2, edgecolor=color, facecolor='#F8F9F9')
    ax.add_patch(box)
    
    # Draw color tab on the left of the box
    tab = mpatches.FancyBboxPatch((2, y-0.5), 0.2, 1, boxstyle="round,pad=0.2,rounding_size=0.2", 
                                  linewidth=0, facecolor=color)
    ax.add_patch(tab)
    
    # Add text
    ax.text(2.5, y+0.15, title, fontsize=12, fontweight='bold', color=color, ha='left')
    ax.text(2.5, y-0.2, desc, fontsize=10, color='black', ha='left', wrap=True)
    
    # Add connecting arrows (except for the last box)
    if i < len(steps) - 1:
        ax.annotate('', xy=(5, y-0.6), xytext=(5, y-1.1),
                    arrowprops=dict(facecolor='#7F8C8D', edgecolor='#7F8C8D', width=3, headwidth=10))

plt.title("Figure 1: CAS-X Computational & Validation Workflow", fontsize=16, fontweight='bold', y=0.95)

fig1_path = os.path.join(PUB_DIR, "Figure1_Polished_Workflow_v6.png")
plt.savefig(fig1_path, dpi=500, bbox_inches='tight')
plt.close()

print(f"SUCCESS: Saved highly polished workflow diagram to: {fig1_path}")
print("======================================================\n")
