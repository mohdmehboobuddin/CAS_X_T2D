import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(16, 9))

ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.axis("off")

# ==========================
# TITLE
# ==========================

ax.text(
    8,
    9.4,
    "CAS-X: Multi-Tissue Prioritization Framework for Type 2 Diabetes Genes",
    ha="center",
    fontsize=22,
    fontweight="bold"
)

# ==========================
# ACT 1
# ==========================

act1 = FancyBboxPatch(
    (0.5, 4.8),
    4.2,
    2.4,
    boxstyle="round,pad=0.05",
    linewidth=2,
    facecolor="#EAF2FB"
)

ax.add_patch(act1)

ax.text(
    2.6,
    6.7,
    "DATA INTEGRATION",
    ha="center",
    fontsize=18,
    fontweight="bold",
    color="#1B4F72"
)

ax.text(
    2.6,
    5.8,
    "GWAS\n50 T2D Loci",
    ha="center",
    fontsize=15,
    fontweight="bold"
)

ax.text(
    2.6,
    4.9,
    "GTEx eQTL\n+\nPancreas | Adipose | Muscle | Blood",
    ha="center",
    fontsize=13
)

# ==========================
# ARROW 1
# ==========================

ax.arrow(
    4.9,
    6.0,
    1.3,
    0,
    width=0.03,
    head_width=0.30,
    head_length=0.25,
    color="black"
)

# ==========================
# ACT 2
# ==========================

act2 = FancyBboxPatch(
    (6.2, 4.5),
    3.6,
    3.0,
    boxstyle="round,pad=0.05",
    linewidth=2,
    facecolor="#FFF4D6"
)

ax.add_patch(act2)

ax.text(
    8.0,
    6.8,
    "CAS-X V5",
    ha="center",
    fontsize=22,
    fontweight="bold",
    color="#7D6608"
)

ax.text(
    8.0,
    5.8,
    "Gene Mapping\n+\nEvidence Integration\n+\nPrioritization",
    ha="center",
    fontsize=14
)

# ==========================
# ARROW 2
# ==========================

ax.arrow(
    10.0,
    6.0,
    1.3,
    0,
    width=0.03,
    head_width=0.30,
    head_length=0.25,
    color="black"
)

# ==========================
# ACT 3
# ==========================

act3 = FancyBboxPatch(
    (11.2, 4.5),
    4.1,
    3.0,
    boxstyle="round,pad=0.05",
    linewidth=2,
    facecolor="#E8F8F5"
)

ax.add_patch(act3)

ax.text(
    13.25,
    6.8,
    "VALIDATION",
    ha="center",
    fontsize=18,
    fontweight="bold",
    color="#117864"
)

ax.text(
    13.25,
    5.9,
    "✓ Open Targets\n✓ CRISPR\n✓ scRNA-seq\n✓ GTEx",
    ha="center",
    fontsize=14
)

# ==========================
# TOP GENES
# ==========================

gene_box = FancyBboxPatch(
    (4.5, 1.4),
    7.0,
    1.3,
    boxstyle="round,pad=0.05",
    linewidth=2,
    facecolor="#FDEDEC"
)

ax.add_patch(gene_box)

ax.text(
    8,
    2.2,
    "TOP PRIORITIZED GENES",
    ha="center",
    fontsize=16,
    fontweight="bold"
)

ax.text(
    8,
    1.75,
    "FTO   |   JAZF1   |   VEGFA   |   NOTCH2   |   TSPAN8",
    ha="center",
    fontsize=15,
    fontweight="bold",
    color="#922B21"
)

# ==========================
# PATHWAYS
# ==========================

ax.text(
    8,
    0.7,
    "Type II Diabetes • Insulin Secretion • Adipogenesis",
    ha="center",
    fontsize=14,
    fontweight="bold",
    color="#145A32"
)

ax.text(
    8,
    0.25,
    "+175% Coverage Improvement Through Multi-Tissue Integration",
    ha="center",
    fontsize=13,
    fontweight="bold"
)

plt.savefig(
    "results/figures/Graphical_Abstract_Final_V2.png",
    dpi=450,
    bbox_inches="tight"
)

plt.close()

print("Saved Graphical_Abstract_Final_V2")
