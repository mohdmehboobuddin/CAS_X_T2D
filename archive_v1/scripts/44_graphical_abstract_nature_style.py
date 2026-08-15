import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# --------------------------------------------------
# Figure setup
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(18, 10))

ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# --------------------------------------------------
# Title
# --------------------------------------------------

plt.title(
    "CAS-X: Multi-Tissue Prioritization Framework for Type 2 Diabetes Genes",
    fontsize=28,
    fontweight="bold",
    pad=20
)

# --------------------------------------------------
# Helper function
# --------------------------------------------------

def add_box(x, y, w, h, color, text, fontsize=18):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=1",
        linewidth=2.5,
        edgecolor="black",
        facecolor=color
    )

    ax.add_patch(box)

    ax.text(
        x + w/2,
        y + h/2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight="bold"
    )

# --------------------------------------------------
# Top workflow
# --------------------------------------------------

add_box(
    3, 62, 18, 14,
    "#DCE6F2",
    "GWAS\n\n50 T2D Loci",
    18
)

add_box(
    27, 62, 22, 14,
    "#D9EAD3",
    "Multi-Tissue Integration\n\nPancreas • Adipose\nMuscle • Blood",
    16
)

add_box(
    55, 58, 18, 22,
    "#FFF2CC",
    "CAS-X V5\n\nEvidence\nIntegration\n+\nGene Prioritization",
    18
)

add_box(
    79, 58, 18, 22,
    "#D9EAD3",
    "Validation\n\n"
    "GTEx (8 genes)\n"
    "scRNA-seq (5 genes)\n"
    "Open Targets (8 genes)\n"
    "CRISPR (5 genes)",
    13
)

# --------------------------------------------------
# Horizontal arrows
# --------------------------------------------------

for start_x, end_x in [(21,27), (49,55), (73,79)]:

    arrow = FancyArrowPatch(
        (start_x,69),
        (end_x,69),
        arrowstyle="simple",
        mutation_scale=35,
        color="black"
    )

    ax.add_patch(arrow)

# --------------------------------------------------
# Vertical arrow
# --------------------------------------------------

ax.add_patch(
    FancyArrowPatch(
        (64,58),
        (64,44),
        arrowstyle="simple",
        mutation_scale=35,
        color="black"
    )
)

# --------------------------------------------------
# Top prioritized genes
# --------------------------------------------------

gene_box = FancyBboxPatch(
    (25,24),
    50,
    12,
    boxstyle="round,pad=0.02,rounding_size=1",
    linewidth=2.5,
    edgecolor="black",
    facecolor="#F4E1E1"
)

ax.add_patch(gene_box)

ax.text(
    50,
    31,
    "TOP PRIORITIZED GENES",
    ha="center",
    va="center",
    fontsize=20,
    fontweight="bold"
)

ax.text(
    50,
    26,
    "FTO   |   JAZF1   |   VEGFA   |   NOTCH2   |   TSPAN8",
    ha="center",
    va="center",
    fontsize=18,
    color="darkred",
    fontweight="bold"
)

# --------------------------------------------------
# Pathway enrichment
# --------------------------------------------------

ax.text(
    50,
    16,
    "Type II Diabetes  •  Insulin Secretion  •  Adipogenesis",
    ha="center",
    fontsize=18,
    color="darkgreen",
    fontweight="bold"
)

# --------------------------------------------------
# Main finding
# --------------------------------------------------

ax.text(
    50,
    8,
    "+175% Coverage Improvement Through Multi-Tissue Integration",
    ha="center",
    fontsize=20,
    fontweight="bold"
)

# --------------------------------------------------
# Save
# --------------------------------------------------

plt.tight_layout()

plt.savefig(
    "results/figures/Graphical_Abstract_Final.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Graphical Abstract")
