import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(16, 9))

ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.axis("off")

# -------------------------
# Helper
# -------------------------

def box(x, y, text, color="#dbe9f4"):

    rect = FancyBboxPatch(
        (x, y),
        2.1,
        0.9,
        boxstyle="round,pad=0.03,rounding_size=0.15",
        linewidth=1.5,
        edgecolor="black",
        facecolor=color
    )

    ax.add_patch(rect)

    ax.text(
        x + 1.05,
        y + 0.45,
        text,
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold"
    )

# -------------------------
# Workflow
# -------------------------

box(0.5, 6.5, "GWAS\n50 T2D Loci")
box(3.0, 6.5, "Gene Mapping\n37 Genes")
box(5.5, 6.5, "GTEx eQTL\nIntegration")
box(8.0, 6.5, "Multi-Tissue\nSupport", "#dff2e1")
box(10.5, 6.5, "External\nValidation", "#f6efd0")
box(13.0, 6.5, "Pathway\nEnrichment", "#f3e0f4")

# arrows

for x in [2.6, 5.1, 7.6, 10.1, 12.6]:

    ax.arrow(
        x,
        6.95,
        0.35,
        0,
        width=0.02,
        head_width=0.20,
        head_length=0.20,
        color="black",
        length_includes_head=True
    )

# -------------------------
# Tissue support
# -------------------------

ax.text(
    8.9,
    5.8,
    "Pancreas  |  Adipose  |  Muscle  |  Blood",
    fontsize=13,
    fontweight="bold",
    color="darkgreen",
    ha="center"
)

# -------------------------
# Validation
# -------------------------

ax.text(
    11.6,
    5.0,
    "GTEx (8 genes)\nscRNA-seq (5 genes)\nOpen Targets (8 genes)\nCRISPR (5 genes)",
    ha="center",
    fontsize=11
)

# -------------------------
# Pathways
# -------------------------

ax.text(
    14.1,
    5.0,
    "Type II Diabetes\nInsulin Secretion\nAdipogenesis",
    ha="center",
    fontsize=11
)

# -------------------------
# CAS-X V5 box
# -------------------------

center_box = FancyBboxPatch(
    (5.2, 2.8),
    5.6,
    1.0,
    boxstyle="round,pad=0.04,rounding_size=0.15",
    linewidth=2,
    edgecolor="black",
    facecolor="#f6ebc8"
)

ax.add_patch(center_box)

ax.text(
    8.0,
    3.3,
    "CAS-X V5 PRIORITIZATION",
    ha="center",
    va="center",
    fontsize=18,
    fontweight="bold"
)

# arrow down

ax.arrow(
    8.0,
    6.2,
    0,
    -2.0,
    width=0.02,
    head_width=0.20,
    head_length=0.20,
    color="black"
)

# -------------------------
# Top genes
# -------------------------

ax.text(
    8.0,
    2.0,
    "Top Prioritized Genes",
    ha="center",
    fontsize=14,
    fontweight="bold"
)

ax.text(
    8.0,
    1.45,
    "FTO | JAZF1 | VEGFA | NOTCH2 | TSPAN8",
    ha="center",
    fontsize=15,
    fontweight="bold",
    color="darkred"
)

# -------------------------
# Coverage
# -------------------------

ax.text(
    8.0,
    0.5,
    "+175% Coverage Improvement Through Multi-Tissue Integration",
    ha="center",
    fontsize=13,
    fontweight="bold"
)

plt.title(
    "Graphical Abstract: CAS-X Framework for Prioritizing Type 2 Diabetes Genes",
    fontsize=20,
    fontweight="bold",
    pad=20
)

plt.savefig(
    "results/figures/Graphical_Abstract_Final.png",
    dpi=450,
    bbox_inches="tight"
)

plt.close()

print("Saved Final Graphical Abstract")
