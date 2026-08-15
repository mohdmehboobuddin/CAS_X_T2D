import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(16, 8))

ax.set_xlim(0, 16)
ax.set_ylim(0, 10)

ax.axis("off")

# -----------------------------
# Main workflow boxes
# -----------------------------

workflow = [

    ("GWAS\n50 T2D Loci", 1.5, 7.5, "#D6EAF8"),

    ("Gene Mapping\n37 Genes", 4.0, 7.5, "#D6EAF8"),

    ("GTEx eQTL\nIntegration", 6.5, 7.5, "#D6EAF8"),

    ("Multi-Tissue\nSupport", 9.0, 7.5, "#D5F5E3"),

    ("Expression\nValidation", 11.5, 7.5, "#FCF3CF"),

    ("External\nValidation", 14.0, 7.5, "#FADBD8")

]

for label, x, y, color in workflow:

    ax.text(
        x,
        y,
        label,
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.6",
            facecolor=color,
            edgecolor="black",
            linewidth=1.5
        )
    )

# arrows

for x in [2.3, 4.8, 7.3, 9.8, 12.3]:

    ax.arrow(
        x,
        7.5,
        0.9,
        0,
        length_includes_head=True,
        head_width=0.15,
        head_length=0.20,
        linewidth=2,
        color="black"
    )

# -----------------------------
# Tissue layer
# -----------------------------

ax.text(
    9,
    5.5,
    "Pancreas   |   Adipose   |   Muscle   |   Blood",
    ha="center",
    fontsize=12,
    fontweight="bold",
    color="darkgreen"
)

# -----------------------------
# Validation layer
# -----------------------------

ax.text(
    11.5,
    4.5,
    "GSE38642\nGTEx Validation",
    ha="center",
    fontsize=10
)

ax.text(
    11.5,
    3.5,
    "GSE81608\nscRNA Validation",
    ha="center",
    fontsize=10
)

ax.text(
    14,
    4.5,
    "Open Targets\n8 Genes",
    ha="center",
    fontsize=10
)

ax.text(
    14,
    3.5,
    "CRISPR Support\n5 Genes",
    ha="center",
    fontsize=10
)

# -----------------------------
# Pathway layer
# -----------------------------

ax.text(
    8,
    2.0,
    "Pathway Enrichment",
    fontsize=14,
    fontweight="bold",
    ha="center"
)

ax.text(
    8,
    1.2,
    "Insulin Signaling   |   Type II Diabetes   |   PPAR Signaling",
    fontsize=11,
    ha="center"
)

# -----------------------------
# Final CAS-X output
# -----------------------------

ax.text(
    8,
    9.2,
    "CAS-X V5 PRIORITIZATION",
    ha="center",
    fontsize=18,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.6",
        facecolor="#FFF2CC",
        edgecolor="black",
        linewidth=2
    )
)

ax.text(
    8,
    8.4,
    "Top Prioritized Genes",
    ha="center",
    fontsize=12,
    fontweight="bold"
)

ax.text(
    8,
    7.9,
    "FTO  |  JAZF1  |  VEGFA  |  NOTCH2  |  TSPAN8",
    ha="center",
    fontsize=12,
    fontweight="bold",
    color="darkred"
)

# -----------------------------
# Coverage improvement
# -----------------------------

ax.text(
    8,
    0.2,
    "+175% Coverage Improvement Through Multi-Tissue Integration",
    ha="center",
    fontsize=12,
    fontweight="bold"
)

plt.title(
    "Graphical Abstract: CAS-X Framework for Prioritizing Type 2 Diabetes Genes",
    fontsize=18,
    pad=25
)

plt.tight_layout()

plt.savefig(
    "results/figures/Graphical_Abstract_CASX_V2.png",
    dpi=450,
    bbox_inches="tight"
)

plt.savefig(
    "results/figures/Graphical_Abstract_CASX_V2.pdf",
    bbox_inches="tight"
)

print("Saved Graphical Abstract V2")
