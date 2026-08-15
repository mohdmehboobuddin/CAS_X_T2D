import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(16, 8))

ax.set_xlim(0, 16)
ax.set_ylim(0, 10)

ax.axis("off")

# ==================================================
# MAIN WORKFLOW
# ==================================================

boxes = [

    ("GWAS\n50 T2D Loci", 1.2, "#D6EAF8"),

    ("Gene Mapping\n37 Genes", 3.8, "#D6EAF8"),

    ("GTEx eQTL\nIntegration", 6.4, "#D6EAF8"),

    ("Multi-Tissue\nSupport", 9.0, "#D5F5E3"),

    ("Expression\nValidation", 11.6, "#FCF3CF"),

    ("External\nValidation", 14.2, "#FADBD8")

]

for label, x, color in boxes:

    ax.text(
        x,
        8,
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

for x in [2.0, 4.6, 7.2, 9.8, 12.4]:

    ax.arrow(
        x,
        8,
        0.9,
        0,
        length_includes_head=True,
        head_width=0.15,
        head_length=0.20,
        linewidth=2,
        color="black"
    )

# ==================================================
# SUPPORTING EVIDENCE
# ==================================================

ax.text(
    9.0,
    6.7,
    "Pancreas | Adipose | Muscle | Blood",
    ha="center",
    fontsize=11,
    color="darkgreen",
    fontweight="bold"
)

ax.text(
    11.6,
    6.8,
    "GSE38642",
    ha="center",
    fontsize=10
)

ax.text(
    11.6,
    6.3,
    "GSE81608",
    ha="center",
    fontsize=10
)

ax.text(
    14.2,
    6.8,
    "Open Targets\n8 Genes",
    ha="center",
    fontsize=10
)

ax.text(
    14.2,
    6.0,
    "CRISPR\n5 Genes",
    ha="center",
    fontsize=10
)

# ==================================================
# PATHWAY ENRICHMENT
# ==================================================

ax.text(
    8,
    4.6,
    "Pathway Enrichment",
    ha="center",
    fontsize=14,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.5",
        facecolor="#EBDEF0",
        edgecolor="black"
    )
)

ax.text(
    8,
    3.9,
    "Insulin Signaling | Type II Diabetes | PPAR Signaling",
    ha="center",
    fontsize=11
)

# arrow downward

ax.arrow(
    8,
    3.3,
    0,
    -0.7,
    length_includes_head=True,
    head_width=0.15,
    head_length=0.15,
    linewidth=2,
    color="black"
)

# ==================================================
# FINAL OUTPUT
# ==================================================

ax.text(
    8,
    2.0,
    "CAS-X V5 PRIORITIZATION",
    ha="center",
    fontsize=17,
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
    1.25,
    "Top Prioritized Genes",
    ha="center",
    fontsize=12,
    fontweight="bold"
)

ax.text(
    8,
    0.75,
    "FTO | JAZF1 | VEGFA | NOTCH2 | TSPAN8",
    ha="center",
    fontsize=12,
    color="darkred",
    fontweight="bold"
)

# ==================================================
# COVERAGE IMPROVEMENT
# ==================================================

ax.text(
    1.3,
    0.6,
    "+175%\nCoverage\nImprovement",
    ha="center",
    fontsize=10,
    fontweight="bold"
)

# ==================================================
# TITLE
# ==================================================

plt.title(
    "Graphical Abstract: CAS-X Framework for Prioritizing Type 2 Diabetes Genes",
    fontsize=18,
    pad=20
)

plt.tight_layout()

plt.savefig(
    "results/figures/Graphical_Abstract_CASX_V3.png",
    dpi=450,
    bbox_inches="tight"
)

plt.savefig(
    "results/figures/Graphical_Abstract_CASX_V3.pdf",
    bbox_inches="tight"
)

print("Saved Graphical Abstract V3")
