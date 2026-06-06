import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(18, 6))

ax.set_xlim(0, 18)
ax.set_ylim(0, 6)

ax.axis("off")

steps = [
    ("GWAS\n50 T2D Loci", 1),
    ("Gene Mapping\n37 Genes", 4),
    ("GTEx eQTL\nIntegration", 7),
    ("Multi-Tissue\nValidation", 10),
    ("Expression\nValidation", 13),
    ("External\nValidation", 16)
]

for text, x in steps:

    ax.text(
        x,
        3,
        text,
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.6",
            facecolor="#D9EAF7",
            edgecolor="black",
            linewidth=1.5
        )
    )

# arrows

for x in [2.3, 5.3, 8.3, 11.3, 14.3]:

    ax.arrow(
        x,
        3,
        1,
        0,
        length_includes_head=True,
        head_width=0.15,
        head_length=0.2,
        linewidth=2,
        color="black"
    )

# Validation notes

ax.text(
    16,
    1.3,
    "Open Targets\n8 Genes",
    ha="center",
    fontsize=10
)

ax.text(
    16,
    0.5,
    "CRISPR Support\n5 Genes",
    ha="center",
    fontsize=10
)

# Final output box

ax.text(
    9,
    5,
    "CAS-X V5 PRIORITIZATION",
    fontsize=16,
    fontweight="bold",
    ha="center",
    bbox=dict(
        boxstyle="round,pad=0.6",
        facecolor="#FFF3C4",
        edgecolor="black",
        linewidth=2
    )
)

# Top genes

ax.text(
    9,
    4.2,
    "FTO  |  JAZF1  |  VEGFA  |  NOTCH2  |  TSPAN8",
    fontsize=12,
    ha="center",
    fontweight="bold"
)

# Coverage improvement

ax.text(
    9,
    0.2,
    "+175% Coverage Improvement via Multi-Tissue Integration",
    fontsize=11,
    ha="center",
    fontweight="bold"
)

plt.title(
    "Graphical Abstract: CAS-X Framework for Prioritization of Type 2 Diabetes Genes",
    fontsize=18,
    pad=20
)

plt.tight_layout()

plt.savefig(
    "results/figures/Graphical_Abstract_CASX.png",
    dpi=450,
    bbox_inches="tight"
)

plt.savefig(
    "results/figures/Graphical_Abstract_CASX.pdf",
    bbox_inches="tight"
)

print("Saved:")
print("results/figures/Graphical_Abstract_CASX.png")
print("results/figures/Graphical_Abstract_CASX.pdf")
