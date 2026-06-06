import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 12))

ax.set_xlim(0, 10)
ax.set_ylim(0, 14)

ax.axis("off")

steps = [
    ("50 T2D GWAS Loci", 13),
    ("Gene Mapping\n(37 Candidate Genes)", 11),
    ("GTEx Pancreas eQTL\nIntegration", 9),
    ("Multi-Tissue GTEx\nValidation", 7),
    ("Expression Validation\n(GSE38642 + GSE81608)", 5),
    ("External Validation\n(Open Targets + CRISPR)", 3),
    ("CAS-X V5 Ranking\n& Prioritization", 1)
]

for text, y in steps:

    ax.text(
        5,
        y,
        text,
        ha="center",
        va="center",
        fontsize=13,
        bbox=dict(
            boxstyle="round,pad=0.6",
            facecolor="#4C72B0",
            edgecolor="black"
        )
    )

for y in [12,10,8,6,4,2]:

    ax.arrow(
        5,
        y + 0.35,
        0,
        -0.7,
        length_includes_head=True,
        head_width=0.25,
        head_length=0.25,
        linewidth=2,
        color="black"
    )

plt.title(
    "Figure 1. CAS-X Framework Workflow",
    fontsize=18,
    pad=25
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure1_CASX_Workflow_V2.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 1")
