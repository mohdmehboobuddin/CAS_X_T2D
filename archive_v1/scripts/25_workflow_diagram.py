import matplotlib.pyplot as plt

plt.figure(figsize=(8,8))

steps = [
    "50 T2D GWAS Loci",
    "37 Unique Genes",
    "GTEx Multi-Tissue eQTL Analysis",
    "Human Islet Expression Validation",
    "CAS-X V4 Scoring",
    "Prioritized Candidate Genes"
]

for i, step in enumerate(steps):

    y = len(steps) - i

    plt.text(
        0.5,
        y,
        step,
        ha="center",
        bbox=dict(boxstyle="round")
    )

    if i < len(steps)-1:

        plt.arrow(
            0.5,
            y-0.2,
            0,
            -0.5,
            length_includes_head=True
        )

plt.axis("off")

plt.title(
    "CAS-X Framework Workflow"
)

plt.savefig(
    "results/figures/Figure1_Workflow.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure1")
