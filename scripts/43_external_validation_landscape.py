import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "DejaVu Sans"

sources = [
    "GTEx\n(GSE38642)",
    "scRNA-seq\n(GSE81608)",
    "Open\nTargets",
    "CRISPR"
]

counts = [8, 5, 8, 5]

fig, ax = plt.subplots(figsize=(12,7))

bars = ax.bar(
    sources,
    counts,
    edgecolor="black",
    linewidth=1.5
)

for bar, value in zip(bars, counts):
    ax.text(
        bar.get_x() + bar.get_width()/2,
        value + 0.15,
        str(value),
        ha="center",
        fontsize=14,
        fontweight="bold"
    )

ax.set_ylim(0,10)

ax.set_ylabel(
    "Supported Genes",
    fontsize=16,
    fontweight="bold"
)

ax.set_title(
    "Figure 7. External Validation of CAS-X Prioritized Genes",
    fontsize=20,
    fontweight="bold",
    pad=20
)

plt.figtext(
    0.5,
    0.02,
    "Genes Supported by Multiple Independent Resources:\n"
    "KCNJ11   |   PPARG   |   TCF7L2",
    ha="center",
    fontsize=15,
    fontweight="bold"
)

plt.tight_layout(rect=[0,0.08,1,1])

plt.savefig(
    "results/figures/Figure7_External_Validation_Landscape.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 7")
