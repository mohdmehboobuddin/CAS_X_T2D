import matplotlib.pyplot as plt

labels = [
    "Pancreas\nOnly",
    "Multi-Tissue\nGTEx"
]

values = [8, 22]

plt.figure(figsize=(8,6))

bars = plt.bar(
    labels,
    values
)

for bar, value in zip(bars, values):

    plt.text(
        bar.get_x() + bar.get_width()/2,
        value + 0.3,
        str(value),
        ha="center",
        fontsize=11,
        fontweight="bold"
    )

# Move annotation away from title
plt.text(
    0.5,
    23.2,
    "+175% Coverage Increase",
    ha="center",
    fontsize=13,
    fontweight="bold"
)

plt.ylabel(
    "Supported Genes",
    fontsize=12
)

plt.title(
    "Figure 2. Improvement in Gene Coverage After Multi-Tissue Integration",
    fontsize=14,
    pad=20
)

plt.ylim(0,26)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure2_Coverage_Improvement_V2.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 2")
