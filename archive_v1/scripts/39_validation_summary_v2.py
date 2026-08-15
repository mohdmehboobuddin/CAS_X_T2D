import matplotlib.pyplot as plt

sources = [
    "GTEx\n(GSE38642)",
    "scRNA-seq\n(GSE81608)",
    "Open\nTargets",
    "CRISPR"
]

supported = [8, 5, 8, 5]
total = 37

percentages = [
    round(x / total * 100, 1)
    for x in supported
]

plt.figure(figsize=(10,7))

bars = plt.bar(
    sources,
    percentages
)

for bar, count, pct in zip(
    bars,
    supported,
    percentages
):
    plt.text(
        bar.get_x()+bar.get_width()/2,
        bar.get_height()+1,
        f"{count}/37\n({pct}%)",
        ha="center",
        fontsize=10
    )

plt.ylabel("Validated Genes (%)")

plt.title(
    "Figure 6. Independent Validation of CAS-X Prioritized Genes",
    fontsize=15
)

plt.ylim(0,30)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure6_Validation_Summary_V2.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 6")
