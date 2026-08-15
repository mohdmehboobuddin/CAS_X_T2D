import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/pathway_enrichment_summary.csv"
)

counts = (
    df["PATHWAY"]
    .value_counts()
    .head(10)
)

plt.figure(figsize=(12,8))

bars = plt.barh(
    counts.index[::-1],
    counts.values[::-1]
)

for bar in bars:
    width = bar.get_width()

    plt.text(
        width + 0.05,
        bar.get_y()+bar.get_height()/2,
        str(int(width)),
        va="center"
    )

plt.xlabel("Number of Supporting Databases")

plt.title(
    "Figure 7. Biological Pathways Enriched Among CAS-X Prioritized Genes",
    fontsize=15
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure7_Pathway_Enrichment_V2.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 7")
