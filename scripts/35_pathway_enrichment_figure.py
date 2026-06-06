import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/pathway_enrichment_summary.csv"
)

counts = (
    df["CATEGORY"]
    .value_counts()
    .reset_index()
)

counts.columns = ["CATEGORY", "COUNT"]

plt.figure(figsize=(9,6))

bars = plt.bar(
    counts["CATEGORY"],
    counts["COUNT"]
)

for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width()/2,
        height + 0.05,
        str(int(height)),
        ha="center"
    )

plt.ylabel("Number of Enriched Pathways")
plt.xlabel("Pathway Database")

plt.title(
    "Figure 7. Pathway Enrichment Across Multiple Databases"
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure7_Pathway_Enrichment.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 7")
