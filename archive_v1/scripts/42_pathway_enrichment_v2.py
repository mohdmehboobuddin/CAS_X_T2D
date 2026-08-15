import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Load enrichment results
# -----------------------------
df = pd.read_csv(
    "results/tables/pathway_enrichment_stats.csv"
)

# -----------------------------
# Compute -log10(adj p)
# -----------------------------
df["NEG_LOG10_P"] = -np.log10(df["ADJ_P"])

# Sort by significance
df = df.sort_values(
    "NEG_LOG10_P",
    ascending=True
)

# -----------------------------
# Database colors
# -----------------------------
color_map = {
    "KEGG": "#4C72B0",
    "WikiPathways": "#55A868",
    "Reactome": "#C44E52"
}

colors = [
    color_map[x]
    for x in df["DATABASE"]
]

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(11, 6))

bars = plt.barh(
    df["PATHWAY"],
    df["NEG_LOG10_P"],
    color=colors,
    edgecolor="black"
)

# annotate p-values
for bar, p in zip(bars, df["ADJ_P"]):

    plt.text(
        bar.get_width() + 0.03,
        bar.get_y() + bar.get_height()/2,
        f"adjP={p:.3f}",
        va="center",
        fontsize=8
    )

plt.xlabel(
    "-log10(Adjusted P-value)",
    fontsize=12,
    fontweight="bold"
)

plt.ylabel(
    "Enriched Pathway",
    fontsize=12,
    fontweight="bold"
)

plt.title(
    "Figure 7. Pathway Enrichment of CAS-X Prioritized Genes",
    fontsize=16,
    fontweight="bold",
    pad=18
)

# legend
from matplotlib.patches import Patch

legend_handles = [
    Patch(
        facecolor=color_map[k],
        edgecolor="black",
        label=k
    )
    for k in color_map
]

plt.legend(
    handles=legend_handles,
    title="Database",
    loc="lower right"
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure7_Pathway_Enrichment.png",
    dpi=450,
    bbox_inches="tight"
)

plt.close()

print("Saved Figure 7")
