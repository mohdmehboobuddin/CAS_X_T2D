import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load pathway enrichment results
df = pd.read_csv(
    "results/tables/pathway_enrichment_stats.csv"
)

# Calculate -log10(adjP)
df["NEGLOG10_P"] = -np.log10(df["ADJ_P"])

# Sort for visualization
df = df.sort_values(
    "NEGLOG10_P",
    ascending=True
)

# Database colors
colors = []

for db in df["DATABASE"]:

    if db == "KEGG":
        colors.append("#4C72B0")

    elif db == "WikiPathways":
        colors.append("#55A868")

    elif db == "Reactome":
        colors.append("#C44E52")

    else:
        colors.append("gray")

# Create figure
plt.figure(figsize=(12, 7))

bars = plt.barh(
    df["PATHWAY"],
    df["NEGLOG10_P"],
    color=colors,
    edgecolor="black"
)

# Add adjusted P-value labels
for i, row in enumerate(df.itertuples()):

    plt.text(
        row.NEGLOG10_P + 0.03,
        i,
        f"adjP={row.ADJ_P:.3f}",
        va="center",
        fontsize=10
    )

# Axis labels
plt.xlabel(
    "-log10(Adjusted P-value)",
    fontsize=14,
    fontweight="bold"
)

plt.ylabel(
    "Enriched Pathway",
    fontsize=14,
    fontweight="bold"
)

# IMPORTANT:
# No figure number in title
plt.title(
    "Pathway Enrichment of CAS-X Prioritized Genes",
    fontsize=18,
    fontweight="bold"
)

# Legend
from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor="#4C72B0", edgecolor="black", label="KEGG"),
    Patch(facecolor="#55A868", edgecolor="black", label="WikiPathways"),
    Patch(facecolor="#C44E52", edgecolor="black", label="Reactome")
]

plt.legend(
    handles=legend_elements,
    title="Database",
    loc="lower right",
    fontsize=11,
    title_fontsize=12
)

plt.tight_layout()

# Save as Figure 8
plt.savefig(
    "results/figures/Figure8_Pathway_Enrichment.png",
    dpi=450,
    bbox_inches="tight"
)

plt.close()

print("Saved Figure 8")
