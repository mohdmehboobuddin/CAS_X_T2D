import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# Load only the columns we actually need
v4 = pd.read_csv(
    "results/tables/casx_v4_rankings.csv"
)[["GENE", "RANK_V4"]]

v5 = pd.read_csv(
    "results/tables/casx_v5_rankings.csv"
)[["GENE", "RANK_V5"]]

# Merge
df = v4.merge(v5, on="GENE")

# Calculate correlation
rho, pval = spearmanr(
    df["RANK_V4"],
    df["RANK_V5"]
)

# Plot
plt.figure(figsize=(8, 8))

plt.scatter(
    df["RANK_V4"],
    df["RANK_V5"],
    s=80
)

# Label top genes
top_genes = [
    "FTO",
    "JAZF1",
    "VEGFA",
    "NOTCH2",
    "TSPAN8"
]

for gene in top_genes:

    row = df[df["GENE"] == gene]

    plt.annotate(
        gene,
        (
            row["RANK_V4"].values[0],
            row["RANK_V5"].values[0]
        ),
        fontsize=9
    )

# Identity line
plt.plot(
    [1, 37],
    [1, 37],
    linestyle="--",
    linewidth=1.5
)

# Correlation text
plt.text(
    18,
    5,
    f"Spearman r = {rho:.2f}",
    fontsize=11,
    bbox=dict(facecolor="white", alpha=0.8)
)

plt.xlabel("CAS-X V4 Rank", fontsize=12)
plt.ylabel("CAS-X V5 Rank", fontsize=12)

plt.title(
    "Figure 8. Robustness of CAS-X Rankings",
    fontsize=15
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure8_Robustness_V2.png",
    dpi=450,
    bbox_inches="tight"
)

print("\nSaved:")
print("results/figures/Figure8_Robustness_V2.png")
print(f"Spearman correlation = {rho:.3f}")
