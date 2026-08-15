import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/casx_v4_rankings.csv"
)

plt.figure(figsize=(7,5))

plt.hist(
    df["CASX_V4_SCORE"],
    bins=10
)

plt.xlabel("CAS-X V4 Score")
plt.ylabel("Number of Genes")

plt.title(
    "Distribution of CAS-X V4 Scores"
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure4_CASX_Distribution.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure4")
