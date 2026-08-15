import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/casx_v4_rankings.csv"
)

top10 = df.head(10)

plt.figure(figsize=(8,6))

bars = plt.barh(
    top10["GENE"][::-1],
    top10["CASX_V4_SCORE"][::-1]
)

for bar in bars:

    width = bar.get_width()

    plt.text(
        width + 0.1,
        bar.get_y() + bar.get_height()/2,
        f"{width:.0f}",
        va="center"
    )

plt.xlabel("CAS-X V4 Score")

plt.title(
    "Top 10 CAS-X Prioritized Genes"
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure4_Top10_CASX_V4.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure4")
