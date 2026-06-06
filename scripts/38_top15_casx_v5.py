import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/casx_v5_rankings.csv"
)

top = df.head(15)

plt.figure(figsize=(10,8))

bars = plt.barh(
    top["GENE"][::-1],
    top["CASX_V5_SCORE"][::-1]
)

for bar in bars:

    width = bar.get_width()

    plt.text(
        width + 0.1,
        bar.get_y() + bar.get_height()/2,
        f"{width:.0f}",
        va="center",
        fontsize=10
    )

plt.xlabel("CAS-X V5 Score")
plt.ylabel("Gene")

plt.title(
    "Figure 4. Top 15 CAS-X V5 Prioritized Genes",
    fontsize=16
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure4_Top15_CASX_V5.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 4")
