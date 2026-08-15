import pandas as pd
import matplotlib.pyplot as plt

v4 = pd.read_csv(
    "results/tables/casx_v4_rankings.csv"
)[["GENE", "RANK_V4"]]

v5 = pd.read_csv(
    "results/tables/casx_v5_rankings.csv"
)[["GENE", "RANK_V5"]]

df = v4.merge(v5, on="GENE")

plt.figure(figsize=(8,7))

plt.scatter(
    df["RANK_V4"],
    df["RANK_V5"],
    s=80
)

for _, row in df.iterrows():
    if row["RANK_V5"] <= 10:
        plt.annotate(
            row["GENE"],
            (row["RANK_V4"], row["RANK_V5"]),
            fontsize=8
        )

plt.plot(
    [1, 37],
    [1, 37],
    linestyle="--"
)

plt.xlabel("CAS-X V4 Rank")
plt.ylabel("CAS-X V5 Rank")

plt.title(
    "Figure 8. Robustness of Gene Rankings Between CAS-X V4 and V5"
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure8_Robustness_V4_V5.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 8")
