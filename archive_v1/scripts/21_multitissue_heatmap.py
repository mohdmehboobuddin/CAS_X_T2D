import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/casx_v4_rankings.csv"
)

top = df.head(20)

heat = top[
    [
        "PANCREAS",
        "ADIPOSE_SUBQ",
        "ADIPOSE_VISC",
        "MUSCLE",
        "BLOOD"
    ]
]

plt.figure(figsize=(9,8))

im = plt.imshow(
    heat,
    aspect="auto"
)

plt.colorbar(im, label="Support")

plt.yticks(
    range(len(top)),
    top["GENE"]
)

plt.xticks(
    range(5),
    [
        "Pancreas",
        "Adipose SubQ",
        "Adipose Visceral",
        "Muscle",
        "Blood"
    ],
    rotation=45
)

plt.title(
    "CAS-X Multi-Tissue GTEx Support Matrix"
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure3_Multitissue_Heatmap.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure3")
