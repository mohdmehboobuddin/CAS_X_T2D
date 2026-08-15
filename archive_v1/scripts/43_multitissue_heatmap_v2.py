import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/casx_multitissue_support.csv"
)

# Top 20 genes
df = df.head(20)

heatmap_data = df[
    [
        "PANCREAS",
        "ADIPOSE_SUBQ",
        "ADIPOSE_VISC",
        "MUSCLE",
        "BLOOD"
    ]
]

plt.figure(figsize=(9,12))

im = plt.imshow(
    heatmap_data,
    aspect="auto",
    cmap="Blues"
)

plt.colorbar(
    im,
    label="Support"
)

plt.yticks(
    range(len(df)),
    df["GENE"],
    fontsize=12
)

plt.xticks(
    range(5),
    [
        "Pancreas",
        "Adipose\nSubQ",
        "Adipose\nVisceral",
        "Muscle",
        "Blood"
    ],
    fontsize=11
)

plt.title(
    "Figure 3. Multi-Tissue Support of Top CAS-X Prioritized Genes",
    fontsize=15,
    pad=20
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure3_Multitissue_Heatmap_V2.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 3")
