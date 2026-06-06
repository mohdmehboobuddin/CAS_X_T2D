import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/casx_multitissue_support.csv"
)

tissues = [
    "PANCREAS",
    "ADIPOSE_SUBQ",
    "ADIPOSE_VISC",
    "MUSCLE",
    "BLOOD"
]

counts = [
    df[t].sum()
    for t in tissues
]

total_genes = len(df)

percentages = [
    round(x/total_genes*100,1)
    for x in counts
]

plt.figure(figsize=(9,6))

bars = plt.bar(
    [
        "Pancreas",
        "Adipose\nSubQ",
        "Adipose\nVisceral",
        "Muscle",
        "Blood"
    ],
    counts
)

for bar, c, p in zip(
    bars,
    counts,
    percentages
):

    plt.text(
        bar.get_x()+bar.get_width()/2,
        c+0.3,
        f"{c}\n({p}%)",
        ha="center",
        fontsize=9
    )

plt.ylabel("Supported Genes")

plt.title(
    "Figure 5. Tissue Contribution to CAS-X Gene Support",
    fontsize=15
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure5_Tissue_Contribution_V2.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 5")
