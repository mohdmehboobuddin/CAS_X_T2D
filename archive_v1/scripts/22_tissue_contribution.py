import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/casx_multitissue_support.csv"
)

counts = {
    "Blood": df["BLOOD"].sum(),
    "Muscle": df["MUSCLE"].sum(),
    "Adipose SubQ": df["ADIPOSE_SUBQ"].sum(),
    "Pancreas": df["PANCREAS"].sum(),
    "Adipose Visceral": df["ADIPOSE_VISC"].sum()
}

names = list(counts.keys())
values = list(counts.values())

plt.figure(figsize=(8,5))

bars = plt.bar(names, values)

for bar, val in zip(bars, values):

    plt.text(
        bar.get_x() + bar.get_width()/2,
        val + 0.2,
        str(val),
        ha="center"
    )

plt.ylabel("Supported Genes")

plt.title(
    "Gene Recovery Across GTEx Tissues"
)

plt.xticks(rotation=25)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure5_Tissue_Contribution.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure5")
