import pandas as pd
import matplotlib.pyplot as plt

# Validation counts

validation = pd.DataFrame({
    "Validation": [
        "GTEx\n(GSE38642)",
        "scRNA-seq\n(GSE81608)",
        "Open\nTargets",
        "CRISPR"
    ],
    "Supported_Genes": [
        8,   # GSE38642
        5,   # GSE81608
        8,   # Open Targets
        5    # CRISPR
    ]
})

plt.figure(figsize=(8,6))

bars = plt.bar(
    validation["Validation"],
    validation["Supported_Genes"]
)

for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width()/2,
        height + 0.1,
        str(int(height)),
        ha="center",
        fontsize=11
    )

plt.ylabel("Number of Supported Genes")
plt.title(
    "Figure 6. Independent Validation of CAS-X Prioritized Genes"
)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure6_Validation_Summary.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure 6")
