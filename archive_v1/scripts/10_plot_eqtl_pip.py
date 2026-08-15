import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/tables/eqtl_candidates.csv"
)

plt.figure(figsize=(8,5))

plt.bar(
    df["gene_name"],
    df["pip"]
)

plt.ylabel("Maximum PIP")
plt.xlabel("Gene")
plt.title("GTEx Pancreas Fine-Mapped eQTL Support")

plt.xticks(rotation=45)

plt.tight_layout()

plt.savefig(
    "results/figures/Figure1_GTEx_PIP.png",
    dpi=300
)

print("Figure saved:")
print("results/figures/Figure1_GTEx_PIP.png")
