import pandas as pd

# CAS-X top genes and Entrez IDs

gene_map = {
    "FTO":79068,
    "IRS1":3667,
    "PPARG":5468,
    "KCNJ11":3767,
    "VEGFA":7422
}

df = pd.read_csv(
    "data/raw/expression/GSE81608_human_islets_rpkm.txt.gz",
    sep="\t"
)

results = []

for gene, gid in gene_map.items():

    row = df[
        df["gene.id"] == gid
    ]

    if len(row) == 0:
        continue

    expr = row.iloc[:,1:].mean(axis=1).values[0]

    results.append(
        [gene, gid, expr]
    )

out = pd.DataFrame(
    results,
    columns=[
        "GENE",
        "ENTREZ_ID",
        "MEAN_EXPRESSION"
    ]
)

out = out.sort_values(
    "MEAN_EXPRESSION",
    ascending=False
)

print(
    "\n=== GSE81608 Validation ===\n"
)

print(out)

out.to_csv(
    "results/tables/gse81608_validation.csv",
    index=False
)

print("\nSaved:")
print(
    "results/tables/gse81608_validation.csv"
)
