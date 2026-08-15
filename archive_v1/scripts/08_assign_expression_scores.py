import pandas as pd

loci = pd.read_csv("data/processed/master_t2d_loci_50.csv")

expression_map = {
    "SLC30A8": 2,
    "TCF7L2": 2,
    "KCNJ11": 1,
    "CDKAL1": 1,
    "IGF2BP2": 1,
    "PPARG": 1,
    "MTNR1B": 0,
    "IRS1": 0,
    "KLF14": 0
}

loci["EXPRESSION_SCORE"] = loci["GENE"].map(expression_map).fillna(1)

print("\nExpression Layer Summary\n")
print(loci["EXPRESSION_SCORE"].value_counts())

loci.to_csv(
    "data/processed/expression_layer.csv",
    index=False
)

print("\nSaved: data/processed/expression_layer.csv")
