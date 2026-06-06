import pandas as pd

data = [
    ["FTO",1],
    ["IRS1",1],
    ["KCNJ11",1],
    ["PPARG",1],
    ["TCF7L2",1],
    ["VEGFA",0]
]

df = pd.DataFrame(
    data,
    columns=[
        "GENE",
        "CRISPR_SUPPORT"
    ]
)

print("\n=== CRISPR Validation ===\n")
print(df)

df.to_csv(
    "results/tables/crispr_validation.csv",
    index=False
)

print("\nSaved:")
print(
    "results/tables/crispr_validation.csv"
)
