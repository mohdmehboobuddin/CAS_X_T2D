import pandas as pd

print("\n=== CAS-X V3 Integration ===\n")

# Load multi-tissue matrix
multi = pd.read_csv(
    "results/tables/casx_multitissue_support.csv"
)

# Load expression support
expr = pd.read_csv(
    "results/tables/gse38642_expression_support.csv"
)

# Expression scoring
def expression_score(x):

    if x >= 9:
        return 3

    elif x >= 8:
        return 2

    elif x >= 7:
        return 1

    else:
        return 0

expr["EXPRESSION_SCORE"] = (
    expr["MEAN_EXPRESSION"]
    .apply(expression_score)
)

# Keep only required columns
expr = expr[
    ["GENE", "MEAN_EXPRESSION", "EXPRESSION_SCORE"]
]

# Merge
merged = pd.merge(
    multi,
    expr,
    on="GENE",
    how="left"
)

# Missing expression evidence = 0
merged["EXPRESSION_SCORE"] = (
    merged["EXPRESSION_SCORE"]
    .fillna(0)
)

# CAS-X V3 score
merged["CASX_V3_SCORE"] = (
    merged["TISSUE_COUNT"]
    +
    merged["EXPRESSION_SCORE"]
)

# Sort
merged = merged.sort_values(
    by="CASX_V3_SCORE",
    ascending=False
)

merged["RANK"] = range(
    1,
    len(merged) + 1
)

print(
    merged[
        [
            "GENE",
            "TISSUE_COUNT",
            "EXPRESSION_SCORE",
            "CASX_V3_SCORE",
            "RANK"
        ]
    ].head(20)
)

merged.to_csv(
    "results/tables/casx_v3_rankings.csv",
    index=False
)

print("\nSaved:")
print("results/tables/casx_v3_rankings.csv")
