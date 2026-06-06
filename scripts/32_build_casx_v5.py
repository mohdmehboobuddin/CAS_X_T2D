import pandas as pd

df = pd.read_csv(
    "results/tables/casx_v4_rankings.csv"
)

# New weighting:
# GWAS = 1
# Tissue count = 1x
# Expression = 2x

df["CASX_V5_SCORE"] = (
    1
    + df["TISSUE_COUNT"]
    + (2 * df["EXPRESSION_SCORE"])
)

df = df.sort_values(
    "CASX_V5_SCORE",
    ascending=False
)

df["RANK_V5"] = range(
    1,
    len(df) + 1
)

print("\n=== CAS-X V5 ===\n")

print(
    df[
        [
            "GENE",
            "CASX_V4_SCORE",
            "CASX_V5_SCORE",
            "RANK_V5"
        ]
    ].head(20)
)

df.to_csv(
    "results/tables/casx_v5_rankings.csv",
    index=False
)

print("\nSaved:")
print(
    "results/tables/casx_v5_rankings.csv"
)
