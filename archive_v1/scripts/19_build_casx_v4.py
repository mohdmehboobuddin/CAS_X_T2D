import pandas as pd

print("\n=== CAS-X V4 ===\n")

df = pd.read_csv(
    "results/tables/casx_v3_rankings.csv"
)

df["GWAS_SCORE"] = 1

df["CASX_V4_SCORE"] = (
    df["GWAS_SCORE"]
    +
    df["TISSUE_COUNT"]
    +
    df["EXPRESSION_SCORE"]
)

df = df.sort_values(
    by="CASX_V4_SCORE",
    ascending=False
)

df["RANK_V4"] = range(
    1,
    len(df) + 1
)

print(
    df[
        [
            "GENE",
            "GWAS_SCORE",
            "TISSUE_COUNT",
            "EXPRESSION_SCORE",
            "CASX_V4_SCORE",
            "RANK_V4"
        ]
    ].head(25)
)

df.to_csv(
    "results/tables/casx_v4_rankings.csv",
    index=False
)

print("\nSaved:")
print("results/tables/casx_v4_rankings.csv")
