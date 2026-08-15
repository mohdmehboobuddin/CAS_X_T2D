import pandas as pd

print("\n=== CAS-X V2 Integration ===\n")

# Load eQTL evidence
eqtl = pd.read_csv(
    "results/tables/eqtl_scores.csv"
)

# Load expression evidence
expr = pd.read_csv(
    "results/tables/gse38642_expression_support.csv"
)

# Assign expression scores
def expression_score(x):
    if x >= 9:
        return 3
    elif x >= 8:
        return 2
    elif x >= 7:
        return 1
    else:
        return 0

expr["EXPRESSION_SCORE"] = expr["MEAN_EXPRESSION"].apply(expression_score)

# Merge datasets
merged = pd.merge(
    eqtl,
    expr[["GENE", "MEAN_EXPRESSION", "EXPRESSION_SCORE"]],
    left_on="gene_name",
    right_on="GENE",
    how="inner"
)

# CAS-X V2 score
merged["CASX_V2_SCORE"] = (
    merged["EQTL_SCORE"] +
    merged["EXPRESSION_SCORE"]
)

# Sort and rank
merged = merged.sort_values(
    by="CASX_V2_SCORE",
    ascending=False
)

merged["RANK"] = range(1, len(merged) + 1)

# Display results
print(
    merged[
        [
            "gene_name",
            "EQTL_SCORE",
            "EXPRESSION_SCORE",
            "CASX_V2_SCORE",
            "RANK"
        ]
    ]
)

# Save results
merged.to_csv(
    "results/tables/casx_v2_rankings.csv",
    index=False
)

print("\nSaved:")
print("results/tables/casx_v2_rankings.csv")
