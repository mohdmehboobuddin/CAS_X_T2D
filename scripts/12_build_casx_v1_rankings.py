import pandas as pd

# Load eQTL scores
eqtl = pd.read_csv(
    "results/tables/eqtl_scores.csv"
)

# Simple CAS-X v1 score
eqtl["CASX_V1_SCORE"] = eqtl["EQTL_SCORE"] + 2

# Rank
eqtl = eqtl.sort_values(
    "CASX_V1_SCORE",
    ascending=False
)

eqtl["RANK"] = range(1, len(eqtl)+1)

# Save
eqtl.to_csv(
    "results/tables/casx_v1_rankings.csv",
    index=False
)

print(eqtl)
