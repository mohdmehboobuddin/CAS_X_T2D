import pandas as pd

df = pd.read_csv(
    "results/tables/eqtl_candidates.csv"
)

def assign_score(pip):

    if pip > 0.5:
        return 3

    elif pip > 0.3:
        return 2

    elif pip > 0.1:
        return 1

    else:
        return 0

df["EQTL_SCORE"] = df["pip"].apply(assign_score)

df.to_csv(
    "results/tables/eqtl_scores.csv",
    index=False
)

print(df)
