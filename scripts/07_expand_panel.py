import pandas as pd

df = pd.read_csv("data/processed/t2d_50_loci_candidate_panel.csv")

print("\n=== CAS-X Expanded Panel ===\n")

print("Total SNPs:", len(df))
print("Unique genes:", df["GENE"].nunique())

print("\nTop genes represented:\n")
print(df["GENE"].value_counts().head(10))
