import pandas as pd

input_file = "data/processed/master_t2d_loci.csv"

df = pd.read_csv(input_file)

print("\n=== CAS-X Seed Dataset ===\n")
print(df.head())

print(f"\nTotal loci: {len(df)}")

