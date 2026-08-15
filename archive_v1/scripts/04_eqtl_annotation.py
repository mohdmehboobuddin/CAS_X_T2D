import pandas as pd

df = pd.read_csv("data/processed/casx_master_dataset.csv")

print("\n=== CAS-X eQTL Module ===")
print(f"Loci loaded: {len(df)}")

print("\nNext step: integrate GTEx / islet eQTL evidence.")
