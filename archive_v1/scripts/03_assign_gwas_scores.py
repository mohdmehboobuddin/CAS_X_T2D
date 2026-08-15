import pandas as pd

file_path = "data/processed/casx_master_dataset.csv"

df = pd.read_csv(file_path)

# All seed loci are replicated genome-wide significant T2D loci

df["GWAS_SCORE"] = 2

df.to_csv(file_path, index=False)

print("GWAS scores assigned.")
print(df[["SNP", "GENE", "GWAS_SCORE"]])
