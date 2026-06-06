import pandas as pd

# 50 CAS-X genes
genes = pd.read_csv(
    "data/processed/master_t2d_loci_50.csv"
)

casx_genes = set(genes["GENE"])

# GTEx-supported genes
gtex = pd.read_csv(
    "results/tables/eqtl_scores.csv"
)

gtex_genes = set(gtex["gene_name"])

present = casx_genes.intersection(gtex_genes)
missing = casx_genes - gtex_genes

print("\n=== CAS-X Gene Coverage ===\n")

print("Total CAS-X genes:", len(casx_genes))
print("GTEx-supported genes:", len(present))
print("Missing genes:", len(missing))

print("\nPresent:")
print(sorted(present))

print("\nMissing:")
print(sorted(missing))
