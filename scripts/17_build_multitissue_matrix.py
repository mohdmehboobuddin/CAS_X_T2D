import pandas as pd

# CAS-X genes
genes = pd.read_csv(
    "data/processed/master_t2d_loci_50.csv"
)

genes = sorted(set(genes["GENE"]))

matrix = pd.DataFrame({
    "GENE": genes
})

tissues = {
    "PANCREAS":
    "data/raw/eqtl/extracted/Pancreas.v11.eQTLs.SuSiE_summary.parquet",

    "ADIPOSE_SUBQ":
    "GTEx_Analysis_v11_eQTL/Adipose_Subcutaneous.v11.eQTLs.SuSiE_summary.parquet",

    "ADIPOSE_VISC":
    "GTEx_Analysis_v11_eQTL/Adipose_Visceral_Omentum.v11.eQTLs.SuSiE_summary.parquet",

    "MUSCLE":
    "GTEx_Analysis_v11_eQTL/Muscle_Skeletal.v11.eQTLs.SuSiE_summary.parquet",

    "BLOOD":
    "GTEx_Analysis_v11_eQTL/Whole_Blood.v11.eQTLs.SuSiE_summary.parquet"
}

for tissue, file in tissues.items():

    print(f"Loading {tissue}")

    df = pd.read_parquet(file)

    gene_set = set(df["gene_name"])

    matrix[tissue] = matrix["GENE"].apply(
        lambda x: 1 if x in gene_set else 0
    )

matrix["TISSUE_COUNT"] = matrix[
    ["PANCREAS",
     "ADIPOSE_SUBQ",
     "ADIPOSE_VISC",
     "MUSCLE",
     "BLOOD"]
].sum(axis=1)

matrix = matrix.sort_values(
    "TISSUE_COUNT",
    ascending=False
)

print("\n=== CAS-X Multi-Tissue Matrix ===\n")
print(matrix.head(20))

matrix.to_csv(
    "results/tables/casx_multitissue_support.csv",
    index=False
)

print("\nSaved:")
print("results/tables/casx_multitissue_support.csv")
