import pandas as pd

missing = [
    "TCF7L2",
    "SLC30A8",
    "MTNR1B",
    "KCNQ1",
    "IGF2BP2",
    "HNF1B",
    "HNF4A"
]

files = {
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

for tissue, file in files.items():

    print("\n", "="*20)
    print(tissue)
    print("="*20)

    df = pd.read_parquet(file)

    genes = set(df["gene_name"])

    hits = []

    for g in missing:

        if g in genes:
            hits.append(g)

    print("Found:", hits)
