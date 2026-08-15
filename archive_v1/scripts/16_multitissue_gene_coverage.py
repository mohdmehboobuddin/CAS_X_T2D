import pandas as pd

genes = [
"ADAMTS9","ADCY5","CDKAL1","CDKN2A","CDKN2B",
"CENTD2","CHCHD9","CRY2","FAF1","GCK","GCKR",
"GRB14","HHEX","HNF1B","HNF4A","IGF2BP2",
"IRS1","KCNK16","KCNQ1","KLF14","LOC646736",
"MTNR1B","PPARG","SLC30A8","SUGP1","TCF7L2",
"TMEM154","ZFAND6","ZMIZ1"
]

files = [
("Adipose_Subcutaneous",
 "GTEx_Analysis_v11_eQTL/Adipose_Subcutaneous.v11.eQTLs.SuSiE_summary.parquet"),

("Adipose_Visceral",
 "GTEx_Analysis_v11_eQTL/Adipose_Visceral_Omentum.v11.eQTLs.SuSiE_summary.parquet"),

("Liver",
 "GTEx_Analysis_v11_eQTL/Liver.v11.eQTLs.SuSiE_summary.parquet"),

("Muscle",
 "GTEx_Analysis_v11_eQTL/Muscle_Skeletal.v11.eQTLs.SuSiE_summary.parquet"),

("Blood",
 "GTEx_Analysis_v11_eQTL/Whole_Blood.v11.eQTLs.SuSiE_summary.parquet")
]

for tissue, file in files:

    print("\n====================")
    print(tissue)
    print("====================")

    df = pd.read_parquet(file)

    gene_set = set(df["gene_name"])

    hits = [g for g in genes if g in gene_set]

    print("Hits:", len(hits))

    if hits:
        print(sorted(hits))
