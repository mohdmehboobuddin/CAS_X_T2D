import pandas as pd
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

OLD_TABLE2 = os.path.join(PROJECT_ROOT, "results", "manuscript_tables", "Table2_Top15_CASX_Genes.csv")
OLD_TABLE4 = os.path.join(PROJECT_ROOT, "results", "manuscript_tables", "Table4_Pathway_Enrichment.csv")

print("======================================================================")
print("  OLD TABLE 2: TOP 15 GENES (PREVIOUS MANUAL/V5 SYSTEM)")
print("======================================================================")
if os.path.exists(OLD_TABLE2):
    df_t2 = pd.read_csv(OLD_TABLE2)
    print(df_t2.head(15).to_string(index=False))
else:
    print("Old Table 2 file not found.")

print("\n======================================================================")
print("  OLD TABLE 4: PATHWAY ENRICHMENT (PREVIOUS MANUAL/V5 SYSTEM)")
print("======================================================================")
if os.path.exists(OLD_TABLE4):
    df_t4 = pd.read_csv(OLD_TABLE4)
    print(df_t4.head(10).to_string(index=False))
else:
    print("Old Table 4 file not found.")
print("======================================================================\n")
