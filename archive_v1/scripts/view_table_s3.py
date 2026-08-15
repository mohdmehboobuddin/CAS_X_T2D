import pandas as pd
from pathlib import Path

# Set pandas options to display all 50 rows and full width
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100) 
pd.set_option('display.width', 1000)

table_path = Path("results/Supplementary_Tables/TableS3_Curated_50_GWAS_Loci.csv")

if table_path.exists():
    print(f"\n{'='*60}")
    print(f"📄 TABLE: {table_path.name}")
    print(f"{'='*60}")
    df = pd.read_csv(table_path)
    # Print the entire table without the index numbers
    print(df.to_string(index=False))
    print(f"{'='*60}")
    print(f"Total loci recovered: {len(df)}\n")
else:
    print(f"❌ Could not find: {table_path}")
