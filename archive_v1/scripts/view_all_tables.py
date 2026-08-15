import pandas as pd
from pathlib import Path

# Set pandas display options for better terminal viewing
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

base_dir = Path("results")

tables_to_view = [
    base_dir / "Main_Manuscript_Tables" / "Table1_Target_Prioritization_Summary.csv",
    base_dir / "Main_Manuscript_Tables" / "Table2_Clinical_Tractability_Metrics.csv",
    base_dir / "Supplementary_Tables" / "TableS1_Pathway_Enrichment_Results.csv",
    base_dir / "Supplementary_Tables" / "TableS2_Complete_Regulatory_Architecture.csv",
    base_dir / "Supplementary_Tables" / "TableS3_Curated_50_GWAS_Loci.csv"
]

for table_path in tables_to_view:
    if table_path.exists():
        print(f"\n\n{'='*80}")
        print(f"📄 TABLE: {table_path.name}")
        print(f"{'='*80}")
        try:
            # Read and display the first 15 rows of each table
            df = pd.read_csv(table_path)
            print(df.head(15).to_string(index=False))
            if len(df) > 15:
                print(f"\n... (and {len(df) - 15} more rows)")
        except Exception as e:
            print(f"Error reading {table_path.name}: {e}")
    else:
        print(f"\n❌ Could not find: {table_path}")

print("\n")
