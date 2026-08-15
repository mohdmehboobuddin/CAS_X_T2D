import pandas as pd
from pathlib import Path

def generate_table_s3():
    # Look for the actual populated data in your processed folder
    potential_sources = [
        "data/processed/master_t2d_loci_50.csv",
        "data/processed/t2d_50_loci_candidate_panel.csv",
        "data/processed/master_t2d_loci.csv"
    ]

    df = None
    source_used = ""
    
    for src in potential_sources:
        try:
            temp_df = pd.read_csv(src)
            # Check if the dataframe actually has data
            if not temp_df.empty and len(temp_df) > 5:
                df = temp_df
                source_used = src
                break
        except FileNotFoundError:
            continue

    if df is not None:
        out_path = Path("results/Supplementary_Tables/TableS3_Curated_50_GWAS_Loci.csv")
        
        # Ensure the directory exists just in case
        out_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save the populated table
        df.to_csv(out_path, index=False)
        print(f"✅ Success! Found populated data in '{source_used}'")
        print(f"✅ Table S3 has been generated at: {out_path}\n")
        
        print("Here is a preview of your recovered Table S3:")
        print("="*80)
        # Display a clean preview of the first 5 rows
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        print(df.head().to_string(index=False))
        print("="*80)
    else:
        print("❌ Could not find a populated dataset. Running the original GWAS collection script...")
        import os
        os.system("python scripts/01_collect_gwas.py")
        print("Please re-run this script (python scripts/69_generate_table_s3.py) after the collection finishes.")

if __name__ == "__main__":
    generate_table_s3()
