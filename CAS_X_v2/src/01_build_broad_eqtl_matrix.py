import pandas as pd
import os
import glob

def extract_broad_gtex_data():
    print("======================================================")
    print(" INITIATING LARGE-SCALE GTEx v11 EXTRACTION (6-TISSUE)")
    print("======================================================")

    candidates_path = "CAS_X_v2/data/processed/broad_t2d_candidates.csv"
    candidates = pd.read_csv(candidates_path)['Target_Gene'].tolist()
    
    # ADDED LIVER TO CORE TISSUES
    core_tissues = {
        'Pancreas': 'Pancreas',
        'Muscle_Skeletal': 'Muscle_Skeletal',
        'Adipose_Subcutaneous': 'Adipose_Subcutaneous',
        'Adipose_Visceral_Omentum': 'Adipose_Visceral_Omentum',
        'Whole_Blood': 'Whole_Blood',
        'Liver': 'Liver'
    }

    gtex_dir = "data/raw/eqtl/extracted/GTEx_Analysis_v11_eQTL/"
    extracted_data = []

    for standard_name, file_keyword in core_tissues.items():
        search_pattern = f"{gtex_dir}*{file_keyword}*.parquet"
        files = glob.glob(search_pattern)
        
        if not files:
            print(f"WARNING: Could not find data for {standard_name}")
            continue
            
        file_path = files[0]
        print(f"Processing {standard_name}...")
        
        try:
            df = pd.read_parquet(file_path, columns=['gene_name', 'variant_id', 'pip', 'afc'])
            df_filtered = df[df['gene_name'].isin(candidates)].copy()
            df_filtered['tissue'] = standard_name
            df_filtered = df_filtered.rename(columns={'gene_name': 'gene_symbol'})
            df_filtered = df_filtered.dropna(subset=['pip', 'afc'])
            extracted_data.append(df_filtered)
            print(f" -> Found {len(df_filtered)} fine-mapped eQTL associations")
        except Exception as e:
            print(f"ERROR reading {file_path}: {e}")

    final_eqtl_df = pd.concat(extracted_data, ignore_index=True)
    output_path = "CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv"
    final_eqtl_df.to_csv(output_path, index=False)
    
    print("\n======================================================")
    print(f"SUCCESS: Built multi-tissue matrix with {len(final_eqtl_df)} total regulatory signals.")
    print("======================================================")

if __name__ == "__main__":
    extract_broad_gtex_data()
