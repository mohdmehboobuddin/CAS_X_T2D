"""
01_data_processing.py
CAS-X Pipeline: Data Extraction & Tissue Filtering
"""
import pandas as pd
import os

def process_raw_gtex(input_path, output_path):
    print(f"Loading raw data from {input_path}...")
    if not os.path.exists(input_path):
        print(f"Warning: {input_path} not found. Please ensure raw data is placed in data/raw/")
        return
        
    df = pd.read_csv(input_path)
    
    # Filter for target metabolic tissues
    target_tissues = [
        'Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 
        'Adipose_Visceral_Omentum', 'Whole_Blood', 'Liver'
    ]
    df_met = df[df['tissue'].isin(target_tissues)].copy()
    
    # Clean and drop null gene symbols or variant IDs
    df_met = df_met.dropna(subset=['gene_symbol', 'variant_id', 'pip', 'afc'])
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_met.to_csv(output_path, index=False)
    print(f"Processed metabolic eQTL data saved to {output_path} ({len(df_met)} rows)")

if __name__ == "__main__":
    process_raw_gtex('../data/raw/broad_continuous_gtex_v11.csv', '../data/processed/clean_metabolic_gtex.csv')
