import pandas as pd
import numpy as np

def execute_max_signal_aggregation():
    file_path = 'CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv'
    target_tissues = [
        'Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 
        'Adipose_Visceral_Omentum', 'Whole_Blood', 'Liver'
    ]
    
    try:
        # Load and filter data
        df = pd.read_csv(file_path)
        df = df[df['tissue'].isin(target_tissues)].copy()
        
        # Calculate expected effect (PIP * |aFC|)
        df['expected_effect'] = df['pip'] * df['afc'].abs()
        
        # Find the MAXIMUM effect size for each gene across any tissue
        max_effects = df.groupby('gene_symbol')['expected_effect'].max().reset_index()
        
        # Get the primary tissue responsible for that max score
        idx = df.groupby('gene_symbol')['expected_effect'].idxmax()
        primary_tissues = df.loc[idx, ['gene_symbol', 'tissue']]
        
        # Merge and rank
        final_df = max_effects.merge(primary_tissues, on='gene_symbol')
        final_df['CASX_Score'] = final_df['expected_effect'].rank(pct=True) * 100
        
        # Sort and display
        final_df = final_df.sort_values('CASX_Score', ascending=False).reset_index(drop=True)
        final_df = final_df.rename(columns={'tissue': 'Primary_Tissue'})
        
        print("===== PHASE 2.7: MAXIMUM SIGNAL AGGREGATION =====")
        print(f"EVALUATED: {final_df.shape[0]} Genes")
        print("\nTOP 15 MASTER REGULATORS:")
        print(final_df[['gene_symbol', 'CASX_Score', 'Primary_Tissue']].head(15).to_string())
        print("=================================================")
        
    except Exception as e:
        print(f"Execution Error: {e}")

if __name__ == "__main__":
    execute_max_signal_aggregation()
