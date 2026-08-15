import pandas as pd
import numpy as np

def execute_threshold_aggregation():
    file_path = 'CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv'
    target_tissues = [
        'Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 
        'Adipose_Visceral_Omentum', 'Whole_Blood', 'Liver'
    ]
    
    try:
        df = pd.read_csv(file_path)
        df = df[df['tissue'].isin(target_tissues)].copy()
        
        # PIP-weighted Expected Effect Size
        df['expected_effect'] = df['pip'] * df['afc'].abs()
        
        agg_df = df.groupby(['gene_symbol', 'tissue'])['expected_effect'].max().reset_index()
        wide_df = agg_df.pivot(index='gene_symbol', columns='tissue', values='expected_effect')
        
        # THE FIX: Minimum Evidence Threshold (Must be in >= 3 tissues)
        tissue_counts = wide_df.notna().sum(axis=1)
        filtered_df = wide_df[tissue_counts >= 3]
        
        # Rank within the filtered, stable set
        rank_df = filtered_df.rank(pct=True, numeric_only=True)
        casx_raw = rank_df.mean(axis=1) * 100
        casx_final = (casx_raw - casx_raw.min()) / (casx_raw.max() - casx_raw.min()) * 100
        
        final_scores = pd.DataFrame({
            'CASX_Score': casx_final,
            'Tissue_Count': tissue_counts[tissue_counts >= 3]
        }).sort_values('CASX_Score', ascending=False)
        
        print("===== PHASE 2.6: THRESHOLD-PROTECTED PIP AGGREGATION =====")
        print(f"MATRIX: {filtered_df.shape[0]} Stable Genes (>=3 Tissues)")
        print("\nTOP 15 TARGETS:")
        print(final_scores.head(15).to_string())
        print("==========================================================")
        
    except Exception as e:
        print(f"Execution Error: {e}")

if __name__ == "__main__":
    execute_threshold_aggregation()
