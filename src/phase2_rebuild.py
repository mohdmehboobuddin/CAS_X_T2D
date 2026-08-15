import pandas as pd
import numpy as np

def execute_systemic_aggregation():
    file_path = 'CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv'
    target_tissues = [
        'Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 
        'Adipose_Visceral_Omentum', 'Whole_Blood', 'Liver'
    ]
    
    try:
        df = pd.read_csv(file_path)
        df = df[df['tissue'].isin(target_tissues)].copy()
        
        # 1. PIP-weighted Expected Effect Size
        df['expected_effect'] = df['pip'] * df['afc'].abs()
        
        agg_df = df.groupby(['gene_symbol', 'tissue'])['expected_effect'].max().reset_index()
        wide_df = agg_df.pivot(index='gene_symbol', columns='tissue', values='expected_effect')
        
        # 2. Base Percentile Rank (Ignores NaNs)
        rank_df = wide_df.rank(pct=True, numeric_only=True)
        base_score = rank_df.mean(axis=1)
        
        # 3. The Systemic Multiplier (Rewards Pan-Metabolic Pleiotropy)
        tissue_counts = wide_df.notna().sum(axis=1)
        systemic_multiplier = tissue_counts / len(target_tissues)
        
        # 4. Final Systemic CAS-X Score
        casx_raw = base_score * systemic_multiplier * 100
        
        casx_final = (casx_raw - casx_raw.min()) / (casx_raw.max() - casx_raw.min()) * 100
        
        final_scores = pd.DataFrame({
            'CASX_Score': casx_final,
            'Tissue_Count': tissue_counts
        }).sort_values('CASX_Score', ascending=False)
        
        print("===== PHASE 2.5: SYSTEMIC PIP-WEIGHTED AGGREGATION =====")
        print(f"MATRIX: {wide_df.shape[0]} Genes x {wide_df.shape[1]} Tissues")
        print("\nTOP 15 SYSTEMIC TARGETS:")
        print(final_scores.head(15).to_string())
        print("========================================================")
        
    except Exception as e:
        print(f"Execution Error: {e}")

if __name__ == "__main__":
    execute_systemic_aggregation()
