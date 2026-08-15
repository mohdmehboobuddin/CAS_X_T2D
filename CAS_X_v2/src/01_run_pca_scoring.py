import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
import os

def run_casx_pca(input_path, output_path):
    print("======================================================")
    print(" INITIATING CAS-X v2.1 MATHEMATICAL ENGINE (SIGNAL-BASED)")
    print("======================================================")

    # 1. Load data
    df = pd.read_csv(input_path)
    
    # 2. Calculate the Biological Association Signal: -log10(P) * abs(slope)
    df['abs_slope'] = df['slope'].abs()
    df['neg_log10_p'] = -np.log10(df['pval_nominal'] + 1e-300) # prevent log(0)
    df['eqtl_signal'] = df['neg_log10_p'] * df['abs_slope']
    
    # Find the variant with the strongest signal for each gene-tissue pair
    idx_max = df.groupby(['gene_symbol', 'tissue'])['eqtl_signal'].idxmax()
    best_signals = df.loc[idx_max, ['gene_symbol', 'tissue', 'eqtl_signal']]
    
    # Pivot to wide format
    matrix = best_signals.pivot(index='gene_symbol', columns='tissue', values='eqtl_signal')
    
    # Identify the 5 core metabolic tissues
    core_tissues = [
        'Pancreas', 
        'Muscle_Skeletal', 
        'Adipose_Subcutaneous', 
        'Adipose_Visceral_Omentum', 
        'Whole_Blood'
    ]
    
    # Filter for core tissues and fill missing with 0.0 (0.0 signal = no confident effect)
    matrix = matrix.reindex(columns=core_tissues).fillna(0.0)
    
    print(f"Constructed Gene x Tissue Signal matrix: {matrix.shape[0]} genes across {matrix.shape[1]} tissues.")

    # 3. Standardization (Z-score within tissue)
    scaler = StandardScaler()
    Z_matrix = scaler.fit_transform(matrix)

    # 4. Unsupervised PCA
    pca = PCA()
    PC_scores = pca.fit_transform(Z_matrix)
    explained_variance = pca.explained_variance_ratio_

    # 5. Component Selection
    K = 5 
    print(f"Retaining all {K} principal components.")

    # 6. Sign-Invariant Variance-Weighted Scoring
    raw_scores = []
    for g_idx in range(len(matrix)):
        gene_pc_scores = PC_scores[g_idx, :K]
        weighted_sq_sum = np.sum(explained_variance[:K] * (gene_pc_scores ** 2))
        raw_scores.append(np.sqrt(weighted_sq_sum))
    
    # 7. Final Min-Max Normalization (0-100)
    min_max = MinMaxScaler(feature_range=(0, 100))
    final_scores = min_max.fit_transform(np.array(raw_scores).reshape(-1, 1)).flatten()

    # Combine into final dataframe (FIXED THE BUG HERE)
    final_df = pd.DataFrame({
        'Target_Gene': matrix.index,
        'Raw_CASX_Score': raw_scores,
        'CASX_Score': final_scores
    })
    
    final_df = final_df.sort_values(by='CASX_Score', ascending=False).reset_index(drop=True)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    final_df.to_csv(output_path, index=False)
    
    print("\n--- TOP 10 PRIORITIZED TARGETS ---")
    print(final_df[['Target_Gene', 'CASX_Score']].head(10).to_string(index=False))
    print("======================================================")

if __name__ == "__main__":
    INPUT_PATH = "CAS_X_v2/data/processed/casx_v6_continuous_gtex.csv"
    OUTPUT_PATH = "CAS_X_v2/results/tables/CASX_Prioritization_Scores.csv"
    run_casx_pca(INPUT_PATH, OUTPUT_PATH)
