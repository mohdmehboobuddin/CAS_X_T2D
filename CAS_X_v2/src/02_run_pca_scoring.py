import pandas as pd
import numpy as np
from sklearn.preprocessing import RobustScaler, MinMaxScaler
from sklearn.decomposition import PCA
import os

def run_casx_pca(input_path, output_path):
    print("======================================================")
    print(" INITIATING CAS-X v2.3 (ROBUST SCALER + DYNAMIC PCA)")
    print("======================================================")

    df = pd.read_csv(input_path)
    
    df['abs_afc'] = df['afc'].abs()
    df['eqtl_signal'] = df['pip'] * df['abs_afc']
    
    idx_max = df.groupby(['gene_symbol', 'tissue'])['eqtl_signal'].idxmax()
    best_signals = df.loc[idx_max, ['gene_symbol', 'tissue', 'eqtl_signal']]
    
    matrix = best_signals.pivot(index='gene_symbol', columns='tissue', values='eqtl_signal')
    
    # 6-TISSUE ARCHITECTURE
    core_tissues = ['Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Whole_Blood', 'Liver']
    matrix = matrix.reindex(columns=core_tissues).fillna(0.0)
    
    print(f"Gene x Tissue Signal matrix: {matrix.shape[0]} genes, {matrix.shape[1]} tissues.")

    # UPGRADE 1: ROBUST SCALER
    scaler = RobustScaler()
    Z_matrix = scaler.fit_transform(matrix)

    pca = PCA()
    PC_scores = pca.fit_transform(Z_matrix)
    explained_variance = pca.explained_variance_ratio_

    # UPGRADE 2: DYNAMIC COMPONENT SELECTION (Targeting >= 85% variance)
    cumulative_variance = np.cumsum(explained_variance)
    K = np.argmax(cumulative_variance >= 0.85) + 1
    print(f"Retaining {K} principal components (explaining {cumulative_variance[K-1]*100:.1f}% of variance).")

    raw_scores = []
    for g_idx in range(len(matrix)):
        gene_pc_scores = PC_scores[g_idx, :K]
        weighted_sq_sum = np.sum(explained_variance[:K] * (gene_pc_scores ** 2))
        raw_scores.append(np.sqrt(weighted_sq_sum))
    
    min_max = MinMaxScaler(feature_range=(0, 100))
    final_scores = min_max.fit_transform(np.array(raw_scores).reshape(-1, 1)).flatten()

    final_df = pd.DataFrame({
        'Target_Gene': matrix.index,
        'Raw_CASX_Score': raw_scores,
        'CASX_Score': final_scores
    }).sort_values(by='CASX_Score', ascending=False).reset_index(drop=True)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    final_df.to_csv(output_path, index=False)
    
    print("\n--- DE NOVO DISCOVERY: TOP 15 NOVEL TARGETS ---")
    print(final_df[['Target_Gene', 'CASX_Score']].head(15).to_string(index=False))
    print("======================================================")

if __name__ == "__main__":
    INPUT_PATH = "CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv"
    OUTPUT_PATH = "CAS_X_v2/results/tables/Broad_CASX_Prioritization_Scores.csv"
    run_casx_pca(INPUT_PATH, OUTPUT_PATH)
