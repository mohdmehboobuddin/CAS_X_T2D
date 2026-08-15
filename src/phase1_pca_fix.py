import pandas as pd
import numpy as np
from sklearn.preprocessing import RobustScaler
from sklearn.decomposition import PCA
from scipy.stats import spearmanr

def run_pca_diagnostics():
    file_path = 'CAS_X_v2/data/processed/casx_v6_continuous_gtex.csv'
    
    try:
        # 1. Load Data and get Absolute Effect Sizes
        df = pd.read_csv(file_path)
        df['abs_slope'] = df['slope'].abs()
        
        # 2. Extract max effect size for each gene per tissue
        agg_df = df.groupby(['gene_symbol', 'tissue'])['abs_slope'].max().reset_index()
        
        # 3. Pivot to Matrix: Rows = Genes, Columns = Tissues
        wide_df = agg_df.pivot(index='gene_symbol', columns='tissue', values='abs_slope')
        
        print("===== PHASE 1: PCA ARCHITECTURE AND SENSITIVITY =====")
        print(f"\n1. CORRECTED MATRIX SHAPE: {wide_df.shape[0]} Genes (rows) x {wide_df.shape[1]} Tissues (columns)")
        
        print("\n2. MISSING VALUES (NaNs) AFTER PIVOTING:")
        print(wide_df.isna().sum().to_string())
        
        # --- MODEL A: ZERO IMPUTATION (Original Pipeline) ---
        X_zero = wide_df.fillna(0)
        scaler = RobustScaler()
        X_zero_scaled = scaler.fit_transform(X_zero)
        
        pca_zero = PCA()
        pca_zero_scores = pca_zero.fit_transform(X_zero_scaled)
        
        print("\n3. TRUE PCA VARIANCE (Zero Imputation):")
        variances = pca_zero.explained_variance_ratio_ * 100
        for i, var in enumerate(variances):
            print(f"   PC{i+1}: {var:.1f}%")
        print(f"   Cumulative PC1-PC3: {sum(variances[:3]):.1f}%")
        
        # --- MODEL B: SENSITIVITY CHECK (Median Imputation) ---
        # This tests if the rankings fall apart if we don't use zero-imputation
        X_med = wide_df.fillna(wide_df.median())
        X_med_scaled = scaler.fit_transform(X_med)
        
        pca_med = PCA()
        pca_med_scores = pca_med.fit_transform(X_med_scaled)
        
        # 4. Calculate Scores properly using absolute transformed gene coordinates
        # Formula: Sum of (|PC Score| * Explained Variance)
        casx_zero = np.sum(np.abs(pca_zero_scores) * pca_zero.explained_variance_ratio_, axis=1)
        casx_med = np.sum(np.abs(pca_med_scores) * pca_med.explained_variance_ratio_, axis=1)
        
        # Scale 0 to 100
        score_zero_scaled = (casx_zero - casx_zero.min()) / (casx_zero.max() - casx_zero.min()) * 100
        score_med_scaled = (casx_med - casx_med.min()) / (casx_med.max() - casx_med.min()) * 100
        
        # 5. Compare the two models
        corr, p_val = spearmanr(score_zero_scaled, score_med_scaled)
        
        print("\n4. MISSING-DATA SENSITIVITY ANALYSIS:")
        print(f"   Spearman Rank Correlation (Zero vs Median Imputation): {corr:.4f}")
        if corr > 0.85:
            print("   -> PASSED: Your model is highly robust to how missing data is handled!")
        else:
            print("   -> FAILED: Rankings are highly sensitive to zero-imputation. We must change the imputation method.")
        
        print("=====================================================")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    run_pca_diagnostics()
