import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# 1. Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

GTEX_FILE = os.path.join(PROCESSED_DIR, "casx_v6_continuous_gtex.csv")
MASTER_FILE = os.path.join(PROCESSED_DIR, "casx_master_dataset.csv")

print("Loading datasets...")
gtex_df = pd.read_csv(GTEX_FILE)
master_df = pd.read_csv(MASTER_FILE)

# 2. Aggregate the continuous eQTL data per gene
print("Aggregating continuous eQTL statistics...")
# Calculate absolute effect sizes to measure magnitude, regardless of direction
gtex_df['abs_slope'] = gtex_df['slope'].abs()

# Group by gene to extract the most significant signals
eqtl_agg = gtex_df.groupby('gene_symbol').agg(
    eQTL_Tissue_Count=('tissue', 'nunique'),
    eQTL_Max_Abs_Slope=('abs_slope', 'max'),
    eQTL_Min_Pval=('pval_nominal', 'min')
).reset_index()

# Convert p-value to -log10 scale for proper weighting (smaller p-value = higher score)
eqtl_agg['eQTL_NegLog10_Pval'] = -np.log10(eqtl_agg['eQTL_Min_Pval'] + 1e-300) # prevent log(0)

# 3. Merge with Master Dataset
# Assuming the master dataset uses 'GENE' as the column name
merged_df = pd.merge(master_df, eqtl_agg, left_on='GENE', right_on='gene_symbol', how='left')

# Drop duplicate gene entries to clean the output
merged_df = merged_df.drop_duplicates(subset=['GENE']).copy()

# Fill NaN values for genes that had no eQTLs
merged_df['eQTL_Tissue_Count'] = merged_df['eQTL_Tissue_Count'].fillna(0)
merged_df['eQTL_Max_Abs_Slope'] = merged_df['eQTL_Max_Abs_Slope'].fillna(0)
merged_df['eQTL_NegLog10_Pval'] = merged_df['eQTL_NegLog10_Pval'].fillna(0)

# Clean GWAS P-value if it exists
if 'P_VALUE' in merged_df.columns:
    merged_df['GWAS_NegLog10_Pval'] = pd.to_numeric(merged_df['P_VALUE'], errors='coerce')
    merged_df['GWAS_NegLog10_Pval'] = -np.log10(merged_df['GWAS_NegLog10_Pval'] + 1e-300)
    merged_df['GWAS_NegLog10_Pval'] = merged_df['GWAS_NegLog10_Pval'].fillna(0)
else:
    merged_df['GWAS_NegLog10_Pval'] = 1.0 # Fallback

# 4. Prepare Features for the Statistical Model
# We select continuous, data-driven features rather than manual categorical points
features = ['eQTL_Tissue_Count', 'eQTL_Max_Abs_Slope', 'eQTL_NegLog10_Pval', 'GWAS_NegLog10_Pval']

# Add expression score if it exists in the master dataset
if 'EXPRESSION_SCORE' in merged_df.columns:
    features.append('EXPRESSION_SCORE')

print(f"\nTraining unsupervised PCA model on features: {features}")
X = merged_df[features].values

# Standardize features (mean=0, variance=1) so they are weighted fairly
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 5. Apply Principal Component Analysis (Optimization Basis)
# PCA finds the mathematical combination of features that explains the most variance
pca = PCA(n_components=1)
merged_df['CASX_v6_Raw_Score'] = pca.fit_transform(X_scaled)

# Normalize the score to a clean 0 to 100 scale for easy reading
min_score = merged_df['CASX_v6_Raw_Score'].min()
max_score = merged_df['CASX_v6_Raw_Score'].max()
merged_df['CASX_v6_Final_Score'] = ((merged_df['CASX_v6_Raw_Score'] - min_score) / (max_score - min_score)) * 100

# Rank the genes
final_df = merged_df.sort_values(by='CASX_v6_Final_Score', ascending=False).reset_index(drop=True)

# 6. Save the results
output_path = os.path.join(PROCESSED_DIR, "casx_v6_probabilistic_rankings.csv")
cols_to_keep = ['GENE', 'CASX_v6_Final_Score'] + features
final_df[cols_to_keep].to_csv(output_path, index=False)

# Format pandas output to remove scientific notation (e.g., 100.00 instead of 1.00e+02)
pd.options.display.float_format = '{:.2f}'.format

print("\n--- NEW CAS-X RANKINGS ---")
print(final_df[['GENE', 'CASX_v6_Final_Score']].head(10).to_string(index=False))
print(f"\nSaved mathematically optimized scores to: {output_path}")
