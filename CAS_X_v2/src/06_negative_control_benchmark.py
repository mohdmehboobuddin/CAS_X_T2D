import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
import os
import glob

def run_negative_control():
    print("======================================================")
    print(" INITIATING TRUE EMPIRICAL NEGATIVE CONTROL BENCHMARK")
    print("======================================================")

    # 1. Load Candidates and Open Targets
    candidates = pd.read_csv("CAS_X_v2/data/processed/broad_t2d_candidates.csv")['Target_Gene'].tolist()
    with open("data/processed/opentargets_t2d_core_genes.txt", "r") as f:
        ot_genes = set(line.strip() for line in f if line.strip())

    # 2. Define Non-Metabolic Negative Control Tissues
    neg_tissues = {
        'Brain_Cortex': 'Brain_Cortex',
        'Skin': 'Skin_Not_Sun_Exposed_Suprapubic',
        'Nerve': 'Nerve_Tibial',
        'Pituitary': 'Pituitary',
        'Spleen': 'Spleen'
    }

    gtex_dir = "data/raw/eqtl/extracted/GTEx_Analysis_v11_eQTL/"
    extracted_data = []

    for standard_name, file_keyword in neg_tissues.items():
        search_pattern = f"{gtex_dir}*{file_keyword}*.parquet"
        files = glob.glob(search_pattern)
        
        if files:
            print(f"Extracting {standard_name}...")
            df = pd.read_parquet(files[0], columns=['gene_name', 'variant_id', 'pip', 'afc'])
            df_filtered = df[df['gene_name'].isin(candidates)].copy()
            df_filtered['tissue'] = standard_name
            df_filtered = df_filtered.rename(columns={'gene_name': 'gene_symbol'})
            df_filtered = df_filtered.dropna(subset=['pip', 'afc'])
            extracted_data.append(df_filtered)

    neg_eqtl = pd.concat(extracted_data, ignore_index=True)
    
    # 3. Process Negative Control Matrix
    neg_eqtl['eqtl_signal'] = neg_eqtl['pip'] * neg_eqtl['afc'].abs()
    idx_max = neg_eqtl.groupby(['gene_symbol', 'tissue'])['eqtl_signal'].idxmax()
    matrix = neg_eqtl.loc[idx_max].pivot(index='gene_symbol', columns='tissue', values='eqtl_signal').fillna(0.0)
    matrix = matrix.reindex(columns=list(neg_tissues.keys())).fillna(0.0)

    # 4. Run CAS-X PCA on Negative Data
    scaler = StandardScaler()
    Z_matrix = scaler.fit_transform(matrix)
    pca = PCA()
    PC_scores = pca.fit_transform(Z_matrix)
    explained_variance = pca.explained_variance_ratio_
    
    raw_scores = []
    for g_idx in range(len(matrix)):
        weighted_sq_sum = np.sum(explained_variance[:5] * (PC_scores[g_idx, :5] ** 2))
        raw_scores.append(np.sqrt(weighted_sq_sum))
        
    final_scores = MinMaxScaler(feature_range=(0, 100)).fit_transform(np.array(raw_scores).reshape(-1, 1)).flatten()
    
    df_eval = pd.DataFrame({'Target_Gene': matrix.index, 'Neg_CASX': final_scores})
    df_eval['Label'] = df_eval['Target_Gene'].apply(lambda x: 1 if x in ot_genes else 0)
    
    neg_auc = roc_auc_score(df_eval['Label'], df_eval['Neg_CASX'])
    print(f"\n-> Negative Control Framework AUC: {neg_auc:.3f}")
    
    # 5. Plot Figure 4B (600 DPI)
    # Hardcoded metabolic AUC from our previous robust validation run for comparison
    metabolic_auc = 0.679 
    
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = sns.barplot(
        x=['Metabolic Tissues\n(Pancreas, Muscle, Adipose, Blood)', 'Negative Control Tissues\n(Brain, Skin, Nerve, Pituitary)'], 
        y=[metabolic_auc, neg_auc], 
        palette=['#2ca25f', '#de2d26'], 
        ax=ax
    )
    
    ax.axhline(0.50, color='gray', linestyle='--', label='Random Expectation (AUC=0.50)')
    ax.set_ylim([0.0, 0.8])
    ax.set_ylabel('Predictive Power (AUROC)', fontweight='bold')
    ax.set_title('Empirical Negative Control Benchmark', fontweight='bold', pad=15)
    ax.legend(loc='upper right')
    
    for i, v in enumerate([metabolic_auc, neg_auc]):
        ax.text(i, v + 0.02, f"AUC = {v:.3f}", color='black', ha='center', fontweight='bold')

    os.makedirs("CAS_X_v2/results/figures", exist_ok=True)
    plt.savefig("CAS_X_v2/results/figures/Fig4B_Negative_Control_Benchmark.png", dpi=600, bbox_inches='tight')
    print("SUCCESS: Saved Negative Control Plot (Fig 4B) at 600 DPI.")
    print("======================================================")

if __name__ == "__main__":
    run_negative_control()
