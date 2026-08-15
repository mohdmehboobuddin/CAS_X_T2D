import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
import os

def evaluate_baselines():
    print("======================================================")
    print(" RUNNING ROBUST BASELINE COMPARISON & VALIDATION")
    print("======================================================")

    # 1. Load Data
    raw_eqtl = pd.read_csv("CAS_X_v2/data/processed/casx_v6_continuous_gtex.csv")
    casx_scores = pd.read_csv("CAS_X_v2/results/tables/CASX_Prioritization_Scores.csv")
    
    # Load Open Targets Core Genes
    try:
        with open("data/processed/opentargets_t2d_core_genes.txt", "r") as f:
            ot_genes = set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        print("ERROR: Could not find Open Targets validation list.")
        return

    # 2. Pivot raw eQTL to wide format for baselines
    raw_eqtl['abs_slope'] = raw_eqtl['slope'].abs()
    raw_eqtl['neg_log10_p'] = -np.log10(raw_eqtl['pval_nominal'] + 1e-300)
    raw_eqtl['eqtl_signal'] = raw_eqtl['neg_log10_p'] * raw_eqtl['abs_slope']
    
    idx_max = raw_eqtl.groupby(['gene_symbol', 'tissue'])['eqtl_signal'].idxmax()
    best_signals = raw_eqtl.loc[idx_max, ['gene_symbol', 'tissue', 'eqtl_signal']]
    matrix = best_signals.pivot(index='gene_symbol', columns='tissue', values='eqtl_signal').fillna(0.0)
    
    # 3. Create evaluation dataframe
    df = casx_scores[['Target_Gene', 'CASX_Score']].copy()
    df['Label'] = df['Target_Gene'].apply(lambda x: 1 if x in ot_genes else 0)
    
    # 4. Calculate Biological Baselines
    abs_matrix = matrix.abs()
    df['Baseline_Max'] = df['Target_Gene'].map(abs_matrix.max(axis=1))
    df['Baseline_Mean'] = df['Target_Gene'].map(abs_matrix.mean(axis=1))
    
    if 'Pancreas' in abs_matrix.columns:
        df['Baseline_Pancreas'] = df['Target_Gene'].map(abs_matrix['Pancreas'])
    else:
        df['Baseline_Pancreas'] = 0.0

    # 5. Evaluate Metrics
    results = []
    methods = {
        'CAS-X (Multi-Tissue PCA)': 'CASX_Score',
        'Max eQTL Signal': 'Baseline_Max',
        'Mean eQTL Signal': 'Baseline_Mean',
        'Pancreas-Only eQTL Signal': 'Baseline_Pancreas'
    }

    baseline_prev = df['Label'].mean()

    for name, col in methods.items():
        auroc = roc_auc_score(df['Label'], df[col])
        auprc = average_precision_score(df['Label'], df[col])
        results.append({'Method': name, 'AUROC': round(auroc, 3), 'AUPRC': round(auprc, 3)})

    # Calculate True Random Expectation over 10,000 iterations
    np.random.seed(42)
    rand_aurocs, rand_auprcs = [], []
    for _ in range(10000):
        rand_scores = np.random.uniform(0, 100, len(df))
        rand_aurocs.append(roc_auc_score(df['Label'], rand_scores))
        rand_auprcs.append(average_precision_score(df['Label'], rand_scores))
        
    results.append({
        'Method': 'Random Expectation (10k iter)', 
        'AUROC': round(np.mean(rand_aurocs), 3), 
        'AUPRC': round(np.mean(rand_auprcs), 3)
    })

    res_df = pd.DataFrame(results).sort_values(by='AUROC', ascending=False)
    
    print(res_df.to_string(index=False))
    
    os.makedirs("CAS_X_v2/results/tables", exist_ok=True)
    res_df.to_csv("CAS_X_v2/results/tables/Validation_Metrics.csv", index=False)
    print("======================================================")

if __name__ == "__main__":
    evaluate_baselines()
