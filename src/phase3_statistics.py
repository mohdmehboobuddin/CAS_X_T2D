import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score

def bootstrap_auc(y_true, y_scores, n_bootstraps=1000, seed=42):
    rng = np.random.RandomState(seed)
    bootstrapped_scores = []
    
    for i in range(n_bootstraps):
        indices = rng.randint(0, len(y_scores), len(y_scores))
        if len(np.unique(y_true[indices])) < 2:
            continue
        score = roc_auc_score(y_true[indices], y_scores[indices])
        bootstrapped_scores.append(score)
        
    sorted_scores = np.array(bootstrapped_scores)
    sorted_scores.sort()
    
    ci_lower = sorted_scores[int(0.025 * len(sorted_scores))]
    ci_upper = sorted_scores[int(0.975 * len(sorted_scores))]
    
    return ci_lower, ci_upper

def run_true_statistical_benchmark():
    file_path = 'CAS_X_v2/data/processed/broad_continuous_gtex_v11.csv'
    # THE FIX: Pointing to the correct main data directory
    validation_path = 'data/processed/opentargets_t2d_core_genes.txt'
    
    try:
        # 1. Load GTEx Data and Calculate Expected Effect
        df = pd.read_csv(file_path)
        df['expected_effect'] = df['pip'] * df['afc'].abs()
        
        target_tissues = ['Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 
                          'Adipose_Visceral_Omentum', 'Whole_Blood', 'Liver']
        
        # Calculate CAS-X Scores (Multi-Tissue Max)
        df_casx = df[df['tissue'].isin(target_tissues)]
        casx_scores = df_casx.groupby('gene_symbol')['expected_effect'].max()
        
        # Calculate Pancreas-Only Scores
        df_panc = df[df['tissue'] == 'Pancreas']
        panc_scores = df_panc.groupby('gene_symbol')['expected_effect'].max()
        
        # Combine into a clean dataframe
        compare_df = pd.DataFrame({'CASX': casx_scores, 'Pancreas': panc_scores}).fillna(0)
        compare_df = compare_df.reset_index()
        
        # 2. Read the OpenTargets Validation List
        with open(validation_path, 'r') as f:
            true_genes = [line.strip() for line in f if line.strip() and line.strip().lower() not in ['gene', 'symbol', 'target', 'gene_symbol']]
            
        # 3. Assign True Labels
        compare_df['True_Label'] = compare_df['gene_symbol'].isin(true_genes).astype(int)
        
        target_count = compare_df['True_Label'].sum()
        if target_count == 0:
            print("CRITICAL ERROR: Still found 0 matching targets. The gene names in the .txt file don't match the GTEx format.")
            return
            
        # 4. Calculate AUROCs
        casx_auc = roc_auc_score(compare_df['True_Label'], compare_df['CASX'])
        panc_auc = roc_auc_score(compare_df['True_Label'], compare_df['Pancreas'])
        
        # 5. Bootstrap Confidence Intervals
        casx_ci = bootstrap_auc(compare_df['True_Label'].values, compare_df['CASX'].values)
        panc_ci = bootstrap_auc(compare_df['True_Label'].values, compare_df['Pancreas'].values)
        
        print("===== PHASE 3: OPEN TARGETS STATISTICAL BENCHMARKING =====")
        print(f"Total GTEx Genes Evaluated: {len(compare_df)}")
        print(f"Open Targets Found in GTEx: {target_count}")
        print("\nPERFORMANCE METRICS:")
        print(f"CAS-X (6-Tissue) AUROC:   {casx_auc:.3f} (95% CI: {casx_ci[0]:.3f} - {casx_ci[1]:.3f})")
        print(f"Pancreas-Only AUROC:      {panc_auc:.3f} (95% CI: {panc_ci[0]:.3f} - {panc_ci[1]:.3f})")
        
        print("\nSTATISTICAL PROOF:")
        if casx_ci[0] > panc_ci[1]:
            print("SUCCESS: CAS-X confidence interval strictly dominates Pancreas-only.")
        elif casx_auc > panc_auc:
            print("MODERATE SUCCESS: CAS-X outperforms the baseline, but confidence intervals overlap.")
        else:
            print("FAILURE: CAS-X does not outperform the baseline.")
        print("==========================================================")
        
    except Exception as e:
        print(f"Execution Error: {e}")

if __name__ == "__main__":
    run_true_statistical_benchmark()
