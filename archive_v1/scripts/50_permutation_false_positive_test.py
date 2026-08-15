import pandas as pd
import numpy as np
import os

# 1. Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

RANKINGS_FILE = os.path.join(PROCESSED_DIR, "casx_v6_probabilistic_rankings.csv")
GTEX_FILE = os.path.join(PROCESSED_DIR, "casx_v6_continuous_gtex.csv")

print("Loading continuous data for Permutation Testing...")
rankings_df = pd.read_csv(RANKINGS_FILE)
gtex_df = pd.read_csv(GTEX_FILE)

# Calculate absolute effect sizes
gtex_df['abs_slope'] = gtex_df['slope'].abs()

# 2. Calculate the REAL continuous signal strength of our Top Genes
# We will use the sum of their absolute eQTL slopes as a metric of signal magnitude
real_top_10_genes = rankings_df.head(10)['GENE'].tolist()
real_agg = gtex_df.groupby('gene_symbol').agg(Total_Effect_Size=('abs_slope', 'sum'))

# What is the average effect size mass of our real top 10 genes?
real_top_10_signal = real_agg.loc[real_agg.index.intersection(real_top_10_genes), 'Total_Effect_Size'].mean()

print(f"Real Top 10 CAS-X Genes Average eQTL Effect Mass: {real_top_10_signal:.2f}")

# 3. Run 1,000 Permutations 
n_permutations = 1000
print(f"\nRunning {n_permutations} permutations to calculate continuous False Discovery Rate (FDR)...")

random_top_10_signals = []

for i in range(n_permutations):
    # Shuffle the gene labels, breaking the biological association with the eQTLs
    shuffled_df = gtex_df.copy()
    shuffled_df['gene_symbol'] = np.random.permutation(shuffled_df['gene_symbol'].values)
    
    # Recalculate the signal strength for this random universe
    shuffled_agg = shuffled_df.groupby('gene_symbol').agg(Total_Effect_Size=('abs_slope', 'sum'))
    
    # Find the top 10 strongest genes in this random universe
    random_top_10 = shuffled_agg.sort_values(by='Total_Effect_Size', ascending=False).head(10)
    random_top_10_signals.append(random_top_10['Total_Effect_Size'].mean())

# 4. Calculate Empirical P-Value
# How many times did random noise create a top 10 list stronger than our real biology?
times_noise_won = sum(1 for random_signal in random_top_10_signals if random_signal >= real_top_10_signal)
empirical_p_value = times_noise_won / n_permutations

print("\n--- CONTINUOUS PERMUTATION TEST RESULTS ---")
print(f"Average Top 10 Signal in Random Simulation: {np.mean(random_top_10_signals):.2f}")
print(f"Empirical P-Value (FDR): {empirical_p_value:.4f}")

if empirical_p_value < 0.05:
    print("\nCONCLUSION: SUCCESS! The CAS-X signal strength is statistically significant (P < 0.05).")
    print("This directly refutes the reviewer's claim that your results are just background noise.")
else:
    print("\nCONCLUSION: The signal may still be influenced by background noise.")
