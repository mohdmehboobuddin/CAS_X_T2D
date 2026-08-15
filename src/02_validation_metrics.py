"""
02_validation_metrics.py
CAS-X Framework:  Specificity and Benchmarking
Runs statistical tests against negative control tissues (e.g., Brain Cortex).
"""
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

def run_specificity_benchmark(metabolic_data, control_data):
    print("Running Specificity Validation against Empirical Negative Controls...")
    
    # Extract expected effect arrays
    metabolic_scores = metabolic_data['expected_effect'].dropna().values
    control_scores = control_data['expected_effect'].dropna().values
    
    # Perform Mann-Whitney U rank test
    stat, p_value = mannwhitneyu(metabolic_scores, control_scores, alternative='greater')
    
    print(f"--- Benchmark Results ---")
    print(f"Metabolic targets (n={len(metabolic_scores)}) Mean Effect: {np.mean(metabolic_scores):.4f}")
    print(f"Control targets (n={len(control_scores)}) Mean Effect: {np.mean(control_scores):.4f}")
    print(f"Mann-Whitney U P-value: {p_value:.2e}")
    
    return p_value

# (In the repository, you would call this function using your loaded dataframes)
