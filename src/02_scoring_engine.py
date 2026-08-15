"""
02_scoring_engine.py
CAS-X Pipeline: Phase 2.7 Target Prioritization Scoring
"""
import pandas as pd
import numpy as np
import os

def run_scoring_engine(input_path, output_path):
    print("Running CAS-X Scoring Engine...")
    df = pd.read_csv(input_path)
    
    # Phase 2.7 Math: Expected Regulatory Effect
    df['expected_effect'] = df['pip'] * df['afc'].abs()
    
    # Group by gene symbol to find maximum regulatory impact across tissues
    gene_scores = df.groupby('gene_symbol')['expected_effect'].max().reset_index()
    
    # Normalize scores to 0 - 100 range
    min_val = gene_scores['expected_effect'].min()
    max_val = gene_scores['expected_effect'].max()
    gene_scores['casx_score'] = ((gene_scores['expected_effect'] - min_val) / (max_val - min_val)) * 100
    
    # Sort descending
    gene_scores = gene_scores.sort_values(by='casx_score', ascending=False)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    gene_scores.to_csv(output_path, index=False)
    print(f"Scoring complete. Saved master rankings to {output_path}")

if __name__ == "__main__":
    run_scoring_engine('../data/processed/clean_metabolic_gtex.csv', '../data/processed/casx_master_rankings.csv')
