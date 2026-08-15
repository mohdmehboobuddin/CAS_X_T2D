"""
01_scoring_engine.py
CAS-X Framework:  Target Prioritization
Calculates the Expected Regulatory Effect (PIP * |aFC|) across metabolic tissues.
"""
import pandas as pd
import numpy as np
import os

def calculate_casx_scores(input_path, output_path):
    print("Loading GTEx v11 continuous data...")
    df = pd.read_csv(input_path)
    
    # Target metabolic tissues
    target_tissues = ['Pancreas', 'Muscle_Skeletal', 'Adipose_Subcutaneous', 
                      'Adipose_Visceral_Omentum', 'Whole_Blood', 'Liver']
    df_met = df[df['tissue'].isin(target_tissues)].copy()
    
    # Phase 2.7 Math: Expected Effect = Posterior Inclusion Probability * Absolute aFC
    print("Computing Expected Regulatory Effect...")
    df_met['expected_effect'] = df_met['pip'] * df_met['afc'].abs()
    
    # Aggregate max effect per gene across all metabolic tissues
    casx_scores = df_met.groupby('gene_symbol')['expected_effect'].max().reset_index()
    
    # Normalize to 0-100 scale for final ranking
    min_score = casx_scores['expected_effect'].min()
    max_score = casx_scores['expected_effect'].max()
    casx_scores['casx_normalized_score'] = ((casx_scores['expected_effect'] - min_score) / (max_score - min_score)) * 100
    
    # Sort and save
    casx_scores = casx_scores.sort_values(by='casx_normalized_score', ascending=False)
    casx_scores.to_csv(output_path, index=False)
    print(f"Scoring complete. Master rankings saved to {output_path}")

if __name__ == "__main__":
    os.makedirs('../data/processed', exist_ok=True)
    calculate_casx_scores('../data/raw/broad_continuous_gtex_v11.csv', '../data/processed/casx_master_rankings.csv')
