"""
Script 70: Generate PCA Scree Plot
Calculates the variance explained by the Principal Component Analysis (PCA)
using the continuous multi-tissue eQTL dataset. Generates a high-resolution
(600 DPI) scree plot for Supplementary Figure S6 to validate dimensionality reduction.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from pathlib import Path

def generate_scree_plot():
    data_path = Path("data/processed/casx_v6_continuous_gtex.csv")
    
    if not data_path.exists():
        print(f"Error: Dataset not found at {data_path}")
        return
        
    df = pd.read_csv(data_path)
    
    # Isolate numeric tissue features, dropping identifiers and GWAS scores
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    features = [c for c in numeric_cols if 'GWAS' not in c.upper() and 'POS' not in c.upper() and 'CHR' not in c.upper()]
    
    if not features:
        print("Error: No numeric features found.")
        return
        
    X = df[features].fillna(0)
    
    # Standardize data and fit PCA
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA()
    pca.fit(X_scaled)
    
    variance_explained = pca.explained_variance_ratio_ * 100
    cumulative_variance = np.cumsum(variance_explained)
    
    # Generate publication-quality 600 DPI plot
    plt.figure(figsize=(10, 6), dpi=600)
    sns.set_theme(style="whitegrid")
    
    x_labels = [f"PC{i+1}" for i in range(len(variance_explained))]
    
    bars = plt.bar(x_labels, variance_explained, color='#2c3e50', alpha=0.8, label='Individual Variance')
    plt.plot(x_labels, cumulative_variance, marker='o', color='#e74c3c', linewidth=2, markersize=8, label='Cumulative Variance')
    
    plt.title('PCA Scree Plot: Multi-Tissue eQTL Variance', fontsize=14, pad=15, fontweight='bold')
    plt.xlabel('Principal Components', fontsize=12)
    plt.ylabel('Percentage of Variance Explained (%)', fontsize=12)
    plt.ylim(0, 115) 
    
    plt.legend(loc='lower right', framealpha=0.9)
    
    # Add offset data labels to prevent overlapping
    for i, bar in enumerate(bars):
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval - 4, f'{yval:.1f}%', ha='center', va='top', fontsize=11, color='white', fontweight='bold')
        
    for i, cum_val in enumerate(cumulative_variance):
        plt.text(i, cum_val + 3, f'{cum_val:.1f}%', ha='center', va='bottom', fontsize=11, color='#c0392b', fontweight='bold')

    plt.tight_layout()
    
    out_dir = Path("results/Supplementary_Figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "FigS6_PCA_Scree_Plot.png"
    
    plt.savefig(out_path)
    print(f"Generated 600 DPI Scree Plot: {out_path}")

if __name__ == "__main__":
    generate_scree_plot()
