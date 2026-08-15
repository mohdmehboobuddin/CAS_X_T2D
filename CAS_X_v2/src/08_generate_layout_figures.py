import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
import seaborn as sns
import networkx as nx
import os

def generate_layout_figures():
    print("======================================================")
    print(" GENERATING LAYOUT FIGURES (FIG 1, 5, 6)")
    print("======================================================")
    
    os.makedirs("CAS_X_v2/results/figures", exist_ok=True)
    
    # ---------------------------------------------------------
    # FIGURE 1: CAS-X Framework Overview Flowchart
    # ---------------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(8, 10))
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 12)
    ax1.axis('off')
    
    steps = [
        "1. T2D GWAS Loci Integration\n(2,900+ Significant Loci)",
        "2. Multi-Tissue GTEx Extractor\n(6 Metabolic Tissues, including Liver)",
        "3. Continuous Feature Mapping\n(PIP * |aFC| Effect Sizes)",
        "4. Unsupervised PCA Engine\n(RobustScaler + Dynamic Variance)",
        "5. Systemic Prioritization Ranking\n(0-100 CAS-X Score)",
        "6. Biological & Clinical Validation\n(Pleiotropy & Tractability)"
    ]
    
    y_pos = np.linspace(10.5, 1.5, len(steps))
    
    for i, step in enumerate(steps):
        # Draw Box
        box = patches.FancyBboxPatch((2, y_pos[i]-0.5), 6, 1, boxstyle="round,pad=0.2", 
                                     edgecolor='black', facecolor='#e0f3db', lw=2)
        ax1.add_patch(box)
        # Add Text
        ax1.text(5, y_pos[i], step, ha='center', va='center', fontsize=11, fontweight='bold')
        
        # Draw Arrows (except for the last box)
        if i < len(steps) - 1:
            ax1.annotate('', xy=(5, y_pos[i+1]+0.7), xytext=(5, y_pos[i]-0.7),
                         arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=8))

    plt.title("CAS-X  Framework Architecture", fontweight='bold', fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/Fig1_Framework_Overview.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated Figure 1 (Framework Flowchart) at 600 DPI.")

    # ---------------------------------------------------------
    # FIGURE 5: Pleiotropy Network
    # ---------------------------------------------------------
    # New top targets from the 6-tissue run
    targets = ['AOC1', 'PSRC1', 'TENM2', 'UGT2B17', 'PSORS1C3', 'AMN', 'DHRS2', 'CELSR2', 'HPR']
    traits = ['Fasting Glucose', 'BMI / Obesity', 'LDL Cholesterol', 'Coronary Artery Disease', 'Triglycerides']
    
    G = nx.Graph()
    G.add_nodes_from(targets, bipartite=0)
    G.add_nodes_from(traits, bipartite=1)
    
    # Define metabolic relationships mapping to the new liver/adipose discoveries
    edges = [
        ('PSRC1', 'LDL Cholesterol'), ('PSRC1', 'Coronary Artery Disease'), 
        ('CELSR2', 'LDL Cholesterol'), ('CELSR2', 'Coronary Artery Disease'),
        ('AOC1', 'Fasting Glucose'), ('AOC1', 'BMI / Obesity'),
        ('TENM2', 'BMI / Obesity'), ('UGT2B17', 'Triglycerides'),
        ('HPR', 'LDL Cholesterol'), ('HPR', 'Triglycerides'),
        ('AMN', 'Fasting Glucose'), ('DHRS2', 'Fasting Glucose')
    ]
    G.add_edges_from(edges)
    
    fig5, ax5 = plt.subplots(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42, k=0.8)
    
    # Draw Traits (Squares, Pink)
    nx.draw_networkx_nodes(G, pos, nodelist=traits, node_color='#dd3497', 
                           node_shape='s', node_size=1500, edgecolors='black', ax=ax5)
    # Draw Genes (Circles, Blue)
    nx.draw_networkx_nodes(G, pos, nodelist=targets, node_color='#41b6c4', 
                           node_size=1200, edgecolors='black', ax=ax5)
    nx.draw_networkx_edges(G, pos, width=2, alpha=0.5, edge_color='gray', ax=ax5)
    
    # Custom Labels
    labels = {node: node for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, font_size=9, font_weight='bold', ax=ax5)
    
    plt.title('Curated Phenome-Wide Pleiotropy Network', fontweight='bold', fontsize=14)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/Fig5_Pleiotropy_Network.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated Figure 5 (Pleiotropy Network) at 600 DPI.")

    # ---------------------------------------------------------
    # FIGURE 6: Clinical Tractability Matrix
    # ---------------------------------------------------------
    # Using the top 15 genes from the new ranking
    top15 = ['AOC1', 'PSRC1', 'TENM2', 'UGT2B17', 'PSORS1C3', 'AMN', 'DHRS2', 
             'AMZ1', 'MIP', 'CELSR2', 'HPR', 'MLIP', 'RAB29', 'TEKT4P2', 'LINC00910']
    
    # 0 = Untractable, 1 = Theoretically Tractable, 2 = Highly Druggable
    tractability_data = {
        'Small Molecule (Pill)': [2, 1, 0, 2, 0, 1, 2, 0, 0, 1, 0, 0, 1, 0, 0],
        'Monoclonal Antibody':   [0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0],
        'PROTAC (Degradation)':  [1, 2, 1, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 1, 0]
    }
    
    df_tract = pd.DataFrame(tractability_data, index=top15)
    
    fig6, ax6 = plt.subplots(figsize=(6, 8))
    
    # Custom discrete colormap: Gray (0), Teal (1), Dark Blue (2)
    cmap_discrete = ListedColormap(['#e0e0e0', '#41b6c4', '#225ea8'])
    
    sns.heatmap(df_tract, cmap=cmap_discrete, cbar=False, linewidths=2, linecolor='white', ax=ax6)
    
    plt.title('Clinical Tractability of Top Targets', fontweight='bold', fontsize=14, pad=15)
    plt.xlabel('Therapeutic Modality', fontweight='bold', fontsize=12)
    plt.ylabel('Top 15 CAS-X Targets', fontweight='bold', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    
    # Create custom legend
    legend_elements = [
        patches.Patch(facecolor='#225ea8', edgecolor='black', label='Highly Druggable (Known Binding)'),
        patches.Patch(facecolor='#41b6c4', edgecolor='black', label='Theoretically Tractable'),
        patches.Patch(facecolor='#e0e0e0', edgecolor='black', label='Currently Untractable / Unknown')
    ]
    ax6.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1))
    
    plt.tight_layout()
    plt.savefig("CAS_X_v2/results/figures/Fig6_Clinical_Tractability.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("-> Generated Figure 6 (Tractability Matrix) at 600 DPI.")

    print("======================================================")

if __name__ == "__main__":
    generate_layout_figures()
