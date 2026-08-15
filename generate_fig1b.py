import os
import matplotlib.pyplot as plt
import seaborn as sns

# Create output directory for TIFF assets
os.makedirs('results/tiff_figures', exist_ok=True)

# Set global publication styling (600 DPI, LZW TIFF)
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'Arial', 'font.family': 'sans-serif'})

# Corrected candidate gene counts with eQTL support
count_pancreas = 12
count_multi = 28
pct_increase = ((count_multi - count_pancreas) / count_pancreas) * 100

print(f"Pancreas-only target count: {count_pancreas}")
print(f"Multi-tissue target count: {count_multi}")
print(f"Computed Coverage Increase: +{pct_increase:.1f}%")

# Generate Figure 1B Bar Chart
fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
categories = ['Single-Tissue\n(Pancreas Only)', 'CAS-X Multi-Tissue\n(5 GTEx Tissues)']
counts = [count_pancreas, count_multi]

bars = ax.bar(categories, counts, color=['#4c72b0', '#55a868'], width=0.55, edgecolor='black')
ax.set_title("Improvement in Target Coverage via Multi-Tissue Integration", fontsize=13, fontweight='bold', pad=15)
ax.set_ylabel("Number of Target Genes with eQTL Support", fontsize=11, fontweight='bold')
ax.set_ylim(0, 35)

for bar in bars:
    y = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2.0, y + 0.8, f"{int(y)}", ha='center', fontweight='bold', fontsize=12)

# Annotate percentage increase
ax.annotate(
    f"+{pct_increase:.1f}% Coverage Increase", 
    xy=(0.5, counts[1] + 1.5), 
    xytext=(0.5, counts[1] + 2.5),
    ha='center', 
    va='bottom', 
    fontsize=13, 
    fontweight='bold', 
    color='#c0392b'
)

plt.tight_layout()
plt.savefig(
    'results/tiff_figures/Figure_1B_Coverage_Improvement.tiff', 
    format='tiff', 
    dpi=600, 
    pil_kwargs={"compression": "tiff_lzw"}
)
plt.close()

print("Figure 1B successfully generated and saved to results/tiff_figures/Figure_1B_Coverage_Improvement.tiff")
