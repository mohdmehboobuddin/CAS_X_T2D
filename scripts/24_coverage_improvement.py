import matplotlib.pyplot as plt

labels = ["Pancreas Only", "Multi-Tissue"]
values = [8, 22]

plt.figure(figsize=(7,5))

bars = plt.bar(labels, values)

for bar, val in zip(bars, values):
    plt.text(
        bar.get_x() + bar.get_width()/2,
        val + 0.3,
        str(val),
        ha="center"
    )

plt.ylabel("Supported Genes")
plt.title("Multi-Tissue Integration Improves CAS-X Gene Coverage")

plt.tight_layout()

plt.savefig(
    "results/figures/Figure2_Coverage_Improvement.png",
    dpi=450,
    bbox_inches="tight"
)

print("Saved Figure2")
