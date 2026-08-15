import pandas as pd

# CAS-X genes
casx = pd.read_csv(
    "results/tables/casx_v4_rankings.csv"
)

casx_genes = set(casx["GENE"])

# Open Targets genes
with open(
    "data/processed/opentargets_t2d_core_genes.txt"
) as f:

    ot_genes = set(
        line.strip()
        for line in f
        if line.strip()
    )

# Overlap
overlap = casx_genes.intersection(
    ot_genes
)

print("\n=== Open Targets Validation ===\n")

print(
    "CAS-X genes:",
    len(casx_genes)
)

print(
    "Open Targets genes:",
    len(ot_genes)
)

print(
    "Overlap:",
    len(overlap)
)

print(
    "\nShared genes:"
)

for g in sorted(overlap):
    print(g)

# Save overlap table

pd.DataFrame(
    {"GENE": sorted(overlap)}
).to_csv(
    "results/tables/opentargets_overlap.csv",
    index=False
)

print("\nSaved:")
print(
    "results/tables/opentargets_overlap.csv"
)
