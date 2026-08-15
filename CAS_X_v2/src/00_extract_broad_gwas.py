import pandas as pd
import zipfile
import os

def extract_all_t2d_genes():
    print("======================================================")
    print(" INITIATING DE NOVO GWAS DISCOVERY MODULE")
    print("======================================================")
    
    gwas_zip_path = "data/raw/gwas/gwas_catalog_v1.0.2-associations_e115_r2026-06-01_full.zip"
    
    if not os.path.exists(gwas_zip_path):
        print(f"ERROR: Cannot find GWAS catalog at {gwas_zip_path}")
        return

    print("Loading full GWAS Catalog (this may take a moment)...")
    
    # Read the TSV directly from the ZIP file
    with zipfile.ZipFile(gwas_zip_path, 'r') as z:
        # Assuming there is only one file in the zip, get its name
        file_name = z.namelist()[0]
        with z.open(file_name) as f:
            # The GWAS catalog is tab-separated and can have mixed dtypes
            df = pd.read_csv(f, sep='\t', low_memory=False)

    print(f"Total genome-wide associations loaded: {len(df)}")
    
    # Filter for Type 2 Diabetes associations
    # Using str.contains to catch variations like "Type 2 diabetes mellitus"
    t2d_mask = df['DISEASE/TRAIT'].str.contains('Type 2 diabetes', case=False, na=False)
    t2d_df = df[t2d_mask].copy()
    print(f"Total T2D associations found: {len(t2d_df)}")
    
    # Filter for genome-wide significance (P < 5x10^-8)
    # Convert P-VALUE to numeric, forcing errors to NaN
    t2d_df['P-VALUE'] = pd.to_numeric(t2d_df['P-VALUE'], errors='coerce')
    sig_df = t2d_df[t2d_df['P-VALUE'] < 5e-8].copy()
    
    # Extract mapped genes, split multiple genes (separated by commas or hyphens), and clean
    genes = set()
    for gene_str in sig_df['MAPPED_GENE'].dropna():
        # Split by comma or hyphen which are common in the GWAS catalog
        for g in gene_str.replace('-', ',').split(','):
            clean_g = g.strip()
            if clean_g != '' and clean_g != 'NR': # Ignore empty and 'NR' (Not Reported)
                genes.add(clean_g)

    print(f"\nSUCCESS: Extracted {len(genes)} unique genome-wide significant candidate genes.")
    
    # Save the expanded candidate list
    output_path = "CAS_X_v2/data/processed/broad_t2d_candidates.csv"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    pd.DataFrame({'Target_Gene': sorted(list(genes))}).to_csv(output_path, index=False)
    
    print(f"Expanded candidate list saved to {output_path}")
    print("======================================================")

if __name__ == "__main__":
    extract_all_t2d_genes()
