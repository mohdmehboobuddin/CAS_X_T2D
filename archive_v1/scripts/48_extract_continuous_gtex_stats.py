import pandas as pd
import os
import urllib.request
import urllib.parse
import json
import time

# 1. Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
CANDIDATE_GENES_FILE = os.path.join(PROCESSED_DIR, "casx_master_dataset.csv")

TISSUES = [
    "Pancreas", 
    "Adipose_Subcutaneous", 
    "Adipose_Visceral_Omentum", 
    "Muscle_Skeletal", 
    "Whole_Blood"
]

print(f"Loading master dataset from: {CANDIDATE_GENES_FILE}")
candidates_df = pd.read_csv(CANDIDATE_GENES_FILE)
gene_col = 'GENE'
raw_genes = set(candidates_df[gene_col].dropna().astype(str).tolist())
print(f"Loaded {len(raw_genes)} candidate genes.")

# 2. Map Gene Symbols to GTEx Versioned GENCODE IDs using GTEx's own API
print("\nMapping genes to GTEx Versioned GENCODE IDs...")
gencode_ids = {}

for gene in raw_genes:
    search_gene = gene.split('.')[0] 
    url = f"https://gtexportal.org/api/v2/reference/gene?geneId={search_gene}"
    
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            if 'data' in data and len(data['data']) > 0:
                # Grab the exact versioned ID GTEx expects
                gtex_gene_info = data['data'][0]
                gencode_id = gtex_gene_info.get('gencodeId')
                if gencode_id:
                    gencode_ids[gene] = gencode_id
        time.sleep(0.1) # Be polite to the server
    except Exception as e:
        print(f"  Warning: Could not map {gene}: {e}")

print(f"Successfully mapped {len(gencode_ids)} genes.")

# 3. Query the GTEx API for significant eQTLs
extracted_data = []
# FIX: Restored the correct v2 endpoint for single tissue eQTLs
base_url = "https://gtexportal.org/api/v2/association/singleTissueEqtl"

print("\nQuerying GTEx API for continuous eQTL statistics...")
for tissue in TISSUES:
    print(f"  Fetching data for {tissue}...")
    
    for original_gene, gencode_id in gencode_ids.items():
        query_params = urllib.parse.urlencode({
            'datasetId': 'gtex_v8',
            'tissueSiteDetailId': tissue,
            'gencodeId': gencode_id
        })
        api_url = f"{base_url}?{query_params}"
        
        try:
            req = urllib.request.Request(api_url)
            with urllib.request.urlopen(req) as response:
                data = json.loads(response.read().decode('utf-8'))
                
                if 'data' in data and len(data['data']) > 0:
                    for record in data['data']:
                        extracted_data.append({
                            'tissue': tissue,
                            'gene_symbol': original_gene,
                            'gencode_id': gencode_id,
                            'variant_id': record.get('variantId'),
                            'pval_nominal': record.get('pValue', record.get('pvalue')),
                            'slope': record.get('nes') # Normalized Effect Size
                        })
            time.sleep(0.1) 
        except urllib.error.HTTPError as e:
            # Catch HTTP errors silently to avoid cluttering the terminal, 
            # as a 404 or 400 often just means no eQTLs exist for that gene in that tissue.
            pass 
        except Exception as e:
            print(f"    Failed to fetch {original_gene} in {tissue}: {e}")

# 4. Save the continuous data
output_df = pd.DataFrame(extracted_data)
output_path = os.path.join(PROCESSED_DIR, "casx_v6_continuous_gtex.csv")

if not output_df.empty:
    output_df.to_csv(output_path, index=False)
    print(f"\nSuccess! Extracted {len(output_df)} continuous eQTL records to {output_path}")
else:
    print("\nExtraction returned 0 records. The genes may not have significant eQTLs in these tissues.")
