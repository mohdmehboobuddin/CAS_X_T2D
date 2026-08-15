import pandas as pd
import requests
import time
import os

# Using the highly stable Open Targets PLATFORM API (v4)
API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Ensembl IDs for your top targets
ensembl_map = {
    "JAZF1": "ENSG00000173273",
    "FTO": "ENSG00000140718",
    "GRB14": "ENSG00000144684",
    "CDKN2B": "ENSG00000147883",
    "VEGFA": "ENSG00000112715",
    "NOTCH2": "ENSG00000134250",
    "IRS1": "ENSG00000169047",
    "GCK": "ENSG00000106633"
}

# GraphQL query to extract DepMap CRISPR essentiality data
query = """
query targetEssentiality($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    depMapEssentiality {
      tissueSummary {
        tissueName
        screensCount
      }
    }
  }
}
"""

print("Executing Dynamic API Queries for CRISPR/DepMap Essentiality...\n")
results = []

for gene, ensembl_id in ensembl_map.items():
    variables = {"ensemblId": ensembl_id}
    
    try:
        response = requests.post(API_URL, json={"query": query, "variables": variables})
        if response.status_code == 200:
            data = response.json()
            
            essentiality_data = data['data']['target'].get('depMapEssentiality')
            
            if essentiality_data and essentiality_data.get('tissueSummary'):
                # Count total tissues where CRISPR knockout showed effect
                tissues_affected = len(essentiality_data['tissueSummary'])
                total_screens = sum([t['screensCount'] for t in essentiality_data['tissueSummary']])
                status = "Functionally Essential"
            else:
                tissues_affected = 0
                total_screens = 0
                status = "No CRISPR Effect Detected"
                
            results.append({
                "Target Gene": gene,
                "Ensembl ID": ensembl_id,
                "CRISPR DepMap Status": status,
                "Affected Tissues (Count)": tissues_affected,
                "Total CRISPR Screens": total_screens
            })
            print(f"Successfully retrieved CRISPR data for {gene}")
        else:
            print(f"API Error for {gene}: {response.status_code}")
            
    except Exception as e:
        print(f"Failed to query {gene}: {e}")
    
    time.sleep(0.5)  # Be polite to the server

# Save the programmatic results
df_results = pd.DataFrame(results)
out_path = "results/tables/automated_crispr_validation.csv"
os.makedirs("results/tables", exist_ok=True)
df_results.to_csv(out_path, index=False)

print("\n=== Automated CRISPR DepMap Validation ===")
print(df_results.to_string())
print(f"\nSaved automated CRISPR validation to: {out_path}")
