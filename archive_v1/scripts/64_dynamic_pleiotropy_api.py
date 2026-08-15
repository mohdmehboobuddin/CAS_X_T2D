import pandas as pd
import requests
import time
import os

# Define the GraphQL endpoint
API_URL = "https://api.genetics.opentargets.org/graphql"

# We map your top genes to their Ensembl IDs for the API
ensembl_map = {
    "FTO": "ENSG00000140718",
    "IRS1": "ENSG00000169047",
    "JAZF1": "ENSG00000173273",
    "CDKN2B": "ENSG00000147883",
    "VEGFA": "ENSG00000112715",
    "GRB14": "ENSG00000144684",
    "NOTCH2": "ENSG00000134250",
    "KLF14": "ENSG00000170454",
    "GCK": "ENSG00000106633",
    "ADCY5": "ENSG00000173175"
}

# The GraphQL Query Structure
query = """
query targetDiseases($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    knownDrugs {
      count
    }
  }
  associatedDiseases(ensemblId: $ensemblId) {
    count
    rows {
      disease {
        name
      }
      score
    }
  }
}
"""

print("Executing Dynamic API Queries to Open Targets Genetics...\n")
results = []

for gene, ensembl_id in ensembl_map.items():
    variables = {"ensemblId": ensembl_id}
    
    try:
        response = requests.post(API_URL, json={"query": query, "variables": variables})
        if response.status_code == 200:
            data = response.json()
            
            # Extract top 3 disease associations based on Open Targets association score
            diseases = data['data']['associatedDiseases']['rows']
            top_traits = [d['disease']['name'] for d in diseases[:3]] if diseases else ["No data"]
            
            results.append({
                "Target Gene": gene,
                "Ensembl ID": ensembl_id,
                "Top Automated API Traits": " | ".join(top_traits)
            })
            print(f"Successfully retrieved data for {gene}")
        else:
            print(f"API Error for {gene}: {response.status_code}")
            
    except Exception as e:
        print(f"Failed to query {gene}: {e}")
    
    # Polite API delay
    time.sleep(0.5)

# Save the programmatic results
df_results = pd.DataFrame(results)
out_path = "results/tables/automated_pleiotropy_network.csv"
os.makedirs("results/tables", exist_ok=True)
df_results.to_csv(out_path, index=False)

print("\n=== Programmatic Pleiotropy Results ===")
print(df_results)
print(f"\nSaved automated network to: {out_path}")
