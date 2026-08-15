import pandas as pd
import os
import urllib.request
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
RANKINGS_FILE = os.path.join(PROJECT_ROOT, "data", "processed", "casx_v6_probabilistic_rankings.csv")

casx_v6 = pd.read_csv(RANKINGS_FILE)
top_15 = casx_v6.head(15)['GENE'].tolist()

print(f"Requesting STRING PPI Network for Top 15 Genes: {', '.join(top_15)}")

string_api_url = "https://version-11-5.string-db.org/api/json/network"
params = urllib.parse.urlencode({
    "identifiers": "%0d".join(top_15), 
    "species": 9606, # Human
    "caller_identity": "CASX_Framework" 
}).encode("utf-8")

try:
    req = urllib.request.Request(string_api_url, data=params)
    with urllib.request.urlopen(req) as response:
        network = json.loads(response.read().decode("utf-8"))
        
    print("\n--- KNOWN BIOLOGICAL INTERACTIONS FOUND ---")
    interactions_found = 0
    for edge in network:
        # Filter for high-confidence interactions (score > 0.4)
        if edge['score'] > 0.4:
            print(f"{edge['preferredName_A']} <---> {edge['preferredName_B']} (Confidence: {edge['score']:.2f})")
            interactions_found += 1
            
    if interactions_found == 0:
        print("No high-confidence interactions found among these specific targets.")
        
except Exception as e:
    print(f"Failed to fetch PPI network: {e}")
