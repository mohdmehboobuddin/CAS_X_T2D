import pandas as pd

def run_diagnostics():
    file_path = 'CAS_X_v2/data/processed/casx_v6_continuous_gtex.csv'
    
    try:
        df = pd.read_csv(file_path)
        
        print("===== PHASE 1: MATRIX & MISSINGNESS DIAGNOSTICS =====")
        print(f"\n1. MATRIX SHAPE: {df.shape[0]} rows x {df.shape[1]} columns")
        
        print("\n2. COLUMN NAMES:")
        print(list(df.columns))
        
        print("\n3. MISSING DATA (NaN) OR ZERO COUNTS PER COLUMN:")
        # Check for NaNs
        nans = df.isna().sum()
        # Check for exact zeros (since you mentioned zero-imputation might have already happened)
        zeros = (df == 0).sum()
        
        summary_df = pd.DataFrame({'NaN_Count': nans, 'Zero_Count': zeros})
        print(summary_df.to_string())
        
        print("\n4. DATA GLIMPSE (FIRST 3 ROWS):")
        print(df.head(3).to_string())
        print("=====================================================")
        
    except Exception as e:
        print(f"Error reading the file: {e}")

if __name__ == "__main__":
    run_diagnostics()
