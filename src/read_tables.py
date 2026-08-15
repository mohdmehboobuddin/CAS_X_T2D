import pandas as pd
import os
import glob

def extract_local_tables(directory_path="."):
    print("----- BEGIN TABLE DATA EXTRACT -----")
    
    # Process all CSV files in the folder
    csv_files = glob.glob(os.path.join(directory_path, "*.csv"))
    for file in csv_files:
        print(f"\n[DATA FROM: {os.path.basename(file)}]")
        try:
            df = pd.read_csv(file)
            print(df.to_csv(index=False))
        except Exception as e:
            print(f"Error reading {file}: {e}")
            
    # Process all Excel files in the folder
    excel_files = glob.glob(os.path.join(directory_path, "*.xlsx"))
    for file in excel_files:
        print(f"\n[DATA FROM: {os.path.basename(file)}]")
        try:
            df = pd.read_excel(file)
            print(df.to_csv(index=False))
        except Exception as e:
            print(f"Error reading {file}: {e}")
            
    print("----- END TABLE DATA EXTRACT -----")

if __name__ == "__main__":
    # If the script is in the same folder as the tables, leave this as "."
    # Otherwise, replace "." with your absolute folder path (e.g., "/home/mehboob/cas_x/tables")
    extract_local_tables(".")
