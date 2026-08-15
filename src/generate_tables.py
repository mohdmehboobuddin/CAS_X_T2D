import os
import pandas as pd

# Create output directory for main manuscript tables
os.makedirs('results/main_tables', exist_ok=True)

print("Generating Main Manuscript Tables...")

# Table 1: Target Prioritization and Scoring Summary
table1_data = {
    'Gene Symbol': [
        'JAZF1', 'FTO', 'GRB14', 'CDKN2B', 'ZBED3', 'VEGFA', 'ZFAND6', 
        'NOTCH2', 'KLF14', 'IRS1', 'GCK', 'ADCY5', 'SUGP1', 'PROX1', 'CDKAL1'
    ],
    'CAS-X Score': [
        100.00, 79.31, 70.81, 67.99, 67.27, 62.79, 60.56, 59.70, 
        59.04, 58.72, 54.83, 53.36, 52.07, 50.79, 50.37
    ],
    'Tissue Count': [5.0, 5.0, 4.0, 3.0, 5.0, 5.0, 3.0, 4.0, 3.0, 3.0, 2.0, 3.0, 1.0, 2.0, 3.0],
    'Max eQTL Effect (|aFC|)': [
        0.638, 0.839, 0.886, 0.650, 0.475, 0.529, 0.773, 0.548, 
        0.729, 0.626, 0.789, 0.745, 0.446, 0.874, 0.685
    ],
    'GWAS -log10(P)': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
}
t1_df = pd.DataFrame(table1_data)
t1_df.to_csv('results/main_tables/Table_1_Target_Prioritization_Summary.csv', index=False)
print("-> Saved Table 1 to results/main_tables/Table_1_Target_Prioritization_Summary.csv")

# Table 2: Clinical Tractability and Therapeutic Modalities
table2_data = {
    'Target Gene': [
        'JAZF1', 'FTO', 'GRB14', 'CDKN2B', 'ZBED3', 'VEGFA', 'ZFAND6', 
        'NOTCH2', 'KLF14', 'IRS1', 'GCK', 'ADCY5', 'SUGP1', 'PROX1', 'CDKAL1'
    ],
    'Small Molecule (Pill)': [0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0],
    'Monoclonal Antibody': [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'PROTAC (Degradation)': [1.0, 1.0, 0.0, 1.0, 0.5, 0.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0]
}
t2_df = pd.DataFrame(table2_data)
t2_df.to_csv('results/main_tables/Table_2_Clinical_Tractability.csv', index=False)
print("-> Saved Table 2 to results/main_tables/Table_2_Clinical_Tractability.csv")

print("All main tables successfully generated!")
