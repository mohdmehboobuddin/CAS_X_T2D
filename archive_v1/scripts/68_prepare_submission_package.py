import os
import shutil
from pathlib import Path

def create_submission_package():
    # 1. Define base paths relative to the project root (CAS_X_T2D)
    project_root = Path.cwd()
    assets_dir = project_root / "results" / "v6_publication_assets"
    data_dir = project_root / "data" / "processed"
    
    # 2. Define the new submission package directories
    submission_dir = project_root / "results" / "submission_package"
    folders = {
        "main_figs": submission_dir / "Main_Manuscript_Figures",
        "main_tabs": submission_dir / "Main_Manuscript_Tables",
        "supp_figs": submission_dir / "Supplementary_Figures",
        "supp_tabs": submission_dir / "Supplementary_Tables"
    }
    
    # Create the directories
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=True)
        print(f"Created directory: {folder.relative_to(project_root)}")

    # 3. Define the file mapping (Source File -> (Destination Folder, New Prefix))
    file_mapping = {
        # --- MAIN MANUSCRIPT FIGURES ---
        assets_dir / "Figure1_CASX_Framework_Clean.png": (folders["main_figs"], "Fig1A_"),
        assets_dir / "Figure2_Coverage_Improvement_v6.png": (folders["main_figs"], "Fig1B_"),
        assets_dir / "Figure3_Multitissue_Heatmap_v6.png": (folders["main_figs"], "Fig2_"),
        assets_dir / "Figure4_Top15_CASX_v6.png": (folders["main_figs"], "Fig3_"),
        assets_dir / "LOTO_Sensitivity_Analysis.png": (folders["main_figs"], "Fig4A_"),
        assets_dir / "Negative_Control_Benchmark.png": (folders["main_figs"], "Fig4B_"),
        assets_dir / "Figure11_Pleiotropy_Network.png": (folders["main_figs"], "Fig5_"),
        assets_dir / "Figure9_Clinical_Tractability_Matrix.png": (folders["main_figs"], "Fig6_"),
        
        # --- MAIN MANUSCRIPT TABLES ---
        assets_dir / "Table2_Top15_CASX_v6.csv": (folders["main_tabs"], "Table1_"),
        assets_dir / "Table5_Clinical_Tractability.csv": (folders["main_tabs"], "Table2_"),
        
        # --- SUPPLEMENTARY FIGURES ---
        assets_dir / "Figure5_Tissue_Contribution_v6.png": (folders["supp_figs"], "FigS1_"),
        assets_dir / "Figure6_Validation_Summary_v6.png": (folders["supp_figs"], "FigS2_"),
        assets_dir / "Figure7_External_Validation_Landscape_v6.png": (folders["supp_figs"], "FigS3_"),
        assets_dir / "Figure8_Pathway_Enrichment_v6.png": (folders["supp_figs"], "FigS4_"),
        assets_dir / "Figure10_Regulatory_Distribution.png": (folders["supp_figs"], "FigS5_"),
        
        # --- SUPPLEMENTARY TABLES ---
        assets_dir / "Table4_Pathway_Enrichment_v6.csv": (folders["supp_tabs"], "TableS1_"),
        assets_dir / "Table6_Regulatory_Architecture.csv": (folders["supp_tabs"], "TableS2_"),
        data_dir / "casx_50_loci_dataset.csv": (folders["supp_tabs"], "TableS3_")
    }

    # 4. Execute the copy and rename operation
    print("\nStarting file organization...")
    for source_path, (dest_folder, prefix) in file_mapping.items():
        if source_path.exists():
            # Construct new filename and path
            new_filename = f"{prefix}{source_path.name}"
            dest_path = dest_folder / new_filename
            
            # Copy file
            shutil.copy2(source_path, dest_path)
            print(f"Copied: {source_path.name} -> {dest_path.relative_to(project_root)}")
        else:
            print(f"WARNING: File not found - {source_path.name}")

    print("\nSuccess! All submission assets have been organized in 'results/submission_package/'")

if __name__ == "__main__":
    create_submission_package()
