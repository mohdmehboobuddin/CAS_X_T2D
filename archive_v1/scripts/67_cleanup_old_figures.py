import os
import glob

# Setup Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
OLD_FIG_DIR = os.path.join(PROJECT_ROOT, "results", "figures")
OLD_TAB_DIR = os.path.join(PROJECT_ROOT, "results", "manuscript_tables")
PUB_DIR = os.path.join(PROJECT_ROOT, "results", "v6_publication_assets")

print("======================================================")
print("  CLEANING UP OLD V1-V5 FIGURES AND TABLES")
print("======================================================\n")

# 1. Clear old figures
old_figs = glob.glob(os.path.join(OLD_FIG_DIR, "*"))
for f in old_figs:
    if os.path.isfile(f):
        os.remove(f)
print(f"Deleted {len(old_figs)} outdated files from results/figures/")

# 2. Clear old manuscript tables
old_tabs = glob.glob(os.path.join(OLD_TAB_DIR, "*"))
for f in old_tabs:
    if os.path.isfile(f):
        os.remove(f)
print(f"Deleted {len(old_tabs)} outdated files from results/manuscript_tables/")

# 3. Remove draft Figure 1s from the publication folder
draft_figs = ["Figure1_CASX_Workflow_v6.png", "Figure1_Polished_Workflow_v6.png"]
deleted_drafts = 0
for draft in draft_figs:
    draft_path = os.path.join(PUB_DIR, draft)
    if os.path.exists(draft_path):
        os.remove(draft_path)
        deleted_drafts += 1
        
print(f"Deleted {deleted_drafts} intermediate drafts from v6_publication_assets/")

print("\n======================================================")
print("SUCCESS: Repository figures are now completely clean!")
print("======================================================\n")
