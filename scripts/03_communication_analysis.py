import scanpy as sc
import squidpy as sq
import pandas as pd
import os

print("--- [SCRIPT 3/3] STARTING: Cell-Cell Communication Analysis ---")

# --- Configuration ---
BASE_DIR = os.getcwd()
RESULTS_DIR = os.path.join(BASE_DIR, 'results/')
PROCESSED_DATA_DIR = os.path.join(RESULTS_DIR, 'processed_data/')
TABLES_DIR = os.path.join(RESULTS_DIR, 'communication_tables/')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures/')
for path in [TABLES_DIR, FIGURES_DIR]: os.makedirs(path, exist_ok=True)
sc.settings.figdir = FIGURES_DIR

# --- 1. Load Data and Annotate ---
print("\nStep 1/3: Loading data and creating annotations...")
adata_full = sc.read_h5ad(os.path.join(PROCESSED_DATA_DIR, "final_annotated_adata.h5ad"))
# Create a simple numerical annotation for Squidpy
adata_full.obs['cell_type_id'] = 'Cluster ' + adata_full.obs['leiden_integrated'].astype(str)
adata_full.obs['cell_type_id'] = adata_full.obs['cell_type_id'].astype('category')

with open(os.path.join(RESULTS_DIR, 'target_cluster_id.txt'), 'r') as f:
    TARGET_CLUSTER_ID = f.read().strip()
print(f"  - Target cluster identified from previous step: {TARGET_CLUSTER_ID}")

# --- 2. Run Squidpy Analysis ---
print("\nStep 2/3: Running Squidpy communication analysis (this may take over an hour)...")
# Subsample for computational feasibility
sc.pp.subsample(adata_full, n_obs=50000, random_state=42)
sq.gr.ligrec(adata_full, n_perms=1000, cluster_key="cell_type_id", use_raw=False, copy=False)
ligrec_results_df = adata_full.uns['cell_type_id_ligrec']
ligrec_results_df.to_csv(os.path.join(TABLES_DIR, 'squidpy_full_ligrec_results.csv'))

# --- 3. Save Key Visuals and Tables ---
print("\nStep 3/3: Saving final plots and tables...")
sq.pl.ligrec(adata_full, cluster_key="cell_type_id", pvalue_threshold=0.05, show=False, save='_ligrec_circle_plot.png')
# Filter for significant interactions targeting our cluster
ligrec_flat = ligrec_results_df.copy()
ligrec_flat.columns = ['_'.join(col).strip() for col in ligrec_flat.columns.values]
ligrec_flat.rename(columns={ligrec_flat.columns[0]: 'source', ligrec_flat.columns[1]: 'target'}, inplace=True)
pvalue_col = ligrec_flat.columns[2]
ligrec_flat[pvalue_col] = pd.to_numeric(ligrec_flat[pvalue_col], errors='coerce')
target_interactions = ligrec_flat[
    (ligrec_flat['target'] == f"Cluster {TARGET_CLUSTER_ID}") &
    (ligrec_flat[pvalue_col] < 0.05)
].sort_values(pvalue_col)
target_interactions.to_csv(os.path.join(TABLES_DIR, f'interactions_targeting_cluster_{TARGET_CLUSTER_ID}.csv'))
print(f"  - Found {len(target_interactions)} significant interactions targeting Cluster {TARGET_CLUSTER_ID}.")

print("--- [SCRIPT 3/3] COMPLETE ---")
