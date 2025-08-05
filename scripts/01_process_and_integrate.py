import scanpy as sc
import pandas as pd
import numpy as np
import harmony
import os

# ============================================================================
# USER CONFIGURATION
# ============================================================================
METADATA_FILE = 'metadata/sample_metadata.csv'
# ============================================================================

print("--- [SCRIPT 1/3] STARTING: Data Integration and Annotation ---")

# --- Configuration ---
BASE_DIR = os.getcwd()
DATA_BASE_PATH = os.path.join(BASE_DIR, 'data/')
RESULTS_DIR = os.path.join(BASE_DIR, 'results/processed_data/')
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- 1. Automated Processing of All Samples ---
print("\nStep 1/4: Processing all raw samples...")
df_meta = pd.read_csv(os.path.join(BASE_DIR, METADATA_FILE))
adata_list = []
for index, row in df_meta.iterrows():
    sample_id = row['Sample ID'].strip()
    library_id = row['Library ID'].strip()
    print(f"  - Loading {sample_id}...")
    try:
        adata_sample = sc.read_10x_mtx(f"{DATA_BASE_PATH}{library_id}/filtered_feature_bc_matrix/", var_names='gene_symbols', cache=False)
        adata_sample.var_names_make_unique()
        # Add all metadata columns to the object
        for col in df_meta.columns:
            adata_sample.obs[col] = row[col]
        
        adata_sample.raw = adata_sample # Store raw counts
        sc.pp.filter_cells(adata_sample, min_genes=500)
        adata_sample.var['mt'] = adata_sample.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata_sample, qc_vars=['mt'], inplace=True)
        adata_sample = adata_sample[adata_sample.obs.pct_counts_mt < 15, :].copy()
        sc.pp.normalize_total(adata_sample, target_sum=1e4)
        sc.pp.log1p(adata_sample)
        adata_list.append(adata_sample)
    except Exception as e:
        print(f"    - WARNING: Failed to process {sample_id}. Error: {e}")

# --- 2. Concatenate and Integrate ---
print("\nStep 2/4: Concatenating and integrating with Harmony...")
adata_full = sc.AnnData.concatenate(*adata_list, join='inner', batch_key='sample_id_internal')
sc.pp.highly_variable_genes(adata_full, n_top_genes=3000, batch_key='sample_id_internal')
adata_full = adata_full[:, adata_full.var.highly_variable].copy()
sc.pp.scale(adata_full, max_value=10)
sc.tl.pca(adata_full, svd_solver='arpack', n_comps=50)
harmony.harmonize(adata_full.obsm['X_pca'], adata_full.obs, 'Sample ID', out_key='X_pca_harmony')
sc.pp.neighbors(adata_full, n_pcs=50, use_rep='X_pca_harmony')
sc.tl.umap(adata_full)
sc.tl.leiden(adata_full, resolution=0.8, key_added='leiden_integrated')

# --- 3. Create Complete Raw Counts Layer ---
print("\nStep 3/4: Assembling complete raw count object...")
adata_full.raw = sc.AnnData.concatenate(*[ad.raw.to_adata() for ad in adata_list], join='inner')[':, adata_full.var_names'].copy()

# --- 4. Save Final Object ---
print("\nStep 4/4: Saving final integrated AnnData object...")
final_adata_path = os.path.join(RESULTS_DIR, "final_integrated_with_raw.h5ad")
adata_full.write(final_adata_path)

print(f"--- [SCRIPT 1/3] COMPLETE. Final object saved to {final_adata_path} ---")
