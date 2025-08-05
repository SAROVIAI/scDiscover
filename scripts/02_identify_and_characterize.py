
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ============================================================================
# USER CONFIGURATION
# ============================================================================
# Define the column in your metadata that separates your groups
CONDITION_KEY = 'Condition'
# Define the label for your "case" group (e.g., the disease)
CASE_GROUP = 'HGSOC'
# Define the label for your "control" group
CONTROL_GROUP = 'CTRL'
# Define a list of genes that represent a positive signature for your target cells
# Example: For a tumor, this would be tumor markers. For a specific immune cell, use its markers.
TARGET_GENE_SIGNATURE = ['PAX8', 'EPCAM', 'KRT18', 'KRT8', 'WT1', 'MUC16', 'WFDC2']
# Set the minimum percentage of cells in a cluster that must come from the CASE_GROUP
# to be considered a candidate. 95.0 means 95% or more.
ENRICHMENT_THRESHOLD = 95.0
# ============================================================================

print("--- [SCRIPT 2/3] STARTING: Statistical Analysis ---")

# --- Configuration ---
BASE_DIR = os.getcwd()
RESULTS_DIR = os.path.join(BASE_DIR, 'results/')
PROCESSED_DATA_DIR = os.path.join(RESULTS_DIR, 'processed_data/')
TABLES_DIR = os.path.join(RESULTS_DIR, 'statistical_tables/')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures/')
for path in [TABLES_DIR, FIGURES_DIR]: os.makedirs(path, exist_ok=True)
sc.settings.figdir = FIGURES_DIR

# --- 1. Load Data & Calculate Signature Score ---
print("\nStep 1/5: Loading integrated data and calculating signature score...")
adata_full = sc.read_h5ad(os.path.join(PROCESSED_DATA_DIR, "final_integrated_with_raw.h5ad"))
signature_present = [gene for gene in TARGET_GENE_SIGNATURE if gene in adata_full.var_names]
if not signature_present:
    raise ValueError("None of the genes in TARGET_GENE_SIGNATURE were found in the data.")
sc.tl.score_genes(adata_full, gene_list=signature_present, score_name='target_score')
sc.pl.umap(adata_full, color=['leiden_integrated', 'target_score'], show=False, save='_target_score.png')

# --- 2. Objectively Identify the Target Cluster ---
print("\nStep 2/5: Objectively identifying the target cluster...")
cluster_scores = adata_full.obs.groupby('leiden_integrated')['target_score'].mean().sort_values(ascending=False)
score_rank_df = pd.DataFrame(cluster_scores).reset_index()
composition_data = []
for cluster_id in score_rank_df['leiden_integrated']:
    adata_cluster = adata_full[adata_full.obs['leiden_integrated'] == cluster_id]
    class_counts = adata_cluster.obs[CONDITION_KEY].value_counts()
    case_count = class_counts.get(CASE_GROUP, 0)
    total_count = class_counts.sum()
    case_percent = (case_count / total_count) * 100 if total_count > 0 else 0
    composition_data.append({'Cluster ID': cluster_id, 'Mean Target Score': cluster_scores[cluster_id], '% Case': case_percent})
composition_df = pd.DataFrame(composition_data).sort_values(by=['% Case', 'Mean Target Score'], ascending=[False, False])
print("Cluster Ranking (by % Case and Target Score):\n", composition_df)
final_candidates = composition_df[composition_df['% Case'] >= ENRICHMENT_THRESHOLD]

if final_candidates.empty:
    print(f"\nWARNING: No cluster passed the {ENRICHMENT_THRESHOLD}% enrichment threshold. Selecting the top-ranked cluster as the target.")
    TARGET_CLUSTER_ID = composition_df['Cluster ID'].iloc[0]
else:
    TARGET_CLUSTER_ID = final_candidates['Cluster ID'].iloc[0]
print(f"==> OBJECTIVELY IDENTIFIED TARGET CLUSTER: {TARGET_CLUSTER_ID} <==")
with open(os.path.join(RESULTS_DIR, 'target_cluster_id.txt'), 'w') as f:
    f.write(str(TARGET_CLUSTER_ID))

# --- 3. Single-Cell DGE and GSEA on Target Cluster ---
print(f"\nStep 3/5: Running DGE and GSEA on Cluster {TARGET_CLUSTER_ID}...")
adata_target = adata_full[adata_full.obs['leiden_integrated'] == TARGET_CLUSTER_ID].copy()
sc.tl.rank_genes_groups(adata_target, groupby=CONDITION_KEY, groups=[CASE_GROUP], reference=CONTROL_GROUP, method='wilcoxon')
dge_results = sc.get.rank_genes_groups_df(adata_target, group=CASE_GROUP)
dge_results.to_csv(os.path.join(TABLES_DIR, f'dge_cluster_{TARGET_CLUSTER_ID}.csv'))
gsea_prerank_list = dge_results[['names', 'scores']].copy().dropna().sort_values(by='scores', ascending=False)
gsea_results = gp.prerank(rnk=gsea_prerank_list, gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human'], outdir=None)
gsea_results.res2d.to_csv(os.path.join(TABLES_DIR, f'gsea_summary_cluster_{TARGET_CLUSTER_ID}.csv'))
sc.pl.rank_genes_groups(adata_target, n_genes=20, show=False, save=f'_dge_cluster_{TARGET_CLUSTER_ID}.png')

# --- 4. Pseudobulk Analysis on Target Cluster ---
print(f"\nStep 4/5: Attempting pseudobulk analysis on Cluster {TARGET_CLUSTER_ID}...")
try:
    counts_df = pd.DataFrame(adata_target.raw.X, index=adata_target.obs.index, columns=adata_target.raw.var_names)
    counts_df['sample_id_col'] = adata_target.obs['Sample ID'].values
    pseudobulk_counts = counts_df.groupby('sample_id_col').sum()
    genes_to_keep = pseudobulk_counts.sum(axis=0) > 10
    if genes_to_keep.sum() < 1:
        print("  - Pseudobulk DGE skipped: Data is too sparse after filtering.")
    else:
        pseudobulk_counts_filtered = pseudobulk_counts.loc[:, genes_to_keep]
        metadata_col_to_use = [c for c in ['Sample ID', CONDITION_KEY] if c != 'sample_id_col']
        pseudobulk_metadata = adata_target.obs[metadata_col_to_use].drop_duplicates().set_index('Sample ID')
        pseudobulk_metadata = pseudobulk_metadata.loc[pseudobulk_counts_filtered.index]
        dds = DeseqDataSet(counts=pseudobulk_counts_filtered, metadata=pseudobulk_metadata, design_factors=CONDITION_KEY)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=(CONDITION_KEY, CASE_GROUP, CONTROL_GROUP))
        pseudobulk_dge_results = stats.results_df
        pseudobulk_dge_results.to_csv(os.path.join(TABLES_DIR, f'pseudobulk_dge_cluster_{TARGET_CLUSTER_ID}.csv'))
        print("  - Pseudobulk DGE analysis complete.")
except Exception as e:
    print(f"  - Pseudobulk DGE failed. Reason: {e}")

# --- 5. Save Final Object ---
print("\nStep 5/5: Saving final annotated object...")
adata_full.write(os.path.join(PROCESSED_DATA_DIR, "final_annotated_adata.h5ad"))

print(f"--- [SCRIPT 2/3] COMPLETE ---")
