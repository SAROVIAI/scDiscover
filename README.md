# scDiscover: A Python Pipeline for Single-Cell Biomarker Analysis

**Version:** 1.0.0
**Contact:** SAROVI Sp. z o.o - hello@sarovi.pl

## Overview

**scDiscover** is an end-to-end bioinformatics pipeline for the analysis of single-cell RNA sequencing data. It is designed to process raw 10x Genomics count data from multiple samples, identify cell populations of interest associated with a specific biological condition (e.g., disease vs. healthy), and perform deep characterization using pathway and cell-cell communication analysis.

This pipeline is built for reproducibility and is designed to be data-agnostic, allowing researchers to apply it to their own single-cell datasets.

## Features

-   **Automated QC & Pre-processing:** Handles quality control and filtering for multiple 10x Genomics samples.
-   **Robust Data Integration:** Uses the Harmony algorithm to correct for batch effects between samples, creating a unified cellular atlas.
-   **Objective Target Cell Identification:** Employs a quantitative gene signature scoring method to objectively identify cell clusters associated with a user-defined condition.
-   **Deep Biological Characterization:**
    -   Performs **Differential Gene Expression (DGE)** analysis.
    -   Runs **Gene Set Enrichment Analysis (GSEA)** to uncover active biological pathways.
    -   Optionally attempts a per-sample **Pseudobulk DGE** for statistical validation.
-   **Interaction Network Mapping:** Infers the cell-cell communication network using Squidpy to understand intercellular signaling.

## Installation

1.  **Clone this repository:**
    ```bash
    git clone [URL_to_this_repository]
    cd scDiscover_Pipeline
    ```
2.  **Install Conda:** If you don't have it, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
3.  **Create and activate the Conda environment:**
    ```bash
    conda env create -f environment.yml
    conda activate scdiscover_env
    ```
    This command creates an environment named `scdiscover_env` and installs all required packages.

## Usage

### 1. Data Input

The user must provide their own data in two parts:

**A. Raw Count Data:**
Place your 10x Genomics `filtered_feature_bc_matrix` directories inside the `/data` folder. The structure should be:

data/
├── [Sample_1_Folder_Name]/
│ └── filtered_feature_bc_matrix/
│ ├── barcodes.tsv.gz
│ ├── features.tsv.gz
│ └── matrix.mtx.gz
├── [Sample_2_Folder_Name]/
│ └── ...
└── ...


**B. Metadata File:**
Create a CSV file and place it in the `/metadata` directory. An example `sample_metadata.csv` is provided. The file **must** have the following columns:
-   `Library ID`: The exact name of the sample folder in the `/data` directory.
-   `Sample ID`: A unique identifier for the sample/patient.
-   `Condition`: A column describing the biological group (e.g., 'Disease', 'Control', 'Treatment', 'Pre-treatment').

### 2. Configuration

Before running, you must configure the key parameters at the top of the `scripts/01_process_and_integrate.py` and `scripts/02_identify_and_characterize.py` files. This includes setting the metadata file path, the condition column name, the "case" vs "control" labels, and the gene signature for your disease of interest.

### 3. Running the Pipeline

Execute the scripts from the main `scDiscover_Pipeline/` directory in order.

**Script 1: Process and Integrate**
This script reads all raw data, performs QC, integrates the samples, performs clustering, and saves a final, processed AnnData object.

```bash
python scripts/01_process_and_integrate.py

Output: results/processed_data/final_integrated_with_raw.h5ad

