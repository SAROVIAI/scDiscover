---

# **scDiscover: A Python Pipeline for Single-Cell Biomarker Analysis**

**Version:** 1.0.0
**Contact:** SAROVI Sp. z o.o — [hello@sarovi.pl](mailto:hello@sarovi.pl)

---

## Overview

**scDiscover** is an end-to-end, reproducible bioinformatics pipeline for the analysis of single-cell RNA sequencing (scRNA-seq) data. It processes raw 10x Genomics count data from multiple samples, identifies disease-associated cell populations, and performs deep biological characterization using differential expression, pathway enrichment, and cell-cell communication analysis.

The pipeline is data-agnostic and can be applied to any suitable single-cell dataset.

---

## Features

* **Automated QC & Pre-processing**
  Quality control and filtering for multiple 10x Genomics samples.

* **Robust Data Integration**
  Batch effect correction using the Harmony algorithm to create a unified cellular atlas.

* **Objective Target Cell Identification**
  Quantitative gene signature scoring to identify disease-associated clusters.

* **Deep Biological Characterization**

  * Differential Gene Expression (DGE)
  * Gene Set Enrichment Analysis (GSEA)
  * Optional per-sample pseudobulk DGE

* **Interaction Network Mapping**
  Cell-cell communication inference using **Squidpy** to explore intercellular signaling.

---

## Installation

### 1. Clone the repository

```bash
git clone [URL_to_this_repository]
cd scDiscover_Pipeline
```

### 2. Install Conda

If not already installed, download [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

### 3. Create and activate the environment

```bash
conda env create -f environment.yml
conda activate scdiscover_env
```

---

## Usage

### 1. Prepare Your Data

#### A. Raw Count Data

Place each 10x Genomics `filtered_feature_bc_matrix` folder under the `/data` directory:

```
data/
├── Sample_1/
│   └── filtered_feature_bc_matrix/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
├── Sample_2/
│   └── ...
```

#### B. Metadata File

Place a CSV file in the `/metadata` directory. Use the provided `sample_metadata.csv` as a template. Required columns:

* **Library ID**: Matches folder name under `/data`
* **Sample ID**: Unique identifier
* **Condition**: Biological group (e.g., "Disease", "Control")

---

### 2. Configuration

Edit the following files before execution:

* `scripts/01_process_and_integrate.py`
* `scripts/02_identify_and_characterize.py`

Update key parameters:

* Metadata path
* Condition column
* "Case" vs "Control" labels
* Disease-specific gene signature

---

### 3. Run the Pipeline

#### Step 1: Process and Integrate

```bash
python scripts/01_process_and_integrate.py
```

**Output:**
`results/processed_data/final_integrated_with_raw.h5ad`

---

#### Step 2: Identify and Characterize

```bash
python scripts/02_identify_and_characterize.py
```

**Output:**
Statistical tables and plots in `/results` and `/figures`
Also generates: `results/target_cluster_id.txt`

---

#### Step 3: Communication Analysis

```bash
python scripts/03_run_communication_analysis.py
```

**Output:**
Squidpy plots and communication tables in `/results` and `/figures`

---

## Sample Dataset

A demonstration dataset is available upon request. Please contact **[hello@sarovi.pl](mailto:hello@sarovi.pl)**.

---

## environment.yml

```yaml
name: scdiscover_env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.9
  - scanpy=1.10.3
  - pandas=2.2.0
  - matplotlib
  - seaborn
  - anndata=0.10.7
  - leidenalg
  - igraph
  - pip
  - pip:
      - pydeseq2==0.4.5
      - squidpy==1.2.3
      - harmony-pytorch==0.1.7
      - gseapy==1.0.3
```
