# Refactored QC Workflow

This repository contains the refactored QC workflow developed during my MSc internship on 4sU-labelled/scNT-seq data from mouse embryonic stem cells (mESCs) (from Qiu et al., 2020).

This workflow serves as the bridge from published processed scNT-seq paired new (`C`) and old (`T`) RNA gene-by-cell count matrices toward plots for biological interpretation. Starting from paired new and old RNA gene-by-cell count matrices, it performs QC (Scanpy), cell-state annotation, check for 4sU amplification bias (grandR based), and exports for downstream variability and burst-kinetics analyses.

Steps of this workflow include:
- constructing an `AnnData` object from paired new and old RNA count matrices
     These matrices are used to build a layered `AnnData` object containing:
        `C`
        `T`
        `total`
        `ntr`
- computing quality-control metrics...
        - genes detected per cell
        - total UMI counts
        - mitochondrial fraction
        - mean NTR per cell
- and filtering low-quality cells and low-informative genes
        - filter low-quality cells
        - filter low-information genes
- performing normalization / log transformation / HVG selection / neighborhood graph / PCA / UMAP / Leiden-based structure analysis
- annotating cells into `Pluripotent`, `Intermediate`, and `2-cell like` states using marker-based scores
- estimating gene-level RNA stability, including half-life, global degradation rate, and global synthesis rate
- measuring 4sU dropout per cell state
- Export processed outputs (state-annotated `.h5ad` files and matrices for burst-kinetics analyses)
- comparing inferred bursting parameters (from separate NASC-Seq2 script) to external reference datasets
- exporting quality-controlled single-cell objects, state-annotated cell populations and processed outputs for downstream analyses

## Repository structure

### `run_qc_report.py`
Main entry point for the workflow.  
Runs the full pipeline from loading data to export of figures and processed outputs.

### `config.py`
Central configuration file.  
Defines paths, filenames, marker sets, plotting settings, state labels, and analysis constants.

### `io_utils.py`
Input/output helper functions.  
Loads required files and saves figures, tables, and `.h5ad` outputs with consistent naming.

### `plotting.py`
Plotting helper functions.  
Contains reusable routines for QC figures, UMAP visualizations, validation plots, and diagnostic plots.

### `pipeline.py`
Core analysis logic.  
Contains the main computational steps for QC, dimensionality reduction, annotation, RNA stability / turnover estimation, validation, and export.

## Example figure: cell-state annotation

<p align="center">
  <img src="examples/umap_cell_states.svg" width="300">
</p>

This UMAP shows the broad state annotation used in the workflow, separating cells into `Pluripotent`, `Intermediate`, and `2-cell like` populations. These state labels are used later to structure downstream analyses and interpret transcriptional heterogeneity in the mESC population.

## Input files

The workflow expects a `data/` folder at the repository root containing the required input files.

Expected filenames:

- `mESC-WT-rep1_C.txt`
- `mESC-WT-rep1_T.txt`
- `41592_2017_BFnmeth4435_MOESM4_ESM.xls`
- `scNTseq_params.xlsx`
- `GSM4671630_CK-TFEA-run1n2_ds3_gene_exonic.intronic_tagged.dge.txt`

These files are not included in the public repository.

### Provenance

The main mESC input matrices used in this workflow were taken from processed supplementary files associated with the scNT-seq study by Qiu et al. This means the workflow does **not** start from raw FASTQ files, but from published gene-by-cell count matrices that were already generated upstream by the original study.

These processed input matrices still required substantial downstream analysis, including:

- construction of the `AnnData` object
- QC metric calculation
- cell and gene filtering
- dimensionality reduction and clustering
- cell-state annotation
- RNA stability / turnover estimation
- comparison to external reference datasets

Additional local reference files used in the workflow were derived from published studies, including:

- scNT-seq reference material from Qiu et al.
- SLAM-seq supplementary material used for half-life comparison
- processed external reference data used in stability/dropout diagnostics

## Main outputs

The workflow writes results to the configured results directory, for example:

`results/_rep1_fix/`

Outputs include:

- QC violin and scatter plots
- PCA variance plot
- UMAP visualizations
- cell-state annotation outputs
- RNA stability estimates, including half-life, global deg/syn rates
- state-annotated `.h5ad` objects
- per-state matrix exports
- HVG-based exports for downstream analyses

## How to run

From the repository root:

```bash
python refactored_qc_workflow/run_qc_report.py
