# Transcriptional dynamics in 4sU-labeled mESCs

This repository contains workflows developed during an MSc research internship to analyze 4sU-labeled single-cell transcriptomic data from mouse embryonic stem cells.

The analyses include quality control, assessment of potential labeling-related quantification bias, transcriptional variability analysis, and inference of transcriptional kinetic parameters.

## Repository structure

### `quality_control_and_4su_bias/`

Quality-control workflow for paired new- and old-RNA count matrices, including:

- cell- and gene-level quality control
- annotation of pluripotent, intermediate, and 2C-like cell states
- assessment of RNA-stability estimates
- checks for potential 4sU-related dropout or quantification bias
- export of processed data for downstream analyses

Parts of the quantification-bias assessment were adapted from approaches implemented in **grandR**.

### `fano_factor_analyses/`

Fano factor-based analyses of transcriptional variability within the pluripotent cell population.

The workflow compares gene-level variability across biological replicates and identifies genes with higher or lower variability than expected from their mean expression.

### `kinetic_inference_workflows/`

Adapted workflows for estimating transcriptional kinetic parameters within the pluripotent cell population.

This directory contains:

- `nascseq2_adapted/`
- `deeptx_adapted/`

These workflows were used for dataset-specific preprocessing, parameter estimation, filtering, visualization, and comparison with external datasets.

## Analysis scope

Annotation of pluripotent, intermediate, and 2C-like cell states was performed as part of the quality-control workflow.

The transcriptional variability and kinetic-inference analyses were restricted to cells classified as pluripotent.

## Adapted methods and attribution

Parts of this repository build on previously published software and analysis workflows:

- **grandR:** adapted for assessment of potential labeling-related quantification bias
- **NASC-seq2:** adapted for transcriptional kinetic inference
- **deeptx:** adapted for kinetic analyses and comparisons with external datasets

Adaptations made for this project include dataset-specific preprocessing, workflow restructuring, parameter selection, quality-control checks, filtering, visualization, and downstream comparisons.

The original publications, repositories, licenses, and details of the modifications should be documented within the corresponding workflow directories.
