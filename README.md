# paired-single-cell-pipeline

[![CI](https://github.com/SecondBook5/paired-single-cell-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/SecondBook5/paired-single-cell-pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)

`paired_sc` is the installable package behind the
`paired-single-cell-pipeline` repository. It provides a CLI and Python API for
paired human tissue single-cell RNA-seq analysis, from manifest-driven 10x H5
ingestion through integration, annotation, donor-aware summaries, and
downstream domain analyses.

## Install And Run

Clone the repository and install it into your environment:

```bash
git clone https://github.com/SecondBook5/paired-single-cell-pipeline.git
cd paired-single-cell-pipeline
pip install -e .
```

Alternatively, install directly from GitHub:

```bash
pip install "git+https://github.com/SecondBook5/paired-single-cell-pipeline.git"
```

## Typical CLI Flow

```bash
paired-sc validate \
  --project examples/skin_paired_human/project.yaml \
  --manifest examples/skin_paired_human/manifest.csv

paired-sc run core \
  --project examples/skin_paired_human/project.yaml \
  --manifest examples/skin_paired_human/manifest.csv \
  --workdir ./demo_run

paired-sc run advanced \
  --project examples/skin_paired_human/project.yaml \
  --manifest examples/skin_paired_human/manifest.csv \
  --workdir ./demo_run

paired-sc report \
  --project examples/skin_paired_human/project.yaml \
  --manifest examples/skin_paired_human/manifest.csv \
  --workdir ./demo_run
```

The short `pts` alias is still available, but `paired-sc` is the preferred
public command.

## Workflow Overview

A typical run follows four stages:

1. validate a `project.yaml` and `manifest.csv` pair
2. run the core workflow to load matrices, perform QC, normalize, integrate,
   cluster, and annotate cells
3. run post-core and advanced domains for cell-cycle scoring, differential
   expression, pathway analysis, intercellular communication, trajectory and
   latent modeling, optimal transport, robustness checks, and integration
   quality assessment
4. build a report and export standardized results, figures, and logs

## What The Package Covers

- manifest-driven 10x H5 ingestion
- configurable Scanpy and Harmony preprocessing
- CellTypist annotation with optional remapping adapters
- donor-aware paired summaries and pseudobulk aggregation
- differential abundance analysis
- callable analysis domains:
  - `cell_cycle`
  - `differential_expression`
  - `pathway_enrichment`
  - `liana`
  - `magic`
  - `trajectory`
  - `latent`
  - `regulatory`
  - `cellot`
  - `pseudobulk_de`
  - `pathway_activity`
  - `robustness`
  - `integration_quality`
  - `target_population`
  - `target_subclustering`
- standardized `results/`, `figures/`, `reports/`, and `logs/` outputs

Any registered domain can be run directly with:

```bash
paired-sc domain run \
  --name differential_expression \
  --project examples/skin_paired_human/project.yaml \
  --manifest examples/skin_paired_human/manifest.csv \
  --workdir ./demo_run
```

## Public API

```python
from paired_sc import (
    WorkflowConfig,
    ManifestTable,
    build_report,
    run_core_pipeline,
    run_advanced_domains,
    run_cell_cycle,
    run_differential_expression,
    run_pathway_enrichment,
    run_liana,
    run_magic,
    run_trajectory,
    run_latent,
    run_regulatory,
    run_cellot,
    run_pseudobulk_de,
    run_pathway_activity,
    run_robustness,
    run_integration_quality,
    run_target_population,
    run_target_subclustering,
)
```

Domain entry points are also available directly:

```python
from paired_sc.domains import run_trajectory, run_regulatory, run_cellot, run_target_subclustering
```

## Input Contracts

### `project.yaml`

Defines project metadata, condition labels, batch and donor keys, QC
thresholds, annotation settings, output folder names, and domain
configuration.

### `manifest.csv`

Required columns:

- `sample_id`
- `donor_id`
- `condition`
- `batch`
- `matrix_h5`

Optional columns:

- `sample_name`
- `replicate_group`
- `metrics_csv`
- `metadata_json`

## Standard Outputs

- `results/adata_raw_with_qc.h5ad`
- `results/adata_filtered.h5ad`
- `results/adata_preprocessed.h5ad`
- `results/adata_annotated.h5ad`
- `results/adata_advanced.h5ad`

## Example Profile

The flagship example lives in
[examples/skin_paired_human](examples/skin_paired_human). It includes an
example annotation remapping adapter while leaving the core package itself
tissue-agnostic.
