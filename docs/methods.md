# Workflow Methods

`paired_sc` implements a standard paired-tissue workflow:

1. 10x H5 ingestion from a manifest
2. QC metric calculation and filtering
3. normalization, HVG selection, PCA, and Harmony integration
4. neighbors, UMAP, Leiden clustering
5. CellTypist-based annotation or fallback cluster labels
6. donor-aware paired summaries, pseudobulk, and differential abundance
7. standardized figure and report generation

The package is tissue-generic. Annotation remapping is handled through adapters
or explicit label maps rather than domain-specific defaults in the core.


