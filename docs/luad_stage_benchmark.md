# LUAD Stage Benchmark

## Question

Do unsupervised strata induced by ARC Institute Stack embeddings recover an already established LUAD pathology label in a small public single-cell cohort?

## First Benchmark Choice

- Cohort: `GSE189357`
- Unit of analysis: single cells embedded with Stack, then aggregated to sample centroids
- Label of interest: `AIS` vs `MIA` vs `IAC`

## Why This Benchmark

- Public and compact
- Explicit pathology labels are available in GEO sample metadata
- Better aligned to an initial exploratory benchmark than inventing private endpoints

## Experimental Shape

1. Download and unpack `GSE189357`.
2. Build a combined AnnData object with sample-level metadata.
3. Downsample to a controlled number of cells per sample for a first-pass CPU run.
4. Generate Stack embeddings from the downloaded checkpoint.
5. Aggregate cell embeddings to sample centroids.
6. Compare unsupervised sample clusters with pathology labels using ARI, NMI, and nearest-neighbor label agreement.
7. Render figures for README inclusion.

## Limits

- This benchmark is intentionally small.
- Sample-level pathology may partly reflect cell-composition shifts rather than purely malignant-state geometry.
- Without cell-type annotations, the first pass should be interpreted as a cohort-level structure check, not a definitive subtype discovery result.
