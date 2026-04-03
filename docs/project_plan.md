# Project Plan

## Objective

Build a reproducible LUAD stratification pipeline that uses ARC Institute Stack embeddings as the primary sample representation for unsupervised subtype discovery and downstream translational benchmarking.

## MVP Focus

- Restrict scope to LUAD
- Start with transcriptomics-derived inputs
- Treat Stack embeddings as the representation layer
- Produce cluster assignments and concise cluster profiles
- Join cluster labels back to clinical and genomic metadata
- Support trial-enrichment-adjacent interpretation while avoiding clinical claims

## Pipeline Shape

### 1. Ingestion

- Standardize LUAD sample identifiers
- Track expression matrix location and metadata provenance
- Freeze a cohort manifest for downstream reproducibility

### 2. Embedding Layer

- Define a stable interface for Stack embedding generation or retrieval
- Store embedding provenance, model version, and run metadata
- Keep backend details abstracted behind a single package module

### 3. Stratification

- Start with straightforward unsupervised methods
- Emphasize reproducible configuration over method breadth
- Save cluster assignments, parameters, and quality diagnostics

### 4. Annotation and Evaluation

- Summarize clinical and genomic enrichment per cluster
- Benchmark stability, separability, and metadata coherence
- Document where interpretation is descriptive rather than predictive

## Immediate Next Steps

1. Decide how Stack embeddings will be produced locally or retrieved from a service.
2. Define the expected transcriptomics input schema.
3. Add one small reference cohort manifest and a dry-run test.
4. Implement the first clustering baseline and artifact format.
