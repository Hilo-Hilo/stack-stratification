# LUAD Stage Benchmark

## Why this benchmark still matters

The original stage benchmark asked whether unsupervised structure in Stack embeddings could recover a known LUAD pathology label in a small public cohort.

That benchmark is still useful, but its role has changed. It should now be treated as an early-screening benchmark and a cautionary baseline, not as the main success story of the repository.

## Current Result

Dataset:

- cohort: `GSE189357`
- labels: `AIS`, `MIA`, `IAC`
- analysis unit: sample centroids built from single-cell embeddings
- cell budget: `300` cells per sample across `9` samples

All-cell result:

- Stack: `ARI 0.071`, `NMI 0.393`, `1-NN 0.333`, `silhouette 0.041`
- PCA: `ARI -0.118`, `NMI 0.218`, `1-NN 0.111`, `silhouette 0.006`

Bootstrap median:

- Stack: `ARI 0.071`, `NMI 0.393`
- PCA: `ARI -0.118`, `NMI 0.218`

## What this benchmark taught us

Three things are now clear:

1. Stack is better than a simple PCA baseline on this cohort.
2. The absolute stage-concordance signal is still weak.
3. The signal does not improve when restricted to a coarse epithelial-only slice.

That means the stage benchmark does not support a strong "Stack recovers LUAD pathology stage" claim. It does support a more cautious claim:

Stack contains more useful cohort structure than PCA here, but the structure is not a clean pathology-stage axis.

## Compartment Follow-Up

The compartment analysis sharpened the interpretation:

- all-cell: `NMI 0.393`
- immune-only: `NMI 0.393`
- epithelial-only: `NMI 0.144`
- endothelial-only: `NMI 0.493`, but only `8` samples and negative silhouette

So the current signal is at least partly driven by mixed-cell or non-epithelial structure.

## Why this matters for the rest of the repo

The stage benchmark helped establish the program logic:

- keep a hard baseline
- keep a stability check
- do not let qualitative embeddings substitute for quantitative evaluation
- use weak results to redirect the benchmark program toward better labels

That redirect has already happened. Radiological phenotype is now the stronger public benchmark and should remain the lead signal until a better mutation- or immune-linked public task is added.
