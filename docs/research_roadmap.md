# Research Roadmap

## One-Month Plan

The next month should answer whether the current radiology-plus-microenvironment signal is real enough to justify a deeper research program.

### Workstream 1: Benchmark expansion

- add a second public LUAD cohort with usable clinical or molecular labels
- add at least one mutation-linked benchmark
- add at least one immune-context benchmark

### Workstream 2: Representation stress tests

- compare all-cell, immune-focused, and malignant-focused views
- replace coarse marker compartments with a stronger malignant-state strategy
- test whether gains survive donor-balanced or sample-balanced perturbations

### Workstream 3: Trial-facing evaluation

- score enrichment-style value rather than only cluster agreement
- define simple cohort-slicing criteria from the embedding space
- check whether those slices are more coherent than PCA-derived slices

### Month-One Decision Gate

Continue into a larger program if most of the following are true:

1. Stack beats PCA on at least two trial-adjacent labels.
2. The gains are stable under bootstrap or balancing checks.
3. The strongest signal is biologically interpretable.
4. The story is closer to enrichment or subgroup definition than to generic visualization.

## Half-Year Plan

If the month-one gate is positive, the next half-year should build a multi-cohort translational benchmark stack rather than a single polished case study.

### Workstream 1: Multi-cohort LUAD benchmark battery

- evaluate across multiple public cohorts
- standardize cohort manifests and label mappings
- track performance by cohort, label class, and cell-view strategy

### Workstream 2: Enrichment-oriented methods

- define cohort-slice scores from Stack space
- compare cluster-based versus neighborhood-based subgroup definitions
- evaluate whether Stack-derived slices improve label enrichment over PCA-derived slices

### Workstream 3: Biology interpretation

- annotate cohort slices with pathway, mutation, and microenvironment context
- separate malignant-state signal from sample-composition effects
- identify which use cases favor mixed-cell embeddings versus malignant-focused embeddings

### Workstream 4: External relevance

- test response- or resistance-adjacent public labels when available
- prioritize endpoints that can plausibly inform trial design or translational study design
- treat any clinical-looking story as hypothesis generation unless external validation exists

## Risks To Manage Over Six Months

- the repo may overfit to small, label-rich public cohorts
- trial-adjacent labels may be noisier than expected
- the strongest embedding gains may reflect composition rather than meaningful subgroup biology
- the model-access pathway may change if Stack interfaces or checkpoints shift

## What success would look like

By six months, a successful version of this project should be able to say something specific and defensible:

Stack-derived LUAD cohort slices recover clinically adjacent public phenotypes more consistently than simple transcriptomic baselines, and the gains are strongest in settings where mixed malignant-plus-microenvironment context matters for translational stratification.

If the evidence does not support that claim, the project should narrow to a smaller methodological or benchmarking contribution rather than pretending to support a stronger translational story.
