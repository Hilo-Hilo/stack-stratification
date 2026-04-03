# Project Plan

## Project Position

This repository has moved past the initial scaffold stage. The first public LUAD benchmark loop is now complete enough to support program decisions.

Current position:

- `GSE189357` is the active public benchmark cohort
- Stack beats a PCA baseline on every all-cell benchmark run completed so far
- direct `AIS/MIA/IAC` recovery is weak and should not be treated as the lead success criterion
- radiological phenotype recovery is materially stronger and more stable
- coarse epithelial-only restriction weakens both the pathology-stage and radiology results

That means the project should now be treated as a translational benchmarking program, not as a generic LUAD clustering repo.

## Working Thesis

A useful thesis for the next phase is:

Foundation-model embeddings may help define LUAD cohort slices that are more clinically adjacent and more microenvironment-aware than what simple transcriptomic baselines recover, but this value needs to be demonstrated through stable benchmark gains on public labels that are closer to trial design than pathology-stage labels alone.

## What Counts As Progress

Progress over the next month should not be measured by adding more plots to the same cohort. It should be measured by whether the repository can support a stronger claim with the same discipline already used here:

- compare Stack against simple baselines
- use cohort- or donor-balanced evaluation where possible
- use bootstrap or perturbation stability checks
- test labels that are plausibly useful for enrichment strategy
- separate descriptive biological structure from anything that sounds clinically predictive

## Immediate Program Priorities

### 1. Public benchmark expansion

- add at least one mutation-linked LUAD benchmark
- add at least one immune-context benchmark
- keep radiology as the current strongest trial-adjacent anchor benchmark

### 2. Better cell-state restriction

- replace coarse marker compartments with a more principled malignant-cell or cell-state restriction
- test whether Stack gains persist after removing obvious compositional shortcuts
- explicitly compare all-cell, malignant-focused, and immune-focused views

### 3. Trial-facing evaluation

- move from generic cluster agreement to enrichment-style comparisons
- ask whether the induced cohort slices would plausibly support trial stratification hypotheses
- quantify whether gains survive bootstrap and sample balancing

### 4. Multi-cohort thinking

- do not let the repo overfit to a single public dataset
- use the current `GSE189357` result as a screen, then look for a second cohort with usable public labels

## Near-Term Deliverables

1. a benchmark battery with stage, radiology, and at least two additional trial-adjacent public labels
2. a cleaner malignant-versus-microenvironment decomposition
3. one GitHub-ready scorecard figure summarizing performance deltas and stability
4. explicit decision gates for whether Stack is earning its complexity over simpler baselines

## Decision Gates

At the end of the next month, the repo should be able to answer these questions:

1. Does Stack show a repeatable gain over PCA on at least two trial-adjacent LUAD labels?
2. Do those gains survive bootstrap or donor-balanced perturbations?
3. Are the strongest gains coming from biologically plausible signal rather than obvious composition artifacts?
4. Is the project converging on an enrichment-style use case, or just on prettier embeddings?

If the answer is mostly no, the project should narrow or be reframed early rather than accumulating more exploratory analyses.
