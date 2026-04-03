# Trial Relevance

## What this repository is trying to learn

This project is testing whether foundation-model embeddings can support trial-adjacent cohort enrichment in LUAD.

That is a stricter objective than showing that a model produces clusters or visually separated embeddings. A useful result here should support at least one of the following:

- biologically coherent cohort slicing
- enrichment-style subgroup definition
- biomarker hypothesis generation
- better recovery of clinically adjacent labels than simpler baselines

## What the current results say

The current public benchmark loop supports four claims:

1. Stack beats PCA on the completed LUAD benchmarks.
2. Direct pathology-stage recovery is weak and should not be oversold.
3. Radiological phenotype is the strongest clinically adjacent signal tested so far.
4. The best signal currently lives in all-cell or immune-rich views, not in a coarse epithelial-only restriction.

This is important for trial relevance. The repo is not converging on a direct stage-prediction story. It is converging on a more specific translational hypothesis:

foundation-model embeddings may be useful for cohort slicing that reflects mixed malignant-plus-microenvironment context, especially around clinically adjacent phenotypes.

## What would make the work more trial-facing

The next benchmarks should move closer to actual trial design choices:

1. mutation-linked subgroups such as `EGFR`, `KRAS`, or other public genomic labels
2. immune-context labels that could support immunotherapy enrichment hypotheses
3. response- or resistance-adjacent public endpoints when clean data exist
4. cohort definitions that stay stable under resampling and donor balancing

## What would count as credible evidence

For this repo, a credible trial-facing result should have all of the following:

- a public and reproducible cohort
- a baseline comparison against PCA or another simple transcriptomic representation
- a stability check
- a label or endpoint that is plausibly useful for cohort enrichment
- a biological interpretation that does not depend on overclaiming clinical utility

## Current risks

The current work also exposes the main risks for a longer project:

- sample-composition confounding can masquerade as a useful signal
- single-cohort wins may not transfer
- epithelial-only restrictions can discard relevant structure if the use case is actually microenvironment-aware
- weak labels can waste time even when the engineering stack is sound

The right response is not to abandon the project. The right response is to tighten the benchmark program so that every new result answers a clearer translational question.
