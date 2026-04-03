# Trial Relevance

## Why this repo exists

This project is meant to test whether foundation-model embeddings can help define trial-relevant cohort structure in LUAD.

That is a higher bar than showing a visually interesting embedding or an unsupervised cluster map. A useful result here should say something about whether the representation can support:

- cohort enrichment
- biomarker hypothesis generation
- translational subgroup definition
- comparison against simpler baselines already available to clinical research teams

## What the current result means

The current `AIS/MIA/IAC` benchmark is useful because it is public, compact, and easy to reproduce. But it is only an initial screening benchmark.

The current evidence says:

- Stack is stronger than a simple PCA baseline on this cohort
- the signal is still weak for direct pathology-stage recovery
- the signal is not rescued by a simple epithelial-only restriction
- some of the structure appears to track non-epithelial or microenvironmental variation

## How to keep the work trial-facing

Future iterations should ask questions like:

1. Does Stack define a cohort slice that could plausibly enrich for an immune context relevant to immunotherapy trials?
2. Does Stack improve separation of mutation-linked biology over simpler expression baselines?
3. Are the identified strata stable enough to support prospective cohort definitions?
4. Is any enrichment signal retained after donor balancing and bootstrap checks?

## Recommended benchmark order

1. Public label sets tied to immune context or mutation status
2. Public labels tied to treatment response or resistance, when available
3. Enrichment-style evaluation rather than pure clustering quality
4. Stability and baseline comparisons on every benchmark

## What to avoid

- Claiming clinical utility from weak unsupervised alignment
- Treating visually separated embeddings as evidence of actionable enrichment
- Using trial language without baseline comparisons and stability checks
