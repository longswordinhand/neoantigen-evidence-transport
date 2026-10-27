# Neoantigen Evidence Transport

Reproducibility repository for **"Benchmarking Neoantigen Prioritization Across Personalized Cancer Vaccine Evidence Domains Reveals Endpoint-Dependent Method Choice"** (CMPB submission, 2026-09-18).

This repository is a benchmark, not a new predictor. It compares five existing HLA-I scoring views across four personalized cancer-vaccine evidence domains, on two endpoints:

- pooled target-level AUROC
- equally weighted macro within-patient AUROC

## What is in this version

- **2,495 exact-common targets across 366 domain-patient instances** (428 responsive)
- Domains: GRANITE 2022 (57 targets / 3 patients / 40 responsive), Ott melanoma (84 / 6 / 14), RCC 2025 (128 / 9 / 74), Zhao 2026 (2,226 / 348 / 300)
- Five scoring views: MHCflurry 2.2 affinity, MHCflurry 2.2 presentation, MixMHCpred 3.0, PRIME 2.1, BigMHC-EL
- Main finding (as reported): no method ordering is stable across domains; in Zhao, the pooled-AUROC point leader and the macro within-patient AUROC point leader differ, and no pairwise 95% interval excludes zero. Do not read a point lead as statistical superiority.

## Repository layout

```text
data/benchmark_features_only.csv    derived, features-only benchmark table (2,495 rows)
scripts/reproduce_primary.py        one-command primary point-metric reproduction
results/reference_metrics.csv       frozen primary point metrics (includes BigMHC-IM rows for supplement transparency; the primary panel is the five views above)
results/reference_rank_probabilities.csv  frozen 5,000-resample method-leader frequencies
results/figure_source_data.csv      exact numeric source data for primary Figures 1 and 2
ENVIRONMENT.txt                     software versions for the core analysis and BigMHC reproduction environment
requirements.txt                    tested Python package versions for the core analysis
```

## Data boundary

`data/benchmark_features_only.csv` contains only: evidence-domain label, within-domain anonymized patient group, within-domain target index, binary response label, and the five primary scalar scoring views. It intentionally omits raw peptide/cassette sequences, HLA strings, source target identifiers, original source patient identifiers, and third-party supplementary files.

## Reproduce

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/reproduce_primary.py
```

The script writes `results/reproduced_metrics.csv` and checks it against `results/reference_metrics.csv` for the five primary views. The frozen check in the submission bundle was: max absolute metric difference about 8.3e-17, PASS.

## Important scope boundaries

- The bootstrap / rank-probability outputs are supplied as frozen reference results. Their deterministic generation scripts and seed family (20260831) are preserved in the full project archive and are not in this compact repository.
- BigMHC-IM is not in the primary panel. Its rows in `reference_metrics.csv` are supplement transparency only; the reproduction script checks the five primary views.
- Response definitions follow the public study-specific assay information and are not assumed to be biologically identical across trials.
- The endpoint is vaccine-induced immune response, not clinical benefit.
- GRANITE (3 patients), Ott (6 patients) and RCC (9 patients) are small domains; their orderings are descriptive only.

## Previous version

The earlier 405-target / 24-patient / 3-study version of this repository is preserved, unchanged, in branch `archive/jbi-405-freeze`. Do not mix its numbers with this version.

## Data provenance

The analysis table was reconstructed from publicly available supplementary materials of published personalized cancer-vaccine studies. Third-party raw supplementary files and raw peptide sequences are not redistributed here. Please obtain original source files from the respective publishers and cite the source studies when reconstructing or extending the benchmark.
