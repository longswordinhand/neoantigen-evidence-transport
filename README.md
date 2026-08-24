# Neoantigen Evidence Transport

Reproducibility repository for **"A Transportability Audit Framework for Neoantigen Ranking Across Personalized Cancer Vaccine Trials."**

## What this repository tests

Personalized cancer-vaccine studies use heterogeneous sequence, expression, clonality, HLA and presentation evidence to prioritize neoantigens. This repository evaluates a separate question: **does predictive evidence retain the same meaning when an entire clinical trial is unseen?**

The manuscript benchmark contains **405 administered vaccine targets from 24 patients across three published human personalized cancer-vaccine studies**. The primary audit separates four properties that are easily conflated by pooled validation:

1. whole-trial held-out discrimination;
2. within-patient ranking;
3. probability calibration;
4. domain/trial identifiability.

The repository is an evaluation framework and does **not** claim a new state-of-the-art neoantigen predictor.

## Key reproducible findings

Using only identically computed sequence-composition features, leave-one-trial-out AUROC is modest but reproducible across the three held-out domains. Peptide length is a strong protocol marker: length alone identifies source trial far above the three-class chance level. Removing length preserves most cross-trial ranking signal while reducing protocol dependence. Study-native epitope-coverage annotations show inconsistent transport behavior across trials; they are therefore treated as trial-native evidence rather than as a harmonized HLA-binding measurement.

## Repository layout

```text
data/harmonized_vaccine_targets_features_only.csv.gz.b64
                                          compressed features-only release benchmark
src/reproduce_core_results.py             one-command core analysis script
results/key_results.json                  fixed manuscript result snapshot
data/README.md                            benchmark provenance and release boundary
requirements.txt                          tested Python package versions
```

The release benchmark intentionally omits raw peptide/cassette sequences and source-specific target identifiers. Patient labels are replaced by deterministic within-trial integers so patient grouping is preserved without redistributing the original labels. The table retains only the derived variables needed to reproduce the core manuscript analyses.

## Reproduce the core results

Python 3.12 was used for the manuscript analyses.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python src/reproduce_core_results.py
```

The script decodes the compressed feature table in memory and writes `results/reproduced_summary.json`.

## Data provenance

The analysis table was reconstructed from publicly available supplementary materials of three published personalized cancer-vaccine studies. Third-party raw supplementary files and raw peptide sequences are **not redistributed** here. Please obtain original source files from the respective publishers and cite the source studies when reconstructing or extending the benchmark.

Source cohorts used in the benchmark:

- Braun et al., *Nature* (2025), personalized neoantigen vaccination in renal cell carcinoma.
- Wongkongkathep et al., *Clinical Cancer Research* (2026), IMPACT personalized synthetic-long-peptide vaccine trial.
- Palmer et al., *Nature Medicine* (2022), individualized heterologous adenovirus/self-amplifying-RNA vaccine trial.

## Important scope boundaries

- Only administered vaccine targets have post-vaccination outcomes; this is a post-selection benchmark, not a reconstruction of each trial's full candidate universe.
- Response definitions follow the public study-specific assay information and are not assumed to be biologically identical across trials.
- `epitope_coverage` is a study-native annotation derived from each study's published design files. It is **not** a version-matched NetMHCpan/MHCflurry score.
- The endpoint is vaccine-induced immune response, not clinical benefit.
- The GRANITE peptide-level subset contains only three patients and is treated as a small sensitivity domain.

## Reproducibility

Randomized analyses use seed `17`. The manuscript environment used Python 3.12.10, pandas 3.0.3, NumPy 2.4.6, scikit-learn 1.9.0 and Matplotlib 3.11.0.

## Citation

A formal citation will be added when the manuscript receives a DOI. Until then, please cite this repository together with the three source clinical studies listed above.
