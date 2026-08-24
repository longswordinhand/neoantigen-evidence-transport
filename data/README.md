# Derived benchmark data

The manuscript benchmark contains 405 administered vaccine targets from 24 patients across three published human personalized cancer-vaccine studies.

The public release is `harmonized_vaccine_targets_features_only.csv.gz.b64`. It is a gzip-compressed CSV encoded as base64 text for repository portability. `src/reproduce_core_results.py` decodes it directly in memory.

To reduce third-party redistribution risk, the release table intentionally omits:

- raw peptide/cassette sequences;
- source-specific target identifiers;
- original patient labels.

Patient labels are replaced by deterministic within-trial integers so grouped validation remains exactly reproducible. The table retains trial label, response variables, sequence-derived numeric features, peptide length, and study-native epitope-density/coverage variables used by the core manuscript analyses.

Third-party raw supplementary files are not redistributed in this repository. Researchers who wish to reconstruct the source-level peptide/HLA data or extend the unified NetMHCpan analysis should obtain the original supplements from the cited publishers and follow the manuscript provenance description.
