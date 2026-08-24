# Derived benchmark data

The manuscript benchmark contains 405 administered vaccine targets from 24 patients across three published human personalized cancer-vaccine studies.

The analysis uses only published, de-identified peptide/mutation-level information reconstructed from public supplementary materials. Third-party raw supplementary files are not redistributed in this repository.

The harmonized analysis table used by `src/reproduce_core_results.py` is named `harmonized_vaccine_targets.csv`. It contains trial, de-identified patient identifier, administered target identifier, peptide/cassette sequence, study-specific immune-response label, identically computed sequence features, and study-native epitope-coverage annotations where reconstructable.

Before the repository is switched from private to public, the derived table will be checked once more for source-license and redistribution compatibility. If redistribution of any derived field is restricted, the repository will instead provide a reconstruction script and provenance manifest while retaining the same analysis interface.
