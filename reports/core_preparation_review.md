# Independent review of the preparation entry point

Reviewed `scripts/run_preparation.py`, `src/preparation/prepare_source_indices.py`,
and the preparation section of `src/preparation/agp_core.py` on 2026-10-06.

The stated scope is accurate: preparation starts from six supplied processed
tables and does not claim to rerun FASTQ processing or original questionnaire
matching. The archived `FROZEN` directory is used for reference comparisons;
it does not supply the newly generated eligible cohort or fitted parameters.
The preparation recreates the 2877 eligible cohort, 2760 descriptive cohort,
2748 analysis cohort, primary design, CLR, target counts, and PC scores.

The core pipeline also read the regenerated sample identifiers explicitly as
strings and verified exact sequence agreement against `raw_predictors.pkl`.
No current numerical or sample-alignment error was identified.

Two hardening suggestions were sent to the parent task: retain string dtype in
`compare_frames`, and check staged input hashes against the archive when a
staged file already exists, to prevent a stale copied input on a later rerun.
These are repeatability safeguards; no discrepancy was observed in this run.
