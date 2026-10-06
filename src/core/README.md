# Core analysis reproduction

Use the project's independent environment. First run `scripts/run_preparation.py`,
then `scripts/run_core.py`. No original workstation drive is required.

The core runner always performs the following steps in order:

1. Reaggregate the supplied BIOM counts under both supplied full-ASV annotation
   releases, using the original 0.80 support and background exclusion rules;
   rebuild target counts, covariates, CLR matrices, and PC scores.
2. Fit the two-background contrasts, independent stacked-model cross-checks,
   and original 5/10/15-PC and detection-threshold analyses.
3. Trace the annotation coverage and refit the Veillonella sensitivity models.
4. Refit all six PC models for each annotation release. The original code only
   loaded historical fitted arrays; this portable version fits fresh models.
5. Compute the 15 genus-pair PC comparisons for each release.
6. Fit all 936 original background-genus models and compute all 156 omnibus
   and 2340 pairwise background-coefficient tests.
7. Refit the three selected contrasts under read-count and annotation conditions.
8. Compare the new results with all columns of submitted workbook C01–C15 and
   the corresponding archived TSV tables.

Primary models read the newly reconstructed matrices from
`runs/preparation/project/03_分析与结果/association_models`.
Alternative taxonomy models read `runs/core/taxonomy_inputs/R07` or `R06`.
Historical outputs are used only as checks after computation; the saved PC
models used to calculate the taxonomy comparisons are newly fitted here.

This module performs 1036 GLM fits, including the four independent stacked
checks. It limits numerical libraries to one thread. Per-step process IDs,
execution times, and exit codes are recorded in `reports/core_run_status.json`
and retained in `reports/core_run_history/`. Numeric checks and individual
column errors are recorded in `reports/core_verification.json`,
`core_table_comparison.tsv`, and `core_column_comparison.tsv`.

The input BIOM and processed genus table are supplied study inputs. The complete
ASV annotation outputs are also supplied inputs: SINTAX, reference database
construction, raw FASTQ processing, and original diet-questionnaire matching are
not rerun by this module. Reaggregation from the supplied labels is distinct
from rerunning sequence classification.

The `--steps` option is only for targeted diagnosis; validation alone does not
constitute a fresh full run. Use the default full entry point for replication.
