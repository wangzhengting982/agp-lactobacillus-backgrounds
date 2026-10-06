# Dietary and extension analyses

Run from the project root with the project's isolated Python environment:

```text
.venv/Scripts/python.exe scripts/run_deepening.py
.venv/Scripts/python.exe src/deepening/validate_deepening.py
```

`run_deepening.py` sequentially executes 13 stages. It fixes numerical libraries
to one thread, writes a separate execution log per stage, and stops on errors.
`--from-stage` and `--through-stage` support resuming an interrupted execution.

## What is regenerated

1. Primary-cohort food and nutrient records are rejoined from the six input
   tables; no archived pickle is loaded as the starting dietary data.
2. Food and nutrient preprocessing, median imputation, log transforms, clipping,
   scaling and PCA are repeated using the original specification. The new
   scores are passed to all downstream extended models.
3. Dietary-adjusted selected associations; joint background/family models;
   four-state multinomial models and standardized probabilities; nonzero
   abundance models; longitudinal conditional logistic models; ASV-level
   models; robustness and sex-interaction models; high-read influence checks.
4. Sex-related additional covariates are rejoined from sample metadata.

Outputs go to `runs/deepening/diet_preprocessing`, `runs/deepening/results`
and `runs/deepening/sex_results`. The validation script reads the submitted
Excel solely for comparison; fitted results never come from that workbook.

## Inputs and scope

The baseline design, standardized background CLR and target-count matrices in
`runs/preparation/project/03_分析与结果/association_models` are inputs to this
module. Run upstream preparation first; these inputs are regenerated from the
input tables. The all-sample target counts likewise come from that fresh
preparation directory. The R06/R07 taxonomy labels and longitudinal sample
index are frozen inputs; running this module does not
re-run taxonomic classification or claim reprocessing from sequencing reads.

`AGP_ASSOC_ROOT` can select an alternative directory containing the same design
and count files; `AGP_PREPARATION_ROOT` can select the preparation project;
`AGP_DEEPENING_OUTPUT` can select a separate output directory.
All defaults are relative to the project directory, not the working directory.

Original R01 numerical failures in male/no-yogurt strata are retained, along
with the specified complete-BMI R04/R05 follow-up analyses. Model warnings are
not suppressed or silently converted to successful results.

`port_sources.py` is the auditable source adapter. It makes exact path-only
substitutions in copied original scripts and records every replacement and
source hash in `reports/deepening_portability_changes.json`. It never changes
the files in `data/archive`.
