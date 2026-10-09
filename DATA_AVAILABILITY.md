# Data and required inputs

[Project overview](README.md) · [中文说明](README.zh-CN.md)

## Included in this repository

- Analysis and figure-generation code.
- Dependency versions and pipeline entry points.
- Aggregate plotting reference values in `src/figures/source_data/`.
- Local validation records, input-file hashes, and table-to-code mappings.
- Current three-figure renderer, aggregate source values and checks under `src/figures/current/`.

## Current manuscript and access

The manuscript updated on 10 October 2026 contains 30 result tables plus an index and three main figures. Figure 3 combines detection-state probabilities and detected-subset abundance; there is no separate supplementary Figure S1. The N01 table has 42 rows, including the archived at-least-3-read sensitivity results. This consolidation did not refit models.

Processed analysis inputs, fixed classification outputs and archived run records can be requested from the corresponding author, Shi Huang (shihuang@hku.hk). The public repository supplies code and aggregate figure values; it does not supply participant-level records or the complete analysis input archive. Aggregate result deposits do not replace those inputs.

The six packages named `06_分析代码与复现文件_01共06.zip` through `_06共06.zip` are a separately maintained local reproduction archive, provided on request when appropriate. They are not default journal attachments. Merge all six into one empty directory. The current archive has 515 unique files, including the 16 files added on 9 October (the current-figure notice and `12_当前图件复现_20261009/`). The original 499-file statistical baseline and its manifests remain unchanged.

## Inputs for the 6 October statistical validation baseline

The full pipeline expects the following directory structure:

```text
data/
└── archive/                  # 499 unique files in the archived input manifest
reference/
├── 01_论文正文.docx
├── 04_补充材料.docx
└── submitted_results.xlsx
```

| Input | Purpose | Integrity manifest |
|---|---|---|
| `data/archive/` | Source tables, BIOM counts, taxonomic annotations, sample indices, archived reference outputs, and figure reference assets | [`reports/input_manifest.json`](reports/input_manifest.json) |
| `reference/01_论文正文.docx` | Check numerical statements and manuscript tables | [`reference_manifest.json`](reference_manifest.json) |
| `reference/04_补充材料.docx` | Check the accompanying supplementary document | [`reference_manifest.json`](reference_manifest.json) |
| `reference/submitted_results.xlsx` | Compare recalculated results with the 40 result sheets frozen on 6 October 2026 | [`reference_manifest.json`](reference_manifest.json) |

Preserve the relative paths listed in the manifests. The archived input bundle contains **499 unique files**, including the six source files under `02_统计输入/六份源数据/` (five TSV tables and one BIOM file). Providing only these six files is insufficient for the current end-to-end validation workflow, which also checks annotations, archived references, and the three files in `reference/`.

The input bundle and reference documents are not distributed in this code repository. There is currently no complete public input-download endpoint documented here. Request the matching 6 October baseline when running the historical validation workflow; substituting the current consolidated workbook would fail its frozen reference checks.

## Check the inputs

After creating the Python environment, run:

```powershell
.\.venv\Scripts\python.exe scripts/verify_inputs.py
```

The checker verifies SHA-256 hashes for the archived input files and all three reference files. Missing or changed inputs are listed in `reports/input_integrity.json`; the full pipeline stops before fitting if this check fails.

## Output and validation scope

New outputs are generated in `runs/`. Archived tables and the submitted workbook are comparison references; they are not substituted for freshly fitted models in the complete run.

The pipeline reconstructs the primary cohort and model matrices from processed source files. Taxonomic annotation outputs and the longitudinal sample index are supplied inputs. Reproducing the upstream FASTQ processing, initial questionnaire linkage, or sequence classifier requires their separate upstream workflows.

The recorded figure checks use the Windows rendering environment described in [`src/figures/README.md`](src/figures/README.md). Numerical reproducibility and rendered-image agreement are reported separately.

Current figures can be rendered independently from the aggregate values in `src/figures/current/`, using that directory's `requirements.txt`. Its SHA-256 manifest identifies the 14 source, profile and validation files copied from the current archive. These plotting checks do not establish a new end-to-end statistical reproduction run.
