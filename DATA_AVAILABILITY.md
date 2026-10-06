# Data and required inputs

[Project overview](README.md) · [中文说明](README.zh-CN.md)

## Included in this repository

- Analysis and figure-generation code.
- Dependency versions and pipeline entry points.
- Aggregate plotting reference values in `src/figures/source_data/`.
- Local validation records, input-file hashes, and table-to-code mappings.

## Inputs supplied separately

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
| `reference/submitted_results.xlsx` | Compare recalculated results with the 40 submitted result sheets | [`reference_manifest.json`](reference_manifest.json) |

Preserve the relative paths listed in the manifests. The archived input bundle contains **499 unique files**, including the six source files under `02_统计输入/六份源数据/` (five TSV tables and one BIOM file). Providing only these six files is insufficient for the current end-to-end validation workflow, which also checks annotations, archived references, and the three files in `reference/`.

The input bundle and reference documents are not distributed in this code repository. There is currently no complete input-download endpoint documented here. Aggregate result deposits are a separate resource and do not replace these pipeline inputs.

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
