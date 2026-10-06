# AGP Lactobacillus Backgrounds

**Analysis code for genus-specific gut microbial backgrounds in the American Gut Project.**

[中文说明](README.zh-CN.md) · [Data and inputs](DATA_AVAILABILITY.md) · [Validation results](reports/final_validation.json) · [Table-to-code map](reports/table_map.csv)

## Overview

This project examines how gut microbial associations differ across six lactic acid bacterial genera, with a focus on *Lactobacillus* and *Lacticaseibacillus*. It includes cohort preparation, background-genus models, dietary adjustment, sensitivity analyses, and reproduction of the manuscript figures.

Companion manuscript: **乳杆菌属与乳酪杆菌属的肠道菌群背景差异**.

| Study component | Scope |
|---|---|
| Primary analysis | 2,748 adults |
| Microbial backgrounds | 156 background genera |
| Target genera | 6 |
| Verified manuscript outputs | 40 result tables and 3 main figures |
| Validated environment | Python 3.12 on Windows |

## Repository structure

| Directory | Contents |
|---|---|
| [`src/preparation/`](src/preparation/) | Cohort selection and model-input construction |
| [`src/core/`](src/core/) | Genus contrasts and taxonomy sensitivity analyses |
| [`src/deepening/`](src/deepening/) | Dietary PCA, adjusted models, detection states, and abundance analyses |
| [`src/sensitivity/`](src/sensitivity/) | Additional robustness analyses |
| [`src/figures/`](src/figures/) | Figure generation and aggregate reference plot values |
| [`scripts/`](scripts/) | Pipeline runners and validation tools |
| [`reports/`](reports/) | Recorded local validation results and input manifests |
| [`docs/`](docs/) | Detailed validation notes |

## Run the analysis

The pipeline uses separately supplied study inputs. Place the complete input archive in `data/archive/` and the three comparison files in `reference/`, following the [input guide](DATA_AVAILABILITY.md).

With Python 3.12 installed, run these commands from the repository root in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.\.venv\Scripts\python.exe scripts/run_all.py
```

The pipeline checks input hashes, reconstructs the cohort and model matrices, fits the models, verifies the result tables, and redraws the figures. Generated outputs are written to `runs/`; execution and comparison records are written to `reports/`.

Package versions are pinned in [`requirements-lock.txt`](requirements-lock.txt). The validated figure-rendering environment uses Arial and YouYuan (幼圆) on Windows; font and rendering details are documented in [`src/figures/README.md`](src/figures/README.md).

## Validation

The recorded local run verified all **40 manuscript result tables**. The largest absolute numerical difference was **8.24 × 10⁻¹³**. All **3 main figures** were regenerated from newly calculated model outputs and matched the submitted PDFs when rendered at 200 dpi.

- [Full pipeline execution](reports/full_run_status.json)
- [Result and figure verification](reports/final_validation.json)
- [Table-to-code mapping](reports/table_map.csv)
- [Detailed validation notes (Chinese)](docs/本地复现验证报告.md)

These records describe the local validation run on 6 October 2026. Workstation-specific interpreter paths in selected reports have been replaced with portable placeholders; numerical results, versions, timestamps, and exit codes are preserved. Repository-status fields in older reports describe their recording date.

## Data availability and scope

This repository contains code, aggregate plotting references, and validation records. Participant-level tables, the complete input archive, manuscript documents, and the validation workbook are maintained separately. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for the required files and integrity checks.

Reproduction starts from processed tables, BIOM counts, supplied taxonomic annotations, and a fixed longitudinal sample index. It covers cohort reconstruction, statistical models, tables, and figures. Raw-read processing, initial questionnaire matching, and sequence classification precede this pipeline. The study evaluates cross-sectional associations.
