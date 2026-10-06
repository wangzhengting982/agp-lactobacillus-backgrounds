# Main-figure reproduction

Run from the project root after preparation, core and deepening analyses finish:

```powershell
.venv/Scripts/python.exe scripts/run_figures.py
```

The entry point refuses to substitute archived model results when a current
statistical output is absent. It reads:

| Panel | Newly calculated source |
|---|---|
| 1a/1b | `runs/cohort/fig1_flow.json`, `fig1_detection.csv` |
| 2a | `runs/core/表S22c_全部15组主成分向量比较.tsv` |
| 2b | `runs/core/表S24i_三背景两目标属效应量.tsv` |
| 2c | `runs/deepening/diet_preprocessing/更广泛饮食调整_全部属间差异.tsv` |
| 3a | `runs/deepening/results/S02_四状态标准化概率.tsv` |
| 3b | `runs/deepening/results/A01_检出后相对丰度.tsv` |

The `source_data` directory beside this README is a **submitted-reference** copy.
Its `original_plotted_values.json` is used only for comparison; it does not supply
the inferential plot values in the normal run. An explicit developer diagnostic
flag can allow frozen Figure 1 input, but that use is reported and was **not used
by the completed normal validation run**.

Outputs go to `runs/figures`: editable SVG, CMYK vector PDF, RGB PNG and CMYK TIFF,
plus rendered previews, source hashes and numerical/export logs. The finalized
main-figure layout is preserved: width 15 cm, black 7 pt labels, 600 dpi raster
exports and embedded color profiles. The script checks 66 plotted rows and 39
interval artists against their input values and compares them with the submitted
reference at relative tolerance 1e-7 and absolute tolerance 1e-9.

The verified rendering environment is Windows with Arial and YouYuan (幼圆).
Fonts are discovered from `%WINDIR%/Fonts`; they are not redistributed here.
Optional `--font-latin` / `--font-cjk` arguments expose font file locations, but
using different typefaces requires a new layout check and is not covered by the
current pixel-match validation. Scientific results do not depend on fonts.

The migrated plotting implementation is derived from the frozen submission
artwork source in `data/archive/11_投稿格式修订_20261006/图件复绘`.
No model is fitted inside the renderer. Rendering success alone is not evidence
of model reproduction; consult the separate statistical execution reports.

Completed local checks are recorded in:

- `reports/figures_provenance.json`
- `reports/figures_visual_QA.json`
- `runs/figures/numeric_plot_audit.json`
- `runs/figures/vector_font_and_raster_audit.json`

On this validated run, all three generated final CMYK PDFs rendered pixel-for-
pixel identically to the archived submitted figure PDFs at 200 dpi. The numerical
differences were floating-point/serialization differences only, maximum about
7.11e-15. Figure provenance is independent of PDF binary equality (PDF creation
timestamps can differ).
