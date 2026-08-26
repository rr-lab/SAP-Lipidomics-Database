# scripts/figures/

One script per figure or table. Run any of them from the **repository root**:

```
Rscript scripts/figures/Fig2.R
```

Every script begins with `source("scripts/figures/_common.R")`, which holds the
paths, palettes, plot theme and save helpers. Change a colour or the theme
there and every figure follows.

Paths can be overridden with environment variables (`SOLD_REPO`, `SOLD_DATA`,
and a per-script `FIGn_OUT`), so the scripts run on a machine where the two
repositories sit elsewhere.

## Main figures

| Script | Output | Status |
|---|---|---|
| `Fig1.R` | `Figure1_Workflow.png` | ready |
| `Fig2.R` | `Figure2_Race_Structure.png` | ready |
| `Fig3.R` | class composition (%TIC) | **not written** -- design not settled |
| `Fig4.R` | `Figure4_GWAS_Manhattan.png` | ready |
| `Fig5.R` | GO-BP LD-aware | **not written** -- must be split out of `53_figure_GO_main_and_table.R` |
| `Fig6.R` | `Figure6_LINEX.png` | ready, panels A and C |
| `Fig7` | Shiny screenshots | no script, assembled by hand |

## Supplementary

Not yet written. The list of surviving supplementary figures has not been
settled, so `SuppFig*.R` and `SuppTable*.R` are pending.

`build_supplementary_workbook.py` already rebuilds the whole workbook from the
CSVs and stays where it is.

## What changed from `scripts/pipeline/`

The pipeline scripts still hold the upstream work (QC, SERRF, SpATS, GWAS,
GO, LD mapping). Only the figure-producing ends moved here. The old
figure scripts remain in `scripts/pipeline/` until this set is complete and
verified, then they should be moved to `scripts/superseded/`.

Two deliberate removals are recorded in the scripts themselves:

- `Fig2.R` drops the old top-variance panels (a cross-trial comparison, and
  judged not to carry biology).
- `Fig6.R` drops the reaction-balance panel (a LIN - CTL contrast).
