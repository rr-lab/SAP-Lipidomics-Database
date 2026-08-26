# scripts/pipeline/ naming

`<run order>_<output>_<description>.R`

The leading number is the execution order and does not change. The middle token
says what the script produces, so a filename alone tells you where its output
appears in the manuscript.

| Token | Meaning |
|---|---|
| `FigN` | main figure N |
| `SuppFigN` | supplementary figure N |
| `SuppTableN` | supplementary table N |
| none | processing or analysis step that produces no numbered output |

Every script that draws something begins with `source("scripts/pipeline/_common.R")`,
which holds the paths, palettes, `plot_theme` and the save helpers.

## Main figures (Option B numbering)

Old Figure 1, the analytical-comparability panel, was removed; its QC content is
already covered in more detail by the supplementary QC figure.

| Script | Output | Was |
|---|---|---|
| `17_Workflow_schematic.R` | `Figure_Workflow.png` | old Fig 1A. **Placement not assigned** -- graphical abstract or supplementary, your call |
| `18_Fig1_race_structure.R` | `Figure1_Race_Structure.png` | old Fig 2, panels C and D only |
| `23_Fig2_class_composition.R` | `Figure2_Class_Composition.png` | replaces old Fig 3 |
| `19_Fig3_gwas_manhattan.R` | `Figure3_GWAS_Manhattan.png` | old Fig 4 |
| `53_Fig4_go_bp_ldaware.R` | GO-BP LD-aware | old Fig 5. **Not yet split out of `53_figure_GO_main_and_table.R`** |
| `21_Fig5_linex.R` | `Figure5_LINEX.png` | old Fig 6, panels A and C only |
| `56_Fig6_shiny_app.R` | Shiny screenshots | old Fig 7. Assembled by hand |

## Removed

- old Fig 1 (analytical comparability) -- QC lives in the supplementary QC figure
- old Fig 2 panels A, B (top-variance species) -- cross-trial comparison, and judged not to carry biology
- old Fig 3 (CLR contrast, chemical space, LION) -- all three panels were LIN - CTL contrasts
- old Fig 6 panel B (reaction-balance score) -- a LIN - CTL contrast

## Still to do

- `53_Fig4_go_bp_ldaware.R` -- split out of script 53
- `SuppFig*.R` and `SuppTable*.R` -- surviving supplementary list not yet settled
- `main.tex` -- untouched by this restructure
