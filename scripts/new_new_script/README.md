# scripts/new_new_script

Every script here produces something that appears in the manuscript: a main
figure, a supplementary figure, or a supplementary table. Nothing else lives in
this folder. Numbering follows the order things appear in the paper. Run each
from the repository root.

| # | Script | Produces | Section |
|---|---|---|---|
| 01 | `01_SuppFig1_workflow.R` | Supp Fig S1 — field-to-database workflow schematic | QC and normalization |
| 02 | `02_SuppFig2_QC_diagnostics.R` | Supp Fig S2 — run-order TIC, QC-RSD, pre/post-SERRF PCA, SpATS residuals | QC and normalization |
| 03 | `03_Fig1_population_structure.R` | **Figure 1** + Supp Table S25 | Botanical race and population structure |
| 04 | `04_SuppFig3_class_pca.R` | Supp Fig S3 + Supp Table S25b | Botanical race and population structure |
| 05 | `05_SuppTableS5G_composition_stability.R` | Supp Table S5G — leave-one-race-out and 80% subsample | Botanical race and population structure |
| 06 | `06_Fig2_class_composition.R` | **Figure 2** + Supp Tables S5D, S5E | Lipidome overview |
| 07 | `07_SuppFig4_chemical_space.R` | Supp Fig S4 + Supp Table S5F | Lipidome overview |
| 08 | `08_SuppFig7_class_correlations.R` | Supp Fig S7 — within-trial CLR correlations | Lipidome overview |
| 09 | `09_Fig3_gwas_manhattan.R` | **Figure 3** — recurring candidate loci | GWAS |
| 10 | `10_SuppTableS16_GO_collapse_and_rank.R` | Supp Table S16 + `tables/chapter2_go_enrichment.tex` | GO enrichment |
| 11 | `11_Fig4_GO_enrichment.R` | **Figure 4** — GO terms by lipid class | GO enrichment |
| 12 | `12_Fig5_linex.R` | **Figure 5** — LINEX subnetwork and GWAS support | LINEX |
| 13 | `13_SuppFig8_ctl_lin_overlap.R` | Supp Fig S8 — CTL/LIN candidate overlap | GWAS overlap |
|  | `_common.R` | Paths, palettes, `plot_theme`, shared data builders. Sourced by all of the above. | |

Figure 6 is a Shiny application screenshot
(`fig/main/Figure5_shiny_triterpenoid_AB.png`) and has no generating script.

Supplementary figures S5 (species PCA) and S6 (species counts) also have no
script in this folder; they were carried over unchanged from the earlier
pipeline.

## Ordering

`10` must run before `11`, which reads its locus table. Everything else is
independent, though all of them read data produced by the preprocessing
pipeline.

## Preprocessing is not in this folder

The scripts that built the data these consume — SERRF normalisation, post-SERRF
annotation and reaggregation, SpATS spatial correction — are in
`_preprocessing_moved/`. Identical copies remain in `scripts/new_script/`.

Two of them are live dependencies rather than history.
`13_postSERRF_annotate_reaggregate.R` writes
`postSERRF_{CTL,LIN}_named_reaggregated.csv`, and `14_run_postSERRF_SpATS.R`
writes `7_CTL_lipid_SpATS_fitted.csv` and `8_LIN_lipid_SpATS_fitted.csv`. Panels
G and H of Supp Fig S2 read those four files directly. Do not delete
`scripts/new_script/` without moving the preprocessing scripts somewhere
permanent first.
