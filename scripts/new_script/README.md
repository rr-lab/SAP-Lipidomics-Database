# scripts/new_script/

`<run order>_<output>_<description>`

The leading number is execution order and does not change. The middle token says
what the script produces, so the filename alone tells you where its output lands
in the manuscript.

| Token | Meaning |
|---|---|
| `FigN` | main figure N |
| `SuppFigN` | supplementary figure N |
| `SuppTable` | a supplementary table (number assigned in the workbook) |
| none | processing or analysis step with no numbered output |

## Order

| Script | Produces |
|---|---|
| `01_check_raw_QC.R` | raw run-order and QC checks |
| `02_SERRF_correction_legacy.R` | SERRF drift/batch normalization |
| `03_wrangle_postSERRF_data.R` | post-SERRF matrices |
| `04_SpATS_correction_legacy.R` | SpATS spatial correction, the analysis-ready matrices |
| `05_audit_raw_lipid_annotations.R` | annotation audit |
| `06_reaggregate_corrected_lipid_names.R` | merges duplicate/renamed features |
| `07_audit_duplicate_intensity_profiles.R` | duplicate-profile audit |
| `08_assign_nonfocused_compound_classes.R` | non-focused class assignment |
| `09_remove_nonlipid_annotations.R` | drops non-lipid annotations |
| `10_prepare_SERRF_inputs.R` | SERRF upload matrices |
| `11_audit_SERRF_feature_traceability.R` | scan-to-lipid mapping audit. Its table is **not cited** in the manuscript, so it carries no SuppTable label |
| `12_build_preSERRF_filtered_inputs.R` | filtered pre-SERRF inputs |
| `13_postSERRF_annotate_reaggregate.R` | re-annotation after SERRF |
| `14_run_postSERRF_SpATS.R` | runs SpATS on post-SERRF data |
| `15_SuppFig1_QC_SERRF_SpATS_diagnostics.R` | **Supp Fig 1** QC, SERRF and SpATS diagnostics |
| `16_high_variance_lipids.R` | analysis only; its figures were removed (see below) |
| `20_SuppTable1to3_ratio_species_stats.R` | **Supp Tables S1-S3** ratio statistics, species jackknife, per-class stability |
| `22_SuppTable16_genomic_heritability.R` | **Supp Table S16** genomic heritability (cited at main.tex line 158) |
| `25`, `29`, `31`, `32`, `33`, `34` | GO enrichment computation, no numbered output |
| `36_Fig4_GO_enrichment_and_table.R` | **Fig 4** GO-BP LD-aware, plus the all-terms supplementary figures |

## Figure scripts

The scripts that draw the remaining main figures live in `scripts/pipeline/`
under the same convention, sharing `_common.R` for paths, palettes and
`plot_theme`:

| Script | Output |
|---|---|
| `17_Workflow_schematic.R` | workflow schematic, placement not yet assigned |
| `18_Fig1_race_structure.R` | **Fig 1** botanical race and genetic cluster |
| `23_Fig2_class_composition.R` | **Fig 2** class composition, %TIC |
| `19_Fig3_gwas_manhattan.R` | **Fig 3** GWAS Manhattan |
| `21_Fig5_linex.R` | **Fig 5** LINEX, panels A and C |

`53_Fig4_go_bp_ldaware.R` is not yet split out of `36_Fig4_GO_enrichment_and_table.R`.

## Moved to `_superseded/`

These produce output that is no longer in the manuscript.

- `17_class_composition_ALR.R`, `20_noncore_lipid_class_composition_CLR.R`, `21_class_CLR_correlations.R` -- CLR/ALR cross-trial contrasts
- `18_lipid_species_counts.R`, `19_lipid_PCA.R` -- superseded by the class-composition figure
- `23_plot_genomic_heritability.R` -- no heritability figure in the manuscript
- `24_plot_GO_enrichment_by_trait_layer.R`, `26_plot_sumratio_GO_by_lipid_class.R`, `27_plot_combined_sumratio_GO.R`, `35_plot_GO_LDaware.R` -- superseded by script 36
- `29_GO_enrichment_locus_collapsed.py` -- its `_LOCUS.tsv` output is read by nothing
- `22_spats_reliability_provisional_DO_NOT_USE.R` -- marked do-not-use by its own filename

## Known issue

`16_high_variance_lipids.R` is a 6,000-line monolith that still writes
`Figure1_Lipidomics_Landscape` and `Figure2_OPLS_DA`, neither of which is in the
manuscript any more, and rebuilds supplementary figures under an older numbering.
Its only unique surviving outputs are the ratio/species statistics tables, which
were extracted into `scripts/pipeline/23_supptables_S1_S3.R`. It should not be
run as-is.
