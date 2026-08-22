# SoLD / Chapter 2 analysis pipeline

Every script used for this study, from raw LC-MS QC through the final figures and
tables, renumbered into execution order. Files were copied here from
`scripts/new_script/`, `scripts/chapter2_addons/`, and the `SoLD_paper` repository;
the originals were left in place. Paths inside the scripts were **not** rewritten, so
check each one's input paths before running.

| # | Script | Purpose | Original location |
|---|---|---|---|
| | **Raw QC, annotation curation** | | |
| 01 | `01_raw_QC_check.R` | Raw run-order and QC-sample checks | `SAP/scripts/new_script/01_check_raw_QC.R` |
| 02 | `02_audit_raw_annotations.R` | Audit raw lipid annotations | `SAP/scripts/new_script/05_audit_raw_lipid_annotations.R` |
| 03 | `03_reaggregate_lipid_names.R` | Merge duplicate/renamed features | `SAP/scripts/new_script/06_reaggregate_corrected_lipid_names.R` |
| 04 | `04_audit_duplicate_profiles.R` | Detect duplicate intensity profiles | `SAP/scripts/new_script/07_audit_duplicate_intensity_profiles.R` |
| 05 | `05_assign_nonfocused_classes.R` | Assign non-focused compound classes | `SAP/scripts/new_script/08_assign_nonfocused_compound_classes.R` |
| 06 | `06_remove_nonlipid_annotations.R` | Drop non-lipid annotations | `SAP/scripts/new_script/09_remove_nonlipid_annotations.R` |
| 07 | `07_prepare_SERRF_inputs.R` | Build SERRF upload matrices | `SAP/scripts/new_script/10_prepare_SERRF_inputs.R` |
| 08 | `08_audit_SERRF_traceability.R` | Scan-to-lipid mapping audit (Supp Table S1) | `SAP/scripts/new_script/11_audit_SERRF_feature_traceability.R` |
| 09 | `09_build_preSERRF_inputs.R` | Filtered pre-SERRF inputs | `SAP/scripts/new_script/12_build_preSERRF_filtered_inputs.R` |
| | **Normalization: SERRF and SpATS** | | |
| 10 | `10_SERRF_correction.R` | SERRF drift/batch normalization | `SAP/scripts/new_script/02_SERRF_correction_legacy.R` |
| 11 | `11_postSERRF_annotate.R` | Re-annotate after SERRF | `SAP/scripts/new_script/13_postSERRF_annotate_reaggregate.R` |
| 12 | `12_wrangle_postSERRF.R` | Wrangle post-SERRF matrices | `SAP/scripts/new_script/03_wrangle_postSERRF_data.R` |
| 13 | `13_SpATS_correction.R` | SpATS spatial field correction | `SAP/scripts/new_script/04_SpATS_correction_legacy.R` |
| 14 | `14_run_postSERRF_SpATS.R` | Run SpATS on post-SERRF data | `SAP/scripts/new_script/14_run_postSERRF_SpATS.R` |
| 15 | `15_QC_SERRF_SpATS_diagnostics.R` | QC diagnostics | `SAP/scripts/new_script/15_QC_SERRF_SpATS_diagnostics.R` |
| 16 | `16_figure_SuppS1_QC.R` | FIGURE Supp S1 QC panel | `SoLD/scripts/numbered_pipeline/28_make_updated_SuppFig_S1_QC.R` |
| | **Lipidome description and figures** | | |
| 20 | `20_variance_species_overview.R` | FIGURE Supp S2 top-variance | `SoLD/scripts/numbered_pipeline/03_qc_variance_and_species_overview.R` |
| 21 | `21_high_variance_lipids.R` | High-variance lipid species | `SAP/scripts/new_script/16_high_variance_lipids.R` |
| 22 | `22_lipidome_class_composition.R` | FIGURE 1, Supp S3 S4 S5 S8 | `SoLD/scripts/numbered_pipeline/04_lipidome_class_composition.R` |
| 23 | `23_class_composition_ALR.R` | Class composition, ALR | `SAP/scripts/new_script/17_class_composition_ALR.R` |
| 24 | `24_lipid_species_counts.R` | Species counts per class | `SAP/scripts/new_script/18_lipid_species_counts.R` |
| 25 | `25_lipid_PCA.R` | PCA of lipid matrices | `SAP/scripts/new_script/19_lipid_PCA.R` |
| 26 | `26_noncore_class_CLR.R` | Non-core class CLR composition | `SAP/scripts/new_script/20_noncore_lipid_class_composition_CLR.R` |
| 27 | `27_class_CLR_correlations.R` | Class-level CLR correlations | `SAP/scripts/new_script/21_class_CLR_correlations.R` |
| 28 | `28_class_logratio_stats.R` | FIGURE Supp S6 CLR correlations | `SoLD/scripts/numbered_pipeline/05_class_logratio_oplsda_and_statistics.R` |
| 29 | `29_top_bottom_genotypes.R` | FIGURE class sums top/bottom, upper-ranked | `SoLD/scripts/33_top_bottom_lipid_genotypes.R` |
| 30 | `30_genomic_heritability.R` | Genomic heritability (Supp Table S5) | `SAP/scripts/new_script/22_genomic_heritability_SERRF_individual_class_sums.R` |
| 31 | `31_plot_genomic_heritability.R` | Heritability plots | `SAP/scripts/new_script/23_plot_genomic_heritability.R` |
| 32 | `32_LION_enrichment.R` | LION ontology enrichment (Fig 1D) | `SoLD/scripts/16_LION_enrichment.R` |
| | **Botanical race and population structure** | | |
| 33 | `33_race_structure_tests.py` | FIGURE 6 + race boxplots, Supp Table S6 | `SAP/scripts/chapter2_addons/race_analysis.py` |
| 34 | `34_stats_lite_helpers.py` | Helper stats functions | `SAP/scripts/chapter2_addons/stats_lite.py` |
| | **GWAS** | | |
| 35 | `35_run_GWAS.R` | GWAS driver | `SoLD/scripts/run_section5_gwas.R` |
| 36 | `36_figure_GWAS_stacked.R` | FIGURE 2 and 3 stacked Manhattan | `SoLD/scripts/34_stacked_gwas_condition_plot.R` |
| 37 | `37_gwas_tables_bonferroni.py` | Bonferroni candidate tables | `SoLD/scripts/rebuild_chapter2_gwas_tables_bonf.py` |
| | **LD candidate mapping and CTL/LIN overlap** | | |
| 38 | `38_extract_LD_candidates.py` | LD-based candidate gene extraction | `SAP/scripts/chapter2_addons/extract.py` |
| 39 | `39_gwas_overlap_gene_level.py` | Gene-level CTL/LIN overlap | `SAP/scripts/chapter2_addons/gwas_overlap.py` |
| 40 | `40_gwas_overlap_locus_level.py` | Locus-collapsed overlap | `SAP/scripts/chapter2_addons/locus_overlap.py` |
| 41 | `41_figure_overlap.py` | FIGURE 7 CTL/LIN overlap | `SAP/scripts/chapter2_addons/fig_overlap.py` |
| | **GO enrichment (BP and MF, LD-aware null)** | | |
| 42 | `42_plot_GO_by_trait_layer.R` | GO by trait layer | `SAP/scripts/new_script/24_plot_GO_enrichment_by_trait_layer.R` |
| 43 | `43_GO_BP_sumratio_by_class.py` | GO-BP sum/ratio by class | `SAP/scripts/new_script/25_GO_enrichment_sumratio_by_lipid_class.py` |
| 44 | `44_plot_GO_sumratio.R` | Plot sum/ratio GO | `SAP/scripts/new_script/26_plot_sumratio_GO_by_lipid_class.R` |
| 45 | `45_plot_GO_combined.R` | Combined GO plot | `SAP/scripts/new_script/27_plot_combined_sumratio_GO.R` |
| 46 | `46_GO_locus_collapsed.py` | Locus-collapsed GO | `SAP/scripts/new_script/29_GO_enrichment_locus_collapsed.py` |
| 47 | `47_GO_BP_LDaware_permutation.py` | GO-BP LD-aware null | `SAP/scripts/new_script/29_GO_enrichment_LDaware_permutation.py` |
| 48 | `48_GO_MF_individual_by_class.py` | GO-MF individual lipids | `SAP/scripts/new_script/31_GO_MF_enrichment_individual_by_lipid_class.py` |
| 49 | `49_GO_MF_sumratio_by_class.py` | GO-MF sum/ratio | `SAP/scripts/new_script/32_GO_MF_enrichment_sumratio_by_lipid_class.py` |
| 50 | `50_GO_MF_LDaware_permutation.py` | GO-MF LD-aware null | `SAP/scripts/new_script/33_GO_MF_LDaware_permutation.py` |
| 51 | `51_combine_MF_LDaware_batches.py` | Combine permutation batches | `SAP/scripts/new_script/34_combine_MF_LDaware_batches.py` |
| 52 | `52_plot_GO_LDaware_layers.R` | Per-layer LD-aware GO figures | `SAP/scripts/new_script/35_plot_GO_LDaware.R` |
| 53 | `53_figure_GO_main_and_table.R` | FIGURE 5, Supp GO BP/MF, GO table | `SAP/scripts/new_script/36_GO_figures_and_table.R` |
| | **LINEX reaction network** | | |
| 54 | `54_linex_reaction_network.R` | FIGURE 4 LINEX network | `SoLD/scripts/numbered_pipeline/13_linex_reaction_network.R` |
| 55 | `55_linex_gwas_support.R` | LINEX-GWAS candidate support | `SoLD/scripts/numbered_pipeline/14_linex_GWAS_candidate_support.R` |
| | **Shiny application** | | |
| 56 | `56_shiny_app.R` | SoLD Shiny application (Figure 8 source) | `SAP/app.R` |

## Figure and table provenance

| Output in the manuscript | Produced by |
|---|---|
| Figure 1 lipidomics landscape | `22_lipidome_class_composition.R` (panel D from `32_LION_enrichment.R`) |
| Figure 2, 3 GWAS Manhattan | `36_figure_GWAS_stacked.R` |
| Figure 4 LINEX subnetwork | `54_linex_reaction_network.R` |
| Figure 5 GO-BP LD-aware | `53_figure_GO_main_and_table.R` (writes `Fig_GO_BP_main.png`, renamed) |
| Figure 6 race and structure | `33_race_structure_tests.py` (writes `FigRace_C_variance_explained.png`, renamed) |
| Figure 7 CTL/LIN overlap | `41_figure_overlap.py` (writes `FigOverlap_CTL_vs_LIN.png`, renamed) |
| Figure 8 Shiny workflow | screenshots from `56_shiny_app.R`, assembled by hand |
| Supp S1 QC | `16_figure_SuppS1_QC.R` |
| Supp S2 top-variance | `20_variance_species_overview.R` |
| Supp S3, S4, S5 | `22_lipidome_class_composition.R` |
| Supp S6 CLR correlations | `28_class_logratio_stats.R` |
| Supp non-focused class context | `22_lipidome_class_composition.R` |
| Supp race boxplots CTL/LIN | `33_race_structure_tests.py` |
| Supp class sums top/bottom, upper-ranked | `29_top_bottom_genotypes.R` |
| Supp GO-BP all, GO-MF all, GO table | `53_figure_GO_main_and_table.R` |
| Supp Table S1 SERRF audit | `08_audit_SERRF_traceability.R` |
| Supp Table S5 heritability | `30_genomic_heritability.R` |
| Supp Table S6 race tests | `33_race_structure_tests.py` |
| Supp Tables S21-S24 overlap | `39_gwas_overlap_gene_level.py`, `40_gwas_overlap_locus_level.py` |

## Known gaps

- **`Genes_in_enriched_GO_terms.xlsx` / `.tsv` (Supp Table S19)** has no script in either
  repository. It was generated ad hoc from the GO enrichment output and should be
  rewritten as a script before submission.
- **Figure 8** panels are screenshots of the running Shiny app, not a scripted figure.
- Several scripts share a number in the original folders (two `22_`, two `29_` in
  `new_script/`). Both were kept and renumbered here.
- `35_run_GWAS.R` is the driver; the GEMMA calls themselves were run outside the repo.

