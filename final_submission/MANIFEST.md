# final_submission

Everything the manuscript needs for submission, named by the number each item
carries in the compiled PDF rather than by its build filename.

    SoLD_manuscript.pdf        compiled from main.tex
    figures/                   the six main figures
    supplementary_figures/     the nine supplementary figures
    supplementary_tables/      the thirty-one supplementary tables

## The PDF is stale

`SoLD_manuscript.pdf` was compiled on 2026-09-21. `main.tex` has changed
substantially since, most recently the supplementary-table renumbering on
2026-09-30. Recompile before submitting. Everything else in this folder is
current as of 2026-09-30.

## Figures

Refreshed 2026-09-30 from `final/fig/`, which is itself byte-identical to the
files `main.tex` compiles against. The copies here had been sitting at the
2026-09-21 build and so predated the LIPID MAPS category rebuild of Figure 2 and
Supplementary Figure S5. Every copy was verified byte-for-byte after renaming.

Three main figures and one supplementary figure have build filenames that do not
match the number they render as. The mapping is the one in main.aux, so every
file is named for the number the text actually cites.

    Figure1.png   <- fig/main/Figure1_Population_Structure.png
    Figure2.png   <- fig/main/Figure2_Class_Composition.png
    Figure3.png   <- fig/main/Figure3_GWAS_Manhattan.png
    Figure4.png   <- fig/main/Figure5_GO_BP_LDaware.png        renamed
    Figure5.png   <- fig/main/Figure5_LINEX.png
    Figure6.png   <- fig/main/Figure7_Shiny_App.png            renamed

    S1_Fig.png    <- fig/supp/SuppFig_S1_Workflow.png
    S2_Fig.png    <- fig/supp/SuppFig_S2_QC_RunOrder_SERRF_PCA_SpATS.png
    S3_Fig.png    <- fig/supp/SuppFig_S3_Class_PCA_Structure.png
    S4_Fig.png    <- fig/supp/SuppFig_S4_Heritability.png
    S5_Fig.png    <- fig/supp/SuppFig_S5_Lipid_Species_Counts.png
    S6_Fig.png    <- fig/supp/SuppFig_S6_CLR_Correlations.png
    S7_Fig.png    <- fig/supp/SuppFig_S7_PCA_Lipids.png
    S8_Fig.png    <- fig/supp/SuppFig_S8_Chemical_Space.png
    S9_Fig.png    <- fig/supp/Figure7_CTL_LIN_Overlap.png      renamed

## Supplementary tables

Thirty-one tables, one per file, renumbered from ten multi-sheet workbooks on
2026-09-30. Each file holds a single sheet, so a citation names a table and
nothing else. Built by `scripts/new_new_script/15_SuppTables_build_flat.R`.

    S1_Table.xlsx   Population structure, group sizes                 12 rows
    S2_Table.xlsx   Population structure, lipid tests                 56 rows
    S3_Table.xlsx   Population structure, PC tests                     8 rows
    S4_Table.xlsx   Heritability, per species structure              430 rows
    S5_Table.xlsx   Heritability, per species paired                 328 rows
    S6_Table.xlsx   Heritability, structure conditioned               26 rows
    S7_Table.xlsx   Heritability, class sums                          34 rows
    S8_Table.xlsx   Heritability, ancestry robustness                 54 rows
    S9_Table.xlsx   Frozen lipid species set                         316 rows
    S10_Table.xlsx  Top variance lipids                               20 rows
    S11_Table.xlsx  Class composition %TIC                            40 rows
    S12_Table.xlsx  Class contrasts CLR and ALR                       25 rows
    S13_Table.xlsx  Class CLR correlation delta                       78 rows
    S14_Table.xlsx  Lipid ratio statistics                            33 rows
    S15_Table.xlsx  Species counts by class and category              25 rows
    S16_Table.xlsx  Composition stability                             26 rows
    S17_Table.xlsx  LION enrichment                                   10 rows
    S18_Table.xlsx  Chemical space                                    26 rows
    S19_Table.xlsx  GWAS candidate genes, CTL individual            1062 rows
    S20_Table.xlsx  GWAS candidate genes, CTL sum ratio               54 rows
    S21_Table.xlsx  GWAS candidate genes, LIN individual            4370 rows
    S22_Table.xlsx  GWAS candidate genes, LIN sum ratio              385 rows
    S23_Table.xlsx  GO enrichment, loci collapsed                    217 rows
    S24_Table.xlsx  GO enrichment, BP all terms                      164 rows
    S25_Table.xlsx  GO enrichment, MF all terms                      180 rows
    S26_Table.xlsx  GO enrichment, genes in enriched terms           919 rows
    S27_Table.xlsx  Overlap by lipid class                            17 rows
    S28_Table.xlsx  Shared candidate genes                           291 rows
    S29_Table.xlsx  LINEX reactions and balance                        4 rows
    S30_Table.xlsx  LINEX GWAS gene support                           15 rows
    S31_Table.xlsx  LINEX branch summary                               5 rows

`S_Table_index.csv` maps each table back to the source CSV or TSV it was built
from, and to the number it carried in the earlier schemes. The full map, with
first-mention line numbers and the reasons four sheets were dropped and four
pairs merged, is `SUPP_TABLE_NUMBERING.md` in the repository root.
