# final_submission

Everything the manuscript needs for submission, named by the number each item
carries in the compiled PDF rather than by its build filename.

    SoLD_manuscript.pdf        31 pages, compiled from main.tex

    figures/                   the six main figures
    supplementary_figures/     the nine supplementary figures
    supplementary_tables/      the ten supplementary workbooks

## Renaming

Three main figures and one supplementary figure had build filenames that did not
match the number they render as. The mapping used here is the one in main.aux,
so every file is named for the number the text actually cites.

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

Every rename was checked by md5, so the bytes are unchanged.

## Supplementary tables

    S1_Table.xlsx   Population structure and ancestry          3 sheets
    S2_Table.xlsx   Genomic heritability                       5 sheets
    S3_Table.xlsx   Species inventory and stability            4 sheets
    S4_Table.xlsx   Class composition and contrasts            6 sheets
    S5_Table.xlsx   Lipid ontology and chemical space          2 sheets
    S6_Table.xlsx   GWAS candidate genes                       4 sheets
    S7_Table.xlsx   GO enrichment                              4 sheets
    S8_Table.xlsx   CTL/LIN candidate overlap                  6 sheets
    S9_Table.xlsx   LINEX reaction mapping                     4 sheets
    S10_Table.xlsx  Frozen lipid species set                   1 sheet

S_Table_index.csv maps every sheet back to the source file it was packed from.

## Before you submit

Six citation keys are used in main.tex but are not in SoLD.bib, so they render
as undefined on pages 14 and 15 of this PDF.

    Matros2017      Zhou2019        Cao2016OsLPR
    Aung2006PHO2    Huang2013PHO2   Park2014NLA

Add them and recompile. Nothing else in the build is broken -- no missing
figures, no undefined labels, two overfull boxes.
