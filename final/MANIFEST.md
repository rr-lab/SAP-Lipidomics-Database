# final/ -- the figures and tables the manuscript compiles against

Built 2026-09-17. Everything `main.tex` includes or cites is in here, and nothing
else. Use this folder, not `fig/` or `table/supp/` at the repository root --
those keep intermediate outputs that are not all current.

    final/fig/main/          6 main figures
    final/fig/supp/          9 supplementary figures
    final/table/workbooks/  10 numbered workbooks plus the index

26 files. There is no `final/table/supp/`. The manuscript cites Supplementary
Tables S1 to S10 and names a sheet inside them -- it never cites a source file --
so only the workbooks belong here. The thirty-nine CSV and TSV files they are
packed from live at `table/supp/` in the repository root and still carry the
pre-consolidation numbers (S5A, S6c, S24b, S31), which is exactly the second
numbering scheme this folder exists to keep out of your way.

All 28 distinct Supplementary Table citations in main.tex were checked against
the workbooks on 2026-09-17. Every one resolves to a workbook and, where the text
names a sheet, to that sheet.

## Figures rebuilt in this pass

Each of these was regenerated from the current data and the current lipid class
file.

- `fig/main/Figure1_Population_Structure.png`  <-  `scripts/new_new_script/03_Fig1_population_structure.R`
- `fig/main/Figure2_Class_Composition.png`  <-  `scripts/new_new_script/06_Fig2_class_composition.R`
- `fig/main/Figure5_GO_BP_LDaware.png`  <-  `scripts/new_new_script/11_Fig4_GO_enrichment.R`
- `fig/supp/Figure7_CTL_LIN_Overlap.png`  <-  `scripts/new_new_script/13_SuppFig8_ctl_lin_overlap.R`
- `fig/supp/SuppFig_S1_Workflow.png`  <-  `scripts/new_new_script/01_SuppFig1_workflow.R`
- `fig/supp/SuppFig_S2_QC_RunOrder_SERRF_PCA_SpATS.png`  <-  `scripts/new_new_script/02_SuppFig2_QC_diagnostics.R`
- `fig/supp/SuppFig_S3_Class_PCA_Structure.png`  <-  `scripts/new_new_script/04_SuppFig3_class_pca.R`
- `fig/supp/SuppFig_S4_Heritability.png`  <-  `scripts/new_script/22g_Fig_heritability.R`
- `fig/supp/SuppFig_S5_Lipid_Species_Counts.png`  <-  `scripts/new_new_script/06c_SuppFig6_species_counts.R`
- `fig/supp/SuppFig_S6_CLR_Correlations.png`  <-  `scripts/new_new_script/08_SuppFig7_class_correlations.R`
- `fig/supp/SuppFig_S8_Chemical_Space.png`  <-  `scripts/new_new_script/07_SuppFig4_chemical_space.R`

## Figures carried over unchanged

These are byte-identical to what was already on disk. Two of them could not be
rebuilt here; two have no generating script at all.

- `fig/main/Figure3_GWAS_Manhattan.png`  <-  `scripts/new_new_script/09_Fig3_gwas_manhattan.R`
  Reads the per-trait GWAS result files under GWAS_result/CTL_ind/, about 500 MB each. Too large to move into the cloud session, and this machine's local shell was unavailable. Nothing about it is stale -- it does not depend on the lipid class file, the species set or the candidate master.

- `fig/main/Figure5_LINEX.png`  <-  `scripts/new_new_script/12_Fig5_linex.R`
  Same reason. Panel B reads LIN_ind GWAS result files directly.

- `fig/main/Figure7_Shiny_App.png`  <-  no script
  A screenshot of the Shiny application. No generating script exists and none is needed.

- `fig/supp/SuppFig_S7_PCA_Lipids.png`  <-  no script
  No generating script exists anywhere in the repository. Carried over from the earlier pipeline, as scripts/new_new_script/README.md records.

To rebuild the two GWAS figures yourself, from the repository root:

    Rscript scripts/new_new_script/09_Fig3_gwas_manhattan.R
    Rscript scripts/new_new_script/12_Fig5_linex.R

Both take `GWAS_ROOT` from the environment and default to
`/Users/nirwantandukar/Documents/Research/data/SAP/GWAS_result`, so on your own
machine they need no arguments. Copy the two PNGs into `final/fig/main/`
afterwards.

## Tables

The ten workbooks are packed from the 39 source files by
`scripts/new_new_script/15_SuppTables_build_workbooks.R`, which asserts that
every source file lands in exactly one sheet.

Four were stale and were rebuilt in this pass. They dated from 24-25 August and
so predated both the 3 September species deduplication and the 16 September
lyso-name normalisation.

| Table | Generator | What moved |
|---|---|---|
| S1 ratio statistics | `scripts/new_script/20_SuppTable1to3_ratio_species_stats.R` | 16 of 33 ratios changed magnitude. No direction flipped. LPC/LPE goes 0.38-fold to 0.50-fold. |
| S5A CLR contrast | `scripts/new_new_script/05b_SuppTableS5A_S5B_class_logratio_contrasts.R` | **LPC flips sign, +0.264 to -0.138.** PG's interval now spans zero. Every other class shifts +0.035 through the CLR reference. |
| S5B ALR contrast | same script | Follows S5A. LPC goes -0.201 to -0.615. |
| S5C correlation delta | `scripts/new_new_script/05c_SuppTableS5C_clr_correlation_delta.R` | LPC-PC in CTL goes -0.412 to -0.130. PC-PS in CTL goes +0.445 to +0.433. 57 of 78 pairs moved. |

S20 and S21 were also rebuilt, by the new
`scripts/new_new_script/13a_SuppTableS20S21_overlap_gene_and_locus.R`. They had
come from the first-pass Python and still carried the LIN individual-lipid set
at 4,323 genes and 971 loci rather than 4,319 and 970.

## Scripts written in this pass

- `05b_SuppTableS5A_S5B_class_logratio_contrasts.R` -- S5A and S5B had **no
  generating script anywhere in the repository**. The definition was recovered by
  fitting against the shipped numbers, which it reproduces for ten of thirteen
  classes to within 0.005 on the effect and on both interval bounds.
- `05c_SuppTableS5C_clr_correlation_delta.R` -- lifted out of
  `_legacy_pipeline/28_class_logratio_stats.R`, which could not be re-run because
  the same pass rewrites `Figure2_OPLS_DA` and restores the OPLS VIP table that
  was removed from the supplement.
- `13a_SuppTableS20S21_overlap_gene_and_locus.R` -- S20, S21 and the
  gene-to-locus inflation table had no R generator.

## Bugs found and fixed along the way

Each of these silently produced wrong or unreproducible output.

1. **`06_Fig2_class_composition.R` and `06b_SuppTableS6_species_inventory.R` read
   the wrong lipid class file.** Both defaulted to
   `data/lipid_class/final_lipid_classes.csv`, the August annotation with the old
   home-grown vocabulary, rather than `data/metadata/final_lipid_classes.csv`
   with the eight LIPID MAPS categories. Run either with no environment variable
   and you got the superseded grouping.
2. **`03_Fig1_population_structure.R` and `04_SuppFig3_class_pca.R` wrote the
   pre-consolidation table numbers** (S25, S25a, S25b). Those are now S24, S24a
   and S24b, and S25 is the heritability table, so every run dropped a duplicate
   S25 into `table/supp` and broke the workbook builder.
3. **Three figure scripts wrote filenames the manuscript does not include.**
   `06c` wrote `SuppFig_S4_Lipid_Species_Counts.png` for a figure that prints as
   S5, `07` wrote `SuppFig_S6_Chemical_Space.png` for one that prints as S8, and
   `22g` wrote `fig/new_figures/Fig_Heritability_CTL_LIN.png` for S4. The
   manuscript compiled against whichever stale copy had last been moved by hand.
4. **`02_SuppFig2_QC_diagnostics.R` and `22g_Fig_heritability.R` had hard-coded
   absolute paths** into one user's home directory, so they ran on one machine
   and nowhere else.
5. **`06a_Fig2C_LION_input.R` asserted against its own output incorrectly.**
   `count.fields()` defaults to `comment.char = "#"` and the sample id row is
   `,#1,#2,...`, so the check saw 2 fields instead of 758 and failed on a file
   that was perfectly well formed.
6. **`13_SuppFig8_ctl_lin_overlap.R` drew from the stale Python overlap tables.**
   It now reads S20, S21 and S30 directly, so the figure cannot drift from the
   tables again.
7. **`20_SuppTable1to3_ratio_species_stats.R` read lyso species as diacyl.**
   `\bPC\b` matches the CTL spelling `PC(18:2/0:0)`, so CTL's lyso-PC species
   were summed into the PC pool while LIN's were summed into LPC. It also loaded
   twelve packages it never calls, any one of which stopped the script on a
   machine that lacked it.
8. **`SuppTable_S1` carried a `jackknife_stability` column** for an analysis the
   manuscript no longer reports, and the script still wrote two jackknife tables
   that had already been deleted.

## Not reproducible from the repository

`fig/supp/SuppFig_S7_PCA_Lipids.png` has no generating script. It is carried over
from the earlier pipeline and cannot currently be rebuilt from the data.
`table/Linex2/LION-enrichment.csv` comes from the LION web tool rather than from
a script; the input written by `06a` is unchanged in content this pass, so the
enrichment result still applies.

## A naming fix worth knowing about

The species inventory shipped a sheet called **By superclass**, with a
`SuperClass` column, and Supplementary Table S4 carried a `Level` value of
`Superclass`. The values in both have been the eight LIPID MAPS categories since
the annotation moved over, and the manuscript contains no occurrence of the word
superclass, because LIPID MAPS has no such tier. Both are now **Category**, and
the source file is `SuppTable_S6c_Species_by_Category.csv`.

One cosmetic mismatch is left deliberately. The text cites "Supplementary
Table S4, composition %TIC" and the sheet is named `Composition pctTIC`, because
a sheet name is easier to live with without a percent sign in it. Rename the
sheet in `15_SuppTables_build_workbooks.R` if you would rather they match
exactly.
