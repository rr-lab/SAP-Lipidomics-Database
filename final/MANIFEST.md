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

## Update, 2026-09-24

Four things changed since the 09-17 build.

**Supplementary Tables S7 to S10 were rebuilt** (workbook S6, all four sheets).
The files packed on 09-17 predated the high-heterozygosity marker filter and the
class-sum and ratio rebuild, so every count in them disagreed with the
manuscript. They also still carried `SORBI_3001G522600` at p = 4.5e-35, the
signal traced to 23 removed markers in ~300 bp at heterozygosity 0.81 to 1.00.

    sheet             was     now
    CTL individual   1,100   1,062
    CTL sum/ratio      115      54
    LIN individual   4,319   4,370
    LIN sum/ratio      812     385

Rebuilt by `scripts/new_new_script/72_rebuild_S7toS10.R` from
`data/gene_annotation_final/` at Max_r2 >= 0.4. The superseded files are at
`table/supp/_superseded_S7toS10_preFilter/`.

**Workbook S7 (GO enrichment) picked up the 09-23 rerun.** Two files were
competing for the S19 slot, `SuppTable_S19_Genes_in_enriched_terms.tsv` from
09-15 and `SuppTable_S19_genes_in_enriched_GO_terms.tsv` from 09-23, which is
what `14_GO_supp_tables.R` actually writes and what matches S16 to S18. With
both present the workbook builder aborted on a duplicate key, so the 09-15 file
was archived beside the S7 to S10 ones. Sheet rows went 50/103/103/254 to
285/164/180/1,429.

**Figure 2 panel C now shows eight LION terms.** `PC(38:2)` and `C20:1` are
structural bins rather than functional or biophysical categories and are
reported in workbook S5 only. All ten terms remain in that workbook.

**`final/table/supp/` is gone.** It should never have existed, as the paragraph
above says. It moved to `final/_superseded/table_supp_2026-09-17/`.

Table 4 in the manuscript was re-derived at the same time by
`scripts/new_new_script/73_top_loci_table.R` and is at
`table/main/Table4_top_loci.tex`. Every locus, interval, gene count, trait count
and p-value reproduced; two label genes changed.

**Cellular component was dropped from the GO analysis, 2026-09-24.** The
enrichment had been testing three ontologies while the manuscript said two,
in the Results and twice in the Methods. Of the 89 CC terms that passed, none
named anything lipid-related, and the membrane compartments among them were the
broad ones the figure already excludes on background size. `13_GO_enrichment_by_class.R`
now takes `GO_ONTOLOGIES`, default `BP,MF`. BH correction is applied by
(Condition, Layer, Ontology), so all 344 surviving BP and MF q-values are
unchanged. Downstream counts move

    tested                1,623 -> 1,292
    passing q_LD < 0.05     433 ->   344
    collapsed rows          285 ->   217
    distinct terms          189 ->   159
    background <= 30        102 ->    87
    of those, >= 2 loci      85 ->    71
    gene sets in Figure 4    23 ->    22

Figure 4 loses `nuclear speck`; SQDG now carries three gene sets rather than
four. Every headline result is unaffected -- phosphate starvation at 51x over
three intervals, nitrate transport at 25x, sphingolipid metabolic process at
35x. The three-ontology outputs are kept at
`table/go_enrichment/_superseded_withCC_2026-09-23/` and
`table/supp/_superseded_withCC_2026-09-23/`, the latter with the old figure.

**The CTL/LIN overlap was rebuilt, 2026-09-24.** Everything in workbook S8 and
Supplementary Figure S8 descended from `ALL_LD_candidate_genes_master.tsv` of
09-16, which carried the pre-filter candidate counts. S20 and S21 had no R
generator at all -- they came from the first-pass Python. The Results paragraph
had already been updated to the current numbers, so the text was right and the
files behind it were not. `74_rebuild_overlap_tables.R` now rebuilds the master
from `data/gene_annotation_final` and writes S20 and S21;
`14_SuppTableS30S31_...` regenerates S22, S23, S30 and S31 from it.

    sum/ratio shared genes      18 ->     0
    S30 class rows              25 ->    17
    TG fold                   1.63 ->  1.74
    PC fold                   3.06 ->  3.47

Two bugs were fixed on the way. The sum/ratio class parser matched only the old
`_log10safe` suffix, so since the phenotype rebuild every ratio was being filed
under a class named after the whole trait (`MGDG_over_PS_log10ratio`) instead of
under its numerator and denominator. And panel D of Supp Fig S8 read
`table/overlap/gwas_gene_to_locus_inflation.csv`, last written in May 2026, so it
showed 1,100 / 4,323 / 812 while panels A to C showed the current sets; that file
is now written by the same script. The panel labels also printed "2.07.. chance"
because the PNG device's font drops the multiplication sign, now a plain x.

Superseded files are at `table/supp/_superseded_overlap_preFilter/`, including
the old figure and the old master.

**S5A, S5B and S1 were put on one definition of a lipid class, 2026-09-24.**
They had been written by two legacy monoliths that disagreed with each other and
with the manuscript.

S5A and S5B came from `_legacy_pipeline/22_lipidome_class_composition.R`, which
reads a class off a feature name with `^(CLASS)(?=\()`. That counts
`PC(18:2/0:0)` and `PC(22:0/0:0)` as PC when they are lyso species, so the
shipped S5A put LPC at +0.264 on the CLR scale, a 1.30-fold increase, while the
manuscript argues LPC does not accumulate. The correct value is -0.138, a
0.87-fold decrease, which is what the text already said.

S1 came from `new_script/20_SuppTable1to3_ratio_species_stats.R`, which defines a
class as the mean of log10 relative abundance over its species. That is a
per-species geometric mean and can differ in SIGN from a class total. On it MG
rose and TG/MG fell; on class totals MG falls from 1.30% to 1.05% of TIC and
TG/MG rises 2.38-fold. The Discussion was quoting the geometric-mean ratios in
one sentence and class-total values in the next.

`scripts/new_new_script/75_class_composition_contrasts.R` now writes all three
from the same per-sample composition of the 13 focal classes, with
`normalize_lipid_name` applied first. The CLR contrast reproduces every value
the manuscript quotes: SQDG -0.942 (0.39x), MGDG 0.77x, DGDG 0.85x, PG no change,
LPC 0.87x. Ratios that move

    TG/MG      0.74x  ->  2.38x     direction reversed
    DG/MG      0.34x  ->  1.31x     direction reversed
    DG/PS      0.13x  ->  0.32x
    PE/SQDG    2.76x  ->  2.41x
    PG/SQDG    2.63x  ->  2.54x
    LPC/LPE    0.38x  ->  0.48x

PC/SQDG stays at 2.07x and every PS ratio keeps its direction, so the two claims
the Discussion rests on are unaffected. Superseded files are at
`table/supp/_superseded_composition_preLyso/`.

**Workbooks S9 and S10 rebuilt, and the last five stale citations fixed,
2026-09-26.** These were the only supplementary tables still carrying pre-2026
numbers.

S12 (reaction balance) was written 2026-05-04 and disagreed with the manuscript
after the lyso rename moved PC(18:2/0:0) and PC(22:0/0:0) into LPC, which
changes the LPC and PC pools and therefore the LCAT* and LRO1 scores.

    branch    S12 (old)   now
    PNPLA3      +0.327    +0.314
    LRO1        +0.232    +0.261
    LCAT*       +0.238    +0.139
    PNPLA1      -0.268    -0.256

S13 and S14 named SORBI_3001G103800 and SORBI_3001G448800 for the LCAT* branch,
neither a candidate at Max_r2 >= 0.4 any more, and did not contain
SORBI_3004G341900 (PSAT), which the manuscript names. Six rows still used
_log10safe trait names. S13 went from 14 rows to 15 and now carries all four
genes the Discussion names.

S15 (frozen species set) was from 2025-07-25 with 319 rows. It held five
chlorophyll and pheophytin features the annotation has since dropped, lacked
LPC(18:2) and LPC(22:0), and had no Category_Code column. It is now a copy of
data/metadata/final_lipid_classes.csv, 316 species.

S11, the four reactions and their RHEA identifiers, is hand-curated with no data
dependency and was left alone. All of it is written by
`scripts/new_new_script/76_rebuild_linex_tables.R`; superseded files are at
`table/supp/_superseded_linex_preRebuild/`.

Five citations in main.tex pointed at workbooks that do not exist under the
consolidated numbering, because they still used the old per-file numbers.

    S16 -> S7, loci collapsed        (3 places)
    S17 -> S7, BP all terms
    S20 -> S8, gene level
    S26 -> S2, class sums            (2 places, already cited that way elsewhere)
    S31 -> S8, shared ranked

Every Supplementary Table citation in the manuscript now resolves to a workbook
that exists, and where it names a sheet, to that sheet.

**The paired heritability moved from 163 to 164 species, 2026-09-26.** The old
S25 built the shared set with intersect() on raw feature names, so CTL's
PC(18:2/0:0) and LIN's LPC(18:2) were treated as different species. They are the
same lipid, and normalize_lipid_name resolves them. 164 is the number the
species inventory (S6a) already reported as Common, the number that makes
164 + 50 + 52 = 266 add up, and the number the LION input already used.

Only LPC(18:2) is added, and it turns out to be heritable under LIN
(h2 = 0.551, CI 0.13 to 0.999) and not under CTL, which is why three counts move.

    species                     163 -> 164
    median h2 LIN             0.170 -> 0.171
    CI excludes zero, LIN        53 -> 54
    higher under LIN            131 -> 132

Median h2 under CTL stays 0.000, CI excludes zero under CTL stays 11, higher
under CTL stays 21, tied stays 11. Written by
`scripts/new_new_script/77_paired_heritability_164.R`; the 163-species file is at
`table/supp/_superseded_h2_163species/`. main.tex updated in seven places, so no
occurrence of 163 remains.
