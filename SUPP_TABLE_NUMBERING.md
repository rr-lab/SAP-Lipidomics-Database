# Supplementary table citations vs. the workbook

Every `Supplementary Table~SN` citation in `main.tex`, checked against the sheets in
`table/supp/SupplementaryTable_SoLD.xlsx` (22 sheets, rebuilt 2026-08-25).

**7 of 21 citations resolve correctly. 7 point at a sheet holding different content.
7 point at nothing at all.**

Nothing has been changed in `main.tex`.

---

## A. Correct as they stand — leave alone

| Line | Citation | Content claimed | Sheet |
|---|---|---|---|
| 158 | S4 | Top-variance lipids | `S4_TopVariance_Lipids` |
| 217 | S5A | Class-level CLR contrasts | `S5A_Class_CLR_Contrast` |
| 217 | S5B | TG-referenced ALR contrasts | `S5B_Class_ALR_Contrast` |
| 217 | S1 | Lyso-lipid ratio behaviour | `S1_Ratio_Statistics` |
| 217 | S6a to S6c | Species counts by class/superclass | `S6a` / `S6b` / `S6c` |
| 233 | S5C | Class-pair correlation deltas | `S5C_Class_CLR_Corr_Delta` |
| 164 | S3 and S4 | Species stability, top-variance | `S3` / `S4` (see note) |

Note on line 164 — S3 and S4 cover the stability and variance halves, but the clause
"genetic tractability of this variation differed sharply between trials" is a
heritability claim, and no heritability table is cited or present. See B2.

---

## B. Wrong sheet — the content exists, the number does not point to it

| Line | Citation | Content claimed | That number holds | Should be |
|---|---|---|---|---|
| 129 | S1 | QC-RSD before/after SERRF | Ratio statistics | **No sheet.** QC-RSD data is not in the workbook. Either add it or drop the citation, since Fig.~1C already shows it |
| 158 | S5 | Genomic heritability, $h^2>0.20$ | CLR/ALR/$\Delta r$ contrasts | New sheet from `table/new_table/SuppTable_Genomic_Heritability_*.csv` |
| 160 | S6 | Kruskal--Wallis race and cluster tests | Species counts | `Race_Structure_Tests` sheet (currently unnumbered) |
| 217 | S6 | PS tracking genome-wide cluster | Species counts | `Race_Structure_Tests` |
| 241 | S13 and S14 | CTL candidate genes, both layers | LINEX GWAS support | **S7 and S8** |
| 241 | S15 and S16 | LIN candidate genes, both layers | Lipid classes; S16 absent | **S9 and S10** |
| 364 | S15 | No candidate for branch (iii) | Lipid classes | **S14** (`LINEX_GWAS_BranchSummary`) |
| 551 | S29 | Lipid classes and subclasses (Lipid Maps) | absent | **S15** (`final_lipid_classes`) |

The S13/S14 and S15/S16 pair at line 241 is the most damaging. A reviewer following
the candidate-gene counts lands on the LINEX tables and the annotation table instead
of the GWAS results those counts came from.

---

## C. No sheet exists — data is on disk but was never added to the workbook

| Line | Citation | Content claimed | Where the data actually is |
|---|---|---|---|
| 274, 276 | S17 and S18 | GO-BP and GO-MF, all terms | `table/go_enrichment/Table_GO_enrichment_all.tsv` |
| 274 | S19 | Genes carrying each enriched term | `table/go_enrichment/genes_in_enriched_GO_terms.tsv` |
| 276 | S20 | Terms recovered by recurrence rather than LD interval | Not located. Needs identifying |
| 300 | S21 to S24 | CTL/LIN overlap, gene and locus level | `table/overlap/gwas_overlap_overall.csv`, `gwas_overlap_locus_level.csv`, `shared_genes_*.csv` |
| 364 | S25 to S28 | LIN GWAS candidate tables for the branches | Duplicates S9/S10, or means S13/S14. Ambiguous |

---

## D. Proposed scheme

The workbook numbering is already coherent and content-grouped, so the cheaper fix is
to keep it and move the citations, adding sheets for section C.

| Sheet | Content | Status |
|---|---|---|
| S1 to S6c | as built | keep |
| S7 to S10 | GWAS candidate genes, four layers | keep |
| S11 to S14 | LINEX reactions, balance, gene support, branches | keep |
| S15 | Final lipid classes | keep |
| **S16** | Genomic heritability | add from `table/new_table/` |
| **S17** | Race and cluster Kruskal--Wallis tests | rename `Race_Structure_Tests` |
| **S18** | GO enrichment, all terms, both ontologies | add from `table/go_enrichment/` |
| **S19** | Genes in enriched GO terms | add from `table/go_enrichment/` |
| **S20** | Recurrence-selected GO terms | add once located |
| **S21** | CTL/LIN overlap, gene level | add from `table/overlap/` |
| **S22** | CTL/LIN overlap, locus level | add from `table/overlap/` |
| **S23** | OPLS-DA VIP scores | rename `OPLS_VIP_ratios`, resolves the S7 filename clash |
| **S24** | SERRF scan-to-lipid audit | add from `table/new_table/` |

### Citation edits that follow

| Line | From | To |
|---|---|---|
| 129 | S1 | drop, or add a QC sheet and cite it |
| 158 | S5 | S16 |
| 160 | S6 | S17 |
| 217 | S6 | S17 |
| 241 | S13 and S14 | S7 and S8 |
| 241 | S15 and S16 | S9 and S10 |
| 274, 276 | S17 and S18 | S18 |
| 274 | S19 | S19 (unchanged, once the sheet exists) |
| 276 | S20 | S20 (unchanged, once located) |
| 300 | S21 to S24 | S21 and S22 |
| 364 | S25 to S28 | S9 and S10, or S13 and S14 — needs your call |
| 364 | S15 | S14 |
| 551 | S29 | S15 |

Adding the section C sheets means extending the `SHEETS` list in
`scripts/pipeline/build_supplementary_workbook.py`. Nothing else in that script changes.

---

## E. Two questions only you can answer

1. **Line 364, S25 to S28.** "Cross-referencing branches with the LIN GWAS candidate
   tables" reads like the raw LIN candidate tables (S9, S10), but `S13_LINEX_GWAS_GeneSupport`
   is literally the output of that cross-reference. Which did you mean?
2. **Line 276, S20.** Terms recovered by recurrence rather than by LD interval. I could
   not find a table producing this. Does it exist, or was the claim made from an
   analysis that was never written out?
