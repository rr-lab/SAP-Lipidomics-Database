# Supplementary tables, S1 to S31

One table per file, built by `scripts/new_new_script/15_SuppTables_build_flat.R` on 2026-09-30.

Files live in `table/supp/tables/` (descriptive names) and are copied to
`final_submission/supplementary_tables/` under the journal's `SN_Table.xlsx` names.

`was` is the number the table carried in the 41-table scheme that preceded the ten
workbooks. `line` is the first line of `main.tex` that cites it.

| # | Title | was | rows | cols | line | cites |
|---|---|---|---|---|---|---|
| S1 | Population structure group sizes | S24a | 12 | 7 | 134 | 1 |
| S2 | Population structure lipid tests | S24 | 56 | 23 | 134 | 1 |
| S3 | Population structure PC tests | S24b | 8 | 9 | 134 | 1 |
| S4 | Heritability per species structure | S28 | 430 | 21 | 134 | 2 |
| S5 | Heritability per species paired | S25 | 328 | 8 | 136 | 1 |
| S6 | Heritability structure conditioned | S27 | 26 | 14 | 136 | 1 |
| S7 | Heritability class sums | S26 | 34 | 8 | 136 | 2 |
| S8 | Heritability ancestry robustness | S29 | 54 | 6 | 136 | 1 |
| S9 | Frozen lipid species set | S15 | 316 | 6 | 144 | 4 |
| S10 | Top variance lipids | S4 | 20 | 5 | 244 | 1 |
| S11 | Class composition pctTIC | S5D | 40 | 4 | 248 | 2 |
| S12 | Class contrasts CLR and ALR | S5A + S5B | 25 | 7 | 248 | 3 |
| S13 | Class CLR correlation delta | S5C | 78 | 10 | 250 | 3 |
| S14 | Lipid ratio statistics | S1 | 33 | 11 | 471 | 2 |
| S15 | Species counts by class and category | S6b + S6c | 25 | 4 | 244 | 1 |
| S16 | Composition stability | S5G | 26 | 10 | 248 | 1 |
| S17 | LION enrichment | S5E | 10 | 7 | 267 | 2 |
| S18 | Chemical space | S5F | 26 | 5 | 274 | 1 |
| S19 | GWAS candidates CTL individual | S7 | 1062 | 15 | 280 | 1 |
| S20 | GWAS candidates CTL sum ratio | S8 | 54 | 15 | 280 | 2 |
| S21 | GWAS candidates LIN individual | S9 | 4370 | 15 | 280 | 2 |
| S22 | GWAS candidates LIN sum ratio | S10 | 385 | 15 | 280 | 2 |
| S23 | GO loci collapsed | S16 | 217 | 15 | 369 | 5 |
| S24 | GO BP all terms | S17 | 164 | 13 | 514 | 2 |
| S25 | GO MF all terms | S18 | 180 | 13 | 440 | 3 |
| S26 | GO genes in enriched terms | S19 | 919 | 17 | 466 | 1 |
| S27 | Overlap by lipid class | S30 | 17 | 12 | 385 | 1 |
| S28 | Shared candidate genes | S31 + S22 | 291 | 9 | 389 | 1 |
| S29 | LINEX reactions and balance | S11 + S12 | 4 | 18 | 436 | 5 |
| S30 | LINEX GWAS gene support | S13 | 15 | 14 | 440 | 1 |
| S31 | LINEX branch summary | S14 | 5 | 8 | 440 | 1 |

## Dropped on the way from 39 sheets to 31 tables

| was | why |
|---|---|
| S23 Shared sum-ratio | 0 rows. The finding is already the sentence "shares no gene at all between trials, against 0.6 expected". |
| S6a Species summary | 6 rows, the counts 214 / 216 / 164 / 50 / 52 / 266. All six are in the Results text and in SuppFig S5A. |
| S20 Overlap gene level | 3 rows, identical to rows 1-3 of main-text Table 3. |
| S21 Overlap locus level | 9 rows, six of which are Table 3's lower block. |

## Merged

| new | from | how |
|---|---|---|
| S12 | S5A + S5B | stacked, `Scale` column holds CLR or ALR |
| S15 | S6b + S6c | stacked, `Tier` column holds Class or Category |
| S28 | S31 + S22 | joined on GeneID. The 282 individual-layer genes are a subset of the 291 shared, so the four per-trial columns are NA for the 9 cross-layer genes. |
| S29 | S11 + S12 | joined on reaction branch, one row per branch |

## Known order exceptions

Numbers follow first mention, except where a table is only cited in the Discussion or
Methods and is kept with the block it belongs to rather than moved to the back.
S3 trails S4 by one sentence; S14, S16, S24 and S26 are cited later than their
block-mates.
