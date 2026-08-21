# LD-aware validation of class-stratified GO-BP enrichment

## Why

LD-based candidate mapping assigns every gene model in a linkage block to one
association. Genes in a block are therefore not independent observations, and a
gene-level hypergeometric test treats one locus as several draws.

Two distinct problems follow, and they need different fixes.

**Problem A - genomic clustering of functionally related genes.** Broad GO terms
(translation, transmembrane transport, cell division) have genes that cluster in
gene-rich regions. Drawing N LD blocks picks them up at a far higher rate than
independent gene sampling predicts. This is a *statistical* problem and the
permutation test below fixes it.

**Problem B - single-locus support.** A term carried by 3 genes of one tandem
array is one observation, not three. It is statistically significant under both
tests and correctly so, but it is not gene-set-level evidence. No statistic fixes
this; it is a *reporting* matter, handled by the `blocks_carrying_term` column.

## What the script does

`scripts/new_script/29_GO_enrichment_LDaware_permutation.py`

The observed statistic stays at gene level, preserving the resolution of the full
GO background. Only the null changes. Each permutation redraws the same number of
linkage blocks as the observed candidate set, each contributing the same number of
genes as the corresponding observed block, sampled from background blocks matched
on gene density (quintiles of annotated genes per 250 kb window). 8,000
permutations; permutation P-values BH-adjusted across terms within layer.

Verified beforehand: the existing gene-level Fisher and BH implementations are
correct (Fisher matches an independent hypergeometric to 1e-11 over 400 random
cases; BH matches exactly).

## Result

All 26 individual-lipid and 22 class-sum/ratio terms that were BH-significant at
gene level remain significant under the LD-aware null. The correction changes the
**ranking and strength**, not the membership. The most affected terms are the
broad ones:

| Term | Genes | Blocks | Gene-level expected | LD-aware null mean | Gene-level q | LD-aware q |
|---|---|---|---|---|---|---|
| CTL TG positive regulation of cytoplasmic translation | 14 | 13 | 4.28 | 8.65 | 4.6e-03 | **0.044** |
| LIN SQDG positive regulation of cytoplasmic translation | 12 | 12 | ~4 | 7.17 | 6.9e-04 | **0.044** |
| LIN FA transmembrane transport (sum/ratio) | 31 | 22 | ~18 | 17.87 | 2.7e-02 | 0.011 |

The two cytoplasmic-translation terms move from the top of the gene-level ranking
to the significance boundary. Terms with narrow definitions and few background
genes are barely affected.

Terms resting on a **single** LD block (report these as one locus, not a gene set):

- individual: LIN PC sterol metabolic process (5 genes); LIN FA nitrate
  assimilation, detection of stimulus, and the three purine-transport terms
- sum/ratio: most of the LIN purine, nitrate and detection-of-stimulus terms;
  12 of the 22 significant results share the identical three genes
  `SORBI_3006G195000 / 195100 / 195200`, a tandem purine-permease array at
  chr6 ~54.79-54.80 Mb
- the CTL PA oxidative-stress terms rest on the tandem peroxidase array
  `SORBI_3005G011200 / 011300 / 011500` (2 blocks)

Terms with genuinely distributed support (>=3 independent blocks) include CTL TG
*L-tyrosine catabolic process* (4 genes / 3 blocks; LD-aware q = 0.007), CTL PC
*proteolysis* (5/4), LIN PC *cell division* (12/11), LIN SQDG *DNA topological
change* (3/3), LIN TG *regulation of salicylic acid biosynthetic process* (4/3),
and LIN AEG *sphingolipid metabolic process* (3/3).

## Outputs

- `GO_BP_enrichment_by_lipid_class_individual_LDaware.tsv`
- `GO_BP_enrichment_by_lipid_class_sumratio_LDaware.tsv`

Added columns: `n_ld_blocks_in_set`, `blocks_carrying_term`, `null_mean`,
`p_perm`, `q_perm_bh`.

## Annotation provenance (resolved)

Canonical GO source for all enrichment in Chapter 2:

    data/annotation/gene_annotation.txt   17,758 genes, 8,464 with >=1 GO_BP

This is the file every enrichment script already pointed at, so no results
required regeneration. A second, older table (`data/metadata/gene_annotation.txt`,
16,557 genes / 7,711 with GO_BP, last modified 2026-05-04) disagreed with it --
only 9,880 genes shared, GO_BP differing for 6,636 of those -- and was archived to
`data/metadata/_to_delete/` on 2026-08-18. The Shiny app (`app.R`) resolves its
gene-annotation path from an ordered candidate list with
`data/annotation/gene_annotation.txt` first, so it is unaffected.
