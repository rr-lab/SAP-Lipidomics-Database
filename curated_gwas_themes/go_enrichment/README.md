# GO_BP enrichment of GWAS candidate genes (CTL & LIN)

Over-representation analysis of GO biological-process terms among the sorghum lipid
GWAS candidate genes, combining the **individual** and **sum/ratio** trait layers
within each condition (CTL, LIN).

## Method
- **Test:** Fisher exact, one-sided (greater), Benjamini-Hochberg FDR.
- **Background:** all annotated genes carrying >=1 GO_BP term (8,464 genes).
- **Thresholds:** term background >=5 genes; term-overlap >=3 candidates (same as the
  original `run_combined_condition_go_enrichment.R`).
- **Recurrence:** per gene per condition = max `N_Phenotypes` across the two trait layers.
- **Loci:** contributing gene midpoints on the same chromosome are merged when their
  +/-25 kb windows overlap (midpoints within 50 kb), approximating the manuscript's
  window-collapse logic, to flag terms driven by a single physical locus (LD artifact).

## Headline result
On the **full** candidate sets there is **no GO_BP enrichment** (0 terms at q<0.05 in
either condition): CTL tests 1,492 genes (17.6% of background), LIN 3,678 (43.5%) — too
close to genome-wide to enrich. Enrichment appears **only after filtering to recurrently
associated genes**; the ~5% recurrence cut is the practical sweet spot (the 0.5% cut is
underpowered). See `significant_terms_locus_breakdown.tsv` for which surviving terms are
real multi-locus signals vs single-locus artifacts.

## Files
- `CTL_GO_BP_enrichment_full.tsv` / `LIN_GO_BP_enrichment_full.tsv` — full-set ranked results
  (all tested terms). Gene-list column dropped for size; regenerate with the script for gene lists.
  (The LIN full table on disk is abbreviated to the top rows + a pointer; run the script for the
  complete 544-row table.)
- `CTL_GO_BP_enrichment_by_recurrence.tsv` / `LIN_..._by_recurrence.tsv` — every tested term at
  each recurrence cut (top10 / top5 / top1), with overlap gene lists.
- `significant_terms_locus_breakdown.tsv` — for every q<0.05 term: #independent loci,
  biggest single-locus cluster, verdict (SINGLE LOCUS / mostly 1 locus / multi-locus), gene list.
- `go_enrichment_recurrence.py` — self-contained; regenerates all of the above. `python3 go_enrichment_recurrence.py`

## Caveats
Most surviving terms rest on only 3 overlapping genes; single-locus terms
(e.g. CTL "carbohydrate transmembrane transport", CTL "sterol metabolic process",
LIN "protein ubiquitination") should be treated as one association, not pathway-level enrichment.
