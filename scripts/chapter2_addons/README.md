# Chapter 2 add-on analyses

Two analyses added to address reviewer-facing gaps in the SoLD chapter.

## 1. Botanical race / population structure  (`race_analysis.py`)
Tests whether lipid-class composition differs by botanical race (`data/race/SAP_race.csv`,
`Original_Race` collapsed to Bicolor / Caudatum / Durra / Guinea / Kafir / Mixed) and by
genetic cluster (`K.Cluster`, K = 6), separately in CTL and LIN.

- Traits: total lipid signal (log10 TIC) + 13 class %TIC values
  (MGDG, PC, DG, DGDG, TG, PE, SQDG, MG, PG, PA, LPC, PS, LPE)
- Test: Kruskal-Wallis (tie-corrected), Benjamini-Hochberg within condition
- Effect size: epsilon-squared (proportion of rank variance explained)

Outputs: `fig/race/`, `table/race/race_structure_lipid_tests.csv`

## 2. CTL vs LIN GWAS candidate overlap  (`gwas_overlap.py`, `locus_overlap.py`, `fig_overlap.py`)
Formal test of whether the two trials' candidate sets overlap more than chance,
at gene level and after collapsing LD blocks into independent loci.

- Input: `data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv`
- Gene universe: 34,027 genes from `data/LD_mapped/genes_ranges/genes.range`
- Locus universe: genome tiled into 100 / 250 / 500 kb windows containing >= 1 gene
- Test: hypergeometric (Fisher right tail); Jaccard index reported alongside

Outputs: `fig/overlap/`, `table/overlap/`

## Running
`python3 race_analysis.py` then `python3 gwas_overlap.py; python3 locus_overlap.py; python3 fig_overlap.py`
(requires pandas, numpy, matplotlib; `stats_lite.py` supplies chi-square / hypergeometric /
BH functions so scipy is not needed. Paths are relative to the repo root.)
