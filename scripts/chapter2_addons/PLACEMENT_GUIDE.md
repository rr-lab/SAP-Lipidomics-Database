# Where the new text goes

All new prose and tables are in `chapter2_new_sections.tex`, split into five labelled blocks.
Below: where each block goes, plus the edits needed to *existing* text so the chapter stays
internally consistent.

---

## New material

| Block | Type | Insert |
|---|---|---|
| 1 | Results subsection + Table + Figure | After **"Lipid variation across genotypes"**, before **"Lipidome overview for CTL and LIN SAP"** |
| 2 | Results subsection + 3 Tables + Figure | After **"LIN recurrent GWAS signals and annotation-based functional enrichment"**, before the **LINEX** subsection |
| 3 | Discussion paragraph | Append to the end of **"Integrated interpretation of lipidome remodeling across contrasting field trials"** |
| 4 | Discussion subsection | After **"Reaction-supported lipid-remodeling candidates under LIN implicate DG--TG cycling and phospholipid turnover"** |
| 5 | Two Methods subsubsections | After **"Genome-wide Association Studies (GWAS)"** |

New figure files (already generated):
- `fig/race/FigRace_C_variance_explained.png` — main figure for Block 1
- `fig/race/FigRace_A_CTL_boxplots.png`, `fig/race/FigRace_B_LIN_boxplots.png` — supplementary
- `fig/overlap/FigOverlap_CTL_vs_LIN.png` — main figure for Block 2

Suggested supplementary tables (CSV already written):
- `table/race/race_structure_lipid_tests.csv` — full KW output, both groupings, both conditions
- `table/overlap/gwas_overlap_by_class.csv` — per-class overlap, both trait layers
- `table/overlap/shared_genes_all_layers_ranked.csv` — the 354 shared genes, ranked
- `table/overlap/shared_genes_sumratio.csv` — the 18 genes / 2 loci recurrent in class sums+ratios
- `table/overlap/gwas_overlap_crossclass_jaccard.csv` — CTL class x LIN class similarity matrix

---

## Edits to existing text

### 1. Intro to the two GWAS results subsections — report loci alongside genes

Both GWAS subsections currently give only gene counts. Add the locus count on first mention so
the reader is not left with an unqualified "4,323".

- CTL, individual-lipid layer: after "identified 1,100 candidate genes" add
  *"(378 independent loci at 250 kb resolution; see Table~\ref{tab:inflation})"*.
- LIN, individual-lipid layer: after "identified 4,323 candidate genes" add
  *"(971 independent loci at 250 kb resolution)"*.
- CTL class-sum/ratio layer: 115 genes -> 36 loci.
- LIN class-sum/ratio layer: 812 genes -> 161 loci.

### 2. Discussion opening paragraph — the functional summary does not match the GO results

Current text reads:

> CTL-associated loci were enriched for genes involved in basal metabolism, chloroplast function,
> development, and cell-wall biology, whereas LIN-associated loci were linked to nutrient-status
> signaling, structural remodeling, defense responses, and cold acclimation.

None of "chloroplast function" or "cell-wall biology" appears in the CTL GO enrichment results,
which were *positive regulation of cytoplasmic translation*, *L-tyrosine catabolic process*,
*DNA replication initiation*, and (class-stratified) *hydrogen peroxide catabolic process* /
*response to oxidative stress*. The cell-wall material is in the commented-out subsection.
Suggested replacement:

> CTL-associated loci were enriched for processes reflecting basal metabolic and developmental
> state -- cytoplasmic translation, amino-acid catabolism, DNA replication initiation, and
> oxidative-stress and peroxidase activity -- whereas LIN-associated loci were enriched for
> nucleotide and purine transport, nitrate assimilation, transmembrane transport, and stimulus
> detection, alongside class-stratified enrichments for defense response, jasmonate signaling,
> and sphingolipid metabolism.

Then append, in the same paragraph:

> A direct comparison of the two candidate sets confirms that this divergence is not merely one
> of emphasis: after collapsing linkage blocks, the loci recovered under the two regimes overlap
> no more than expected by chance.

The same mismatch appears in the thesis Conclusion (`\section{Conclusion}` in the thesis build,
not present in `main.tex`); apply the same correction there, and replace the sentence beginning
"GWAS integration further showed that CTL-associated loci were enriched for genes linked to lipid
and isoprenoid metabolism, chloroplast function, growth, and cell-wall processes..." with the
corrected wording above plus:

> Formal comparison of the two candidate sets showed that, once linkage blocks were collapsed,
> the loci recovered under the two regimes overlapped no more than expected by chance, while
> class-level composite traits retained a modest but significant excess -- indicating that the
> genetic architecture of the sorghum leaf lipidome is strongly environment-dependent and that
> class-level traits are the more reproducible unit for cross-environment comparison. Lipid
> composition was, however, essentially independent of botanical race and of genome-wide genetic
> structure, supporting the interpretation that the reported associations are not artefacts of
> population stratification.

### 3. Heritability paragraph — reconcile with the new structure result

The current text reports near-zero genomic heritability in CTL (one feature above $h^2=0.20$)
immediately before reporting 1,100 CTL candidate genes, without addressing the tension.
Add after the heritability sentences:

> These heritability estimates are trial-specific and were computed on SERRF-normalised molecular
> features rather than on the class-level traits used for most downstream analyses, so they are
> not directly comparable to the trait set entering GWAS. They do, however, set an expectation of
> limited cross-trial reproducibility that the candidate-locus comparison below confirms directly.

If the $h^2$ method is not yet in Methods, it must be added there; it is currently absent.

### 4. Statement of scale for PS / LPC / LPE

Block 3 supplies this for the Discussion. Consider also softening the Results text where the four
signatures are first enumerated, e.g. by adding after the list:

> The SQDG and TG axes involve classes contributing several percent of TIC, whereas the PS and
> lysophospholipid axes describe redistribution within pools accounting for well under 0.2\% of TIC.

### 5. Repetition to trim

With Block 4 added, the four signatures are now stated in six places. Cut the restatements in the
LINEX subsection and in the Discussion opening; keep the enumerations in the lipidome overview,
in the "Integrated interpretation" subsection, and in the Conclusion.

---

## Numbers used, for checking

Race: 222 CTL / 211 LIN accessions with race; groups 17/16, 76/71, 22/23, 19/17, 7/7, 81/77.
Cluster: 387 CTL / 355 LIN. No race effect survives BH; min raw p = 0.0397 (PG, LIN, q = 0.468)
and 0.0648 (TG, CTL, q = 0.492). PS x cluster in LIN: H = 43.43, p = 3.02e-8, q = 4.23e-7,
eps^2 = 0.110; PS x cluster in CTL q = 0.206, eps^2 = 0.024. All other traits q > 0.20,
eps^2 < 0.021.

Mean %TIC: PS 0.014 (CTL) / 0.037 (LIN); LPE 0.0064 / 0.0114; LPC 0.096 / 0.127;
MGDG 34.5 / 32.1.

Overlap background: 34,027 gene models; 4,844 / 2,279 / 1,239 gene-containing windows at
100 / 250 / 500 kb. All other overlap numbers are in Tables 2 and 3 of the new text and in
`RESULTS_SUMMARY.md`.
