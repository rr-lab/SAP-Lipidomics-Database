# Proposed responses to the Reviewer 2 pass

Nothing below has been applied to `main.tex`. Line numbers are from the current file.

Items are grouped by what they need from you. Three still need an answer only you have
— **A2** (what the LD mapping actually did), **D2** (four permutation details), and
**D3** (field, extraction, and LC-MS numbers). The GEMMA driver and LD-mapping code are
not in either repo, so I cannot read those off the scripts.

A1 is resolved — you confirmed three PCs.

---

## A. Hard contradictions — fix before anyone else reads it

### A1. GWAS covariates — RESOLVED, three PCs

**Confirmed by you: three genotype PCs were fitted.** So Results (line 237) is correct
and Methods (line 591) is the sentence that is wrong.

**Delete line 591's sentence**

> GWAS analysis was conducted on the first two PCs for each class.

**Replace with**

```latex
GWAS was performed in GEMMA using a univariate linear mixed model with a centered
kinship matrix and the first three genotype principal components fitted as fixed
covariates.
```

Line 237 needs no change. This one is ready to apply whenever you say go.

### A2. Candidate-gene mapping, three different definitions (lines 237, 605, and the 250 kb intervals)

**Results, line 237** maps significant SNPs to genes "using the LD-based cutoff of 0.4."

**Methods, line 605**

> For each marker, a 50~kb segment was designated, spanning 25~kb on either side, and all gene models within this area were retrieved.

**Enrichment** then uses non-overlapping 250 kb intervals.

That is three rules, and line 605 looks like text left over from before you moved to
the LD approach. This is the reviewer's strongest catch and I agree with it completely.

**You need to confirm what the LD mapping actually did.** The most likely version,
given how the candidate tables look, is below — but do not let me write this into the
paper unless it matches your code.

```latex
Bonferroni-significant SNPs were grouped into association regions by linkage
disequilibrium with the lead SNP ($r^2 \geq 0.4$). Candidate genes were defined as
all annotated \textit{S. bicolor} v3.1 gene models overlapping an association region.
For locus-level summaries and LD-aware enrichment, candidate genes were then collapsed
into consecutive, non-overlapping 250~kb genomic intervals, which serve as the unit of
independence throughout.
```

If a $\pm$25 kb flank was genuinely applied on top of the LD regions, add to the second
sentence — "...overlapping an association region, extended by 25~kb on either side."

Whatever the answer, line 605 as written has to go. Three questions I need from you.

1. Was the r² 0.4 threshold computed against the lead SNP, or pairwise among all significant SNPs?
2. Was a fixed bp window applied at all in the final pipeline, or is the 25 kb sentence a leftover?
3. Where does the LD-mapping script live? It is not in `SAP-Lipidomics-Database` or `SoLD_paper`, and the README concedes GEMMA was run outside the repo. For a resource paper this needs to be archived.

---

### A3. SpATS output, residuals or residuals plus mean (line 544)

The two sentences contradict each other back to back.

> We utilized the residuals, defined as the difference between observed intensity and the fitted spatial trend, as our analysis-ready spatially corrected intensity data, which will be used for \%TIC, CLR, and log-ratio. We used SpATS-adjusted intensities (residuals + fitted mean) to preserve scale and avoid negative values.

The second sentence is almost certainly the true one, since %TIC and CLR both require
positive values and raw residuals are signed. Replace both with —

```latex
For downstream analyses we used SpATS-adjusted intensities, calculated as the model
residuals plus the fitted mean. This preserves the original intensity scale and avoids
negative values, which is required for the \%TIC, CLR, and log-ratio transforms.
```

---

### A4. LINEX branches, one is not from LINEX (lines 353 and 615)

`SuppTable_S11_LINEX_Reactions.csv` says, for the LCAT row, verbatim —

```
"R1","LCAT","L_FAdelete","PC, DG","LPC, TG","RHEA:32843, RHEA:44236","Hypothesized based on reaction specificity (not from LINEX)"
```

Meanwhile the manuscript calls all four branches "LINEX-supported" in the Figure 6
caption (353) and in Methods (615). Your own supplement contradicts your main text.
This is the item I would fix first, because a reviewer who opens S11 loses trust in
everything else.

**Figure 6 caption, line 353**

```latex
\caption{\textbf{Candidate lipid reaction remodeling between LIN and CTL conditions.}
```

Then, in the caption body, after the four branches are introduced, add —

```latex
Three branches were retrieved from the LINEX2 reaction network; the LCAT-like branch
(marked LCAT*) was defined manually from reaction specificity and is not a LINEX
output.
```

**Methods, line 615**

```latex
For branch-level analysis we extracted candidate reactions linking DG/MG, TG, and
lysophospholipid pools. Three were taken from the LINEX2 network. The fourth, an
LCAT-like acyl-editing branch denoted LCAT*, was defined manually from reaction
specificity (RHEA:32843, RHEA:44236) rather than retrieved from LINEX2, and the
asterisk marks it as such throughout.
```

That last clause also solves the reviewer's separate complaint that `LCAT*` is never
defined and reads like a footnote marker.

---

## B. Overclaims

### B1. The worst one (line 497)

> ...these jointly implicate phospholipid head-group hydrolysis, TG synthesis and TG turnover among the LIN candidates, the pathway by which nitrogen limitation drives TG accumulation at the expense of phospholipids.

You never isolated nitrogen from phosphorus, planting date, year, or LC-MS run context.
As written this claims a demonstrated causal mechanism for a single nutrient.

```latex
...these jointly implicate phospholipid head-group hydrolysis, TG synthesis and TG
turnover among the LIN candidates, a plausible route for TG accumulation at the
expense of phospholipids under the combined LIN field context.
```

### B2. Sulfur was never imposed (lines 379, 411, 413)

You reduced N and P. Sulfur was not a treatment variable, so "acquisition of the three
macronutrients likely limiting the low-input trial" asserts something you did not test.

**Section header, line 411**

```latex
\subsection*{Under LIN, candidates implicate nitrate transport, phosphate signaling and plastid sulfur supply}
```

**Line 413, first sentence**

```latex
The main theme among LIN candidates is nutrient transport and status signaling,
spanning nitrate transport, phosphate signaling, and plastid sulfate supply relevant
to SQDG metabolism.
```

**Line 379, Discussion roadmap** — change "converge on nitrogen, phosphorus, and sulfur
acquisition" to "converge on nutrient transport and status signaling."

### B3. SULTR3 wording (line 417)

> ...yet SQDG was depleted under LIN, consistent with restricted plastidic sulfate supply. Both loci remain correlative and untested in sorghum.

You did not measure sulfate, cysteine, glutathione, or SULTR3 expression. The existing
caveat sentence helps, but the first clause still states the mechanism.

```latex
...yet SQDG was depleted under LIN, which makes restricted plastidic sulfate supply a
plausible hypothesis for follow-up rather than a demonstrated mechanism. Both loci
remain correlative and untested in sorghum.
```

### B4. Reaction-balance p-values (line 348)

I do **not** agree with the reviewer that this analysis is circular — the score is a
directional contrast between substrate and product classes, not a re-test of the
quantity that defined them. But leading with $p<10^{-56}$ invites the accusation for no
gain, and direction is the actual result.

```latex
All four branches shifted strongly and consistently in the expected direction (paired
genotype-matched Wilcoxon, BH-adjusted $p<10^{-56}$), indicating that the
reaction-balance summaries recapitulate class-level lipid remodeling rather than
measuring flux.
```

---

## C. Typos and small wording — no judgment needed

| Line | Current | Change to |
|---|---|---|
| 222 | `Arrows indicates direction` | `Arrows indicate direction` |
| 495 | `PS-centerd phospholipid redistribution` | `PS-centered phospholipid redistribution` |

The other formatting bugs the reviewer listed — `Table S19 Table`, `S12 Fig ; S13 Fig`,
`Supplementary S10 FigA`, `bioRχiv` — **do not exist in the current source.** Line 270
reads `Table~S19` correctly. They were reading an older PDF. No action.

### C1. The "seven GO-BP terms" count (line 274)

The reviewer is right that this reads as a miscount, because *transmembrane transport*
is then listed three times and a reader counting entries gets nine.

```latex
Seven unique GO-BP terms met both criteria, represented by nine condition/class-specific
enrichments because \textit{transmembrane transport} was enriched in three sets.
```

Then the existing list follows unchanged. Check my count of nine against the figure
before applying — I counted the entries in the sentence, not the underlying table.

---

## D. Gaps that need new text from you

### D1. GO-MF provenance (line 609)

The Methods paragraph describes GO-BP only, but the Results lean on GO-MF just as hard.
Add after the first sentence —

```latex
GO molecular-function terms were taken from the same PANTHER-derived annotation.
GO-BP and GO-MF were analyzed separately, each against a background restricted to
sorghum genes carrying at least one annotation in that ontology, and $p$-values were
BH-corrected within each ontology and analysis.
```

That single addition also delivers the point your Results already make about the two
ontologies being corrected separately.

### D2. LD-aware permutation detail (line 609)

Your paragraph is better than the reviewer credits — it already states 8,000 draws,
250 kb intervals, gene-density matching, and interval-count preservation. Four things
are genuinely missing, and I cannot invent them.

1. How were gene-density bins defined (quantiles, fixed counts, how many bins)?
2. Sampling with or without replacement within a draw?
3. What happened when a density bin had too few eligible intervals to match?
4. Why 250 kb, and why the ≥3-interval reporting threshold?

Give me those four and I will write two sentences that close it. Items 3 and 4 are the
ones a methods-minded reviewer will actually press on.

### D3. Field, extraction, and LC-MS Methods

I agree with the reviewer here and it is the largest writing task in the list. For a
paper whose central claim is cross-year analytical comparability, "standard extraction
protocols explained in Barnes et al. 2022" is not defensible. The minimum additions are
planting and sampling dates for both trials, N and P rates and fertilizer form,
pesticide/herbicide status (currently in the Introduction but absent from Methods),
field design and replication, tissue and growth stage, plus tissue mass, solvent system,
and internal standards for extraction, and ionization mode, column, mobile phases,
MS/MS acquisition, mass range, and QC injection frequency for the LC-MS. I can draft
this once you give me the numbers.

---

## E. Not a text problem — stale figures

The reviewer's headline "fatal" finding, 243 vs 244 species, is not an inconsistency in
your writing. `main.tex` says 243/152/49/42 at lines 88, 141, 151, and 495 and never
says otherwise. The 244 they saw is inside `SuppFig_S4_Lipid_Species_Counts.png`, which
is dated **Aug 21 17:10**. You removed the two artifact features from the master data on
**Aug 23 14:52**. I counted the current matrices — CTL 201, LIN 194, union 243, common
152, which is exactly what the text claims.

**Do not change the text to 244.** Re-run the figures.

These were all generated before the data edit and still show the old lipid set.

| Figure | File | Generated |
|---|---|---|
| Figure 3 | `fig/main/Figure1_Lipidomics_Landscape.png` | Aug 21 17:10 |
| Supp S3 | `fig/supp/SuppFig_S3_Compositional_Contrasts.png` | Aug 21 17:10 |
| Supp S4 | `fig/supp/SuppFig_S4_Lipid_Species_Counts.png` | Aug 21 17:10 |
| Supp S5 | `fig/supp/SuppFig_S5_PCA_Lipids.png` | Aug 21 17:10 |
| Supp S6 | `fig/supp/SuppFig_S6_CLR_Correlations.png` | Aug 21 17:10 |
| Supp non-focused | `fig/supp/SuppFig_S8_NonFocused_Lipid_Class_Context.png` | Aug 21 17:26 |

Regenerating them means

```
Rscript scripts/pipeline/22_lipidome_class_composition.R
Rscript scripts/pipeline/28_class_logratio_stats.R
```

Figure 4 and Figure 6 are already current (Aug 24 12:15 and 11:23), so the reviewer's
"Figure 4 is nearly unreadable" was aimed at the pre-restyle version and should be
re-judged.

---

## F. Where I would push back on the reviewer

**Cutting Figure 4 to eight panels.** The figure's only job is showing that one locus
recurs across multiple phenotypes. Two panels per locus makes that an anecdote rather
than a pattern. Keep sixteen and give it a full page, or enlarge the panels — do not
halve the evidence to fix a legibility problem that is really a sizing problem.

**"Circular" reaction-balance scores.** Addressed in B4. Fix the emphasis, reject the framing.

**Moving the Shiny section to supplement.** Depends entirely on where you submit. Plant
Genome publishes resource papers and the app is a genuine contribution. I would keep it
in the main text and cut its length rather than relocate it. The triterpenoid example
being off-theme is a fair hit though — trimming that example is the cheaper fix.

---

## G. Things the reviewer missed

- **Staleness as a class of problem.** Audit every figure against the Aug 23 data edit before submitting, not just the ones flagged above.
- **`SuppTable_S12_ReactionBalance.csv` has no producing script.** The numbers are recomputed inside `21_figure6_linex_full.R` for panel B but never written out.
- **Supplementary figures are numbered by citation order, tables by chapter order.** Two different principles in one supplement.
- **Filenames do not match compiled numbers.** `Figure1_Lipidomics_Landscape.png` is Figure 3; `Figure5_shiny_triterpenoid_AB.png` is Figure 7.
- **The Methods-before-Results move for Plant Genome is still pending**, along with the unresolved `main.tex` merge conflict.
- **Data availability.** The reviewer flagged this but understated it. The LD-mapping and GEMMA code being outside the repo (see A2) is a bigger problem for a resource paper than the missing Zenodo DOI.
