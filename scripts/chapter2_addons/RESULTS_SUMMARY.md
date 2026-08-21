# Results summary — chapter 2 add-on analyses

## 1. Botanical race and population structure

n with race assignment: CTL 222 / 394, LIN 211 / 363.
n with K.Cluster assignment: CTL 387, LIN 355.
Groups (CTL/LIN): Bicolor 17/16, Caudatum 76/71, Durra 22/23, Guinea 19/17, Kafir 7/7, Mixed 81/77.

**Headline: no lipid trait differs significantly by botanical race in either condition
after BH correction.** Smallest raw p was PG in LIN (p = 0.040, q = 0.47). Effect sizes
were uniformly tiny (epsilon-squared 0 to 0.033, i.e. race explains 0-3% of variance).

Genetic cluster (K = 6) gave the same picture with one striking exception:

| Trait | Condition | H | p | q (BH) | epsilon^2 |
|---|---|---|---|---|---|
| PS | LIN | 43.43 | 3.0e-08 | 4.2e-07 | 0.110 |
| PS | CTL | 14.15 | 0.015 | 0.21 | 0.024 |

All other traits: q > 0.2, epsilon-squared < 0.021.

Interpretation for the chapter:
- Lipid-class variation in the SAP is essentially independent of botanical race and of
  genome-wide genetic structure. This *supports* the GWAS: the signals are unlikely to be
  simple structure artefacts, and the 3 genotype PCs in the model are adequate.
- PS is the exception. It is the one class whose LIN abundance tracks genetic cluster
  (11% of variance). Since PS is also one of the four headline remodeling signatures,
  this needs to be stated explicitly.

Relevant abundance context (mean %TIC, computed here):
PS  = 0.014% (CTL) / 0.037% (LIN)   — a ~2.6x shift, but on a class at ~1/100 of a percent of TIC
LPE = 0.0064% / 0.0114%
LPC = 0.096% / 0.127%
(compare MGDG 34.5% / 32.1%). No zeros in any of these, so detection is reliable — the issue
is prominence, not artefact.

## 2. CTL vs LIN GWAS candidate overlap

Gene universe N = 34,027.

### Gene level
| Layer | n CTL | n LIN | shared | expected | fold | Jaccard | P |
|---|---|---|---|---|---|---|---|
| All layers | 1165 | 4323 | 354 | 148 | 2.39 | 0.069 | 4.7e-59 |
| Individual lipids | 1100 | 4323 | 289 | 140 | 2.07 | 0.056 | 4.3e-35 |
| Class sums / ratios | 115 | 812 | 18 | 2.7 | 6.56 | 0.020 | 2.9e-10 |

### Locus level (genome tiled into windows containing >= 1 gene)
| Window | Layer | n CTL | n LIN | shared | expected | fold | P |
|---|---|---|---|---|---|---|---|
| 100 kb | Individual lipids | 521 | 1564 | 193 | 168 | 1.15 | 0.008 |
| 250 kb | Individual lipids | 378 | 971 | 168 | 161 | **1.04** | **0.23 (n.s.)** |
| 500 kb | Individual lipids | 297 | 656 | 162 | 157 | **1.03** | **0.29 (n.s.)** |
| 100 kb | Class sums / ratios | 50 | 247 | 7 | 2.5 | 2.75 | 0.013 |
| 250 kb | Class sums / ratios | 36 | 161 | 6 | 2.5 | 2.36 | 0.038 |
| 500 kb | Class sums / ratios | 28 | 121 | 5 | 2.7 | 1.83 | 0.13 |

**The gene-level enrichment is an LD artefact.** Once genes are collapsed into independent
loci, the individual-lipid overlap is indistinguishable from chance. Class-level composite
traits (sums / ratios) replicate somewhat better than individual species, which is the one
positive reproducibility result here.

### Gene-count inflation (250 kb loci)
| Condition | Layer | genes | loci | genes/locus |
|---|---|---|---|---|
| CTL | individual lipids | 1100 | 378 | 2.9 |
| CTL | class sums/ratios | 115 | 36 | 3.2 |
| LIN | individual lipids | 4323 | 971 | 4.5 |
| LIN | class sums/ratios | 812 | 161 | 5.0 |

The "4,323 candidate genes" in LIN correspond to ~971 independent loci.

### Per-class (individual-lipid GWAS, gene level)
| Class | n CTL | n LIN | shared | fold | q |
|---|---|---|---|---|---|
| TG | 866 | 890 | 37 | 1.63 | 0.009 |
| PC | 58 | 1324 | 8 | 3.55 | 0.009 |
| Fatty acid | 109 | 1897 | 7 | 1.15 | 0.95 (n.s.) |
| Terpenoid | 37 | 39 | **0** | - | 1 |
| DG | 13 | 34 | **0** | - | 1 |
| SQDG | 9 | 396 | **0** | - | 1 |
| AEG | 7 | 119 | **0** | - | 1 |

SQDG — one of the four headline remodeling signatures — shares **zero** candidate genes
between the two trials.

### Class concordance of the shared genes
Of the 354 genes shared across all layers, only 90 (25.4%) are shared *for the same lipid
class*. For the individual-lipid layer alone, 52 / 289 (18.0%).

### The genuinely recurrent loci
Only **two** independent loci recur between CTL and LIN in the class-sum/ratio GWAS:
- chr10 ~ SORBI_3010G150325-SORBI_3010G151325 (10 gene models in one LD block; includes
  CMP-sialic acid transporter 1 = SORBI_3010G151000, the gene the chapter already discusses)
- chr1 ~ SORBI_3001G168800-SORBI_3001G171600 (8 gene models in one LD block)

These are single loci, not 18 genes. Full lists in `table/overlap/shared_genes_*.csv`.

## What this means for the chapter text
1. The Conclusion's claim of environment-dependent architecture is *strengthened* — quantify it.
2. The candidate-gene counts (1,100 / 4,323) must be reported alongside locus counts.
3. Add the race/structure result as a positive control for the GWAS.
4. The SQDG signature has zero genetic support shared across trials; say so.
5. PS: flag the population-structure association and the very low absolute abundance.
