#!/usr/bin/env python3
"""Locus-collapsed GO-BP enrichment for LD-mapped GWAS candidates.

Motivation
----------
LD-based candidate mapping assigns every gene model within a linkage block to a
single association. Genes in a block - especially tandem arrays - share GO terms
by descent, so a gene-level hypergeometric test treats one locus as several
independent observations and inflates significance.

This script repeats the class-stratified GO-BP enrichment with the independent
locus, not the gene, as the unit of observation:

  * the genome is tiled into non-overlapping WIN-bp windows;
  * a window is a *background* window if it contains >=1 GO-BP-annotated gene;
  * a window is a *term* window if it contains >=1 gene annotated with that term;
  * a window is a *candidate* window if it contains >=1 candidate gene;
  * enrichment is a one-sided hypergeometric test on window counts.

Terms passing BH are additionally checked with a gene-density-matched
permutation test, because candidate windows are drawn preferentially from
gene-rich regions (more genes -> more SNPs -> more chance of a hit).

Outputs both trait layers and a gene-level vs locus-level comparison.
"""

import csv, math, os, re, random
from collections import defaultdict

REPO   = os.path.expanduser("~/mnt/SAP-Lipidomics-Database")
ANNOT  = f"{REPO}/data/annotation/gene_annotation.txt"      # canonical GO source
RANGES = f"{REPO}/data/LD_mapped/genes_ranges/genes.range"
MASTER = f"{REPO}/data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv"
OUT    = f"{REPO}/data/LD_mapped/go_enrichment"
WIN        = 250_000    # locus width (bp)
MIN_TERM   = 5          # min background windows carrying the term
MIN_LOCI   = 3          # min candidate loci overlapping the term
N_PERM     = 10_000
SEED       = 1
os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10**8)


def parse_go_bp(value):
    out = []
    for token in (value or "").split(";"):
        token = token.strip()
        if not token:
            continue
        m = re.match(r"^(.*)\((GO:[0-9]+)\)$", token)
        out.append((m.group(1).strip(), m.group(2)) if m else (token, ""))
    return out


def bh_adjust(pvalues):
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adjusted = [0.0] * len(pvalues)
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = len(pvalues) - reverse_rank + 1
        running = min(running, pvalues[index] * len(pvalues) / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def hyper_sf(k, N, K, n):
    """P(X >= k), X ~ Hypergeometric(N population, K successes, n draws)."""
    lo, hi = max(0, n + K - N), min(K, n)
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0
    denom = math.comb(N, n)
    return min(1.0, sum(math.comb(K, i) * math.comb(N - K, n - i) / denom
                        for i in range(int(k), int(hi) + 1)))


def individual_class(ph):
    m = re.match(r"^([A-Za-z0-9]+)\(", ph.strip())
    return [m.group(1)] if m else []


def sumratio_classes(ph):
    x = ph.strip().replace("Sum_", "", 1).replace("_log10safe", "")
    return [v for v in x.split("_over_") if v]


# ---------------------------------------------------------------- gene -> window
gene_window, window_genes = {}, defaultdict(set)
for line in open(RANGES):
    p = line.split()
    if len(p) >= 4:
        w = (int(p[0]), int(p[1]) // WIN)
        gene_window[p[3]] = w
        window_genes[w].add(p[3])

# ---------------------------------------------------------------- GO annotation
gene_terms = defaultdict(set)
for row in csv.DictReader(open(ANNOT), delimiter="\t"):
    gene_terms[row["GeneID"].strip()].update(parse_go_bp(row.get("GO_BP", "")))
annotated = {g for g, t in gene_terms.items() if t and g in gene_window}

background = {gene_window[g] for g in annotated}
N_BG = len(background)

term_windows, term_ids = defaultdict(set), {}
for g in annotated:
    for term, go_id in gene_terms[g]:
        term_windows[term].add(gene_window[g])
        if go_id:
            term_ids[term] = go_id

# gene-density stratum for each background window (quintiles of annotated genes)
dens = {w: sum(1 for g in window_genes[w] if g in annotated) for w in background}
cuts = sorted(dens.values())
q = [cuts[int(len(cuts) * f)] for f in (0.2, 0.4, 0.6, 0.8)]
def stratum(w):
    d = dens[w]
    return sum(d > c for c in q)
strata = defaultdict(list)
for w in background:
    strata[stratum(w)].append(w)

# ---------------------------------------------------------------- candidate sets
layers = {"individual": individual_class, "sumratio": sumratio_classes}
cand = {L: defaultdict(lambda: defaultdict(set)) for L in layers}
cand_genes = {L: defaultdict(lambda: defaultdict(set)) for L in layers}
for row in csv.DictReader(open(MASTER), delimiter="\t"):
    L = row["layer"].strip().lower()
    if L not in layers:
        continue
    cond, gene = row["condition"].strip(), row["GeneID"].strip()
    if gene not in gene_window:
        continue
    for ph in filter(None, map(str.strip, row.get("Phenotypes", "").split(";"))):
        for cl in layers[L](ph):
            cand[L][cond][cl].add(gene_window[gene])
            cand_genes[L][cond][cl].add(gene)

random.seed(SEED)
print(f"Locus width {WIN//1000} kb | background windows with >=1 GO-annotated gene: {N_BG}")

for L in ("individual", "sumratio"):
    rows = []
    for cond in sorted(cand[L]):
        for cl in sorted(cand[L][cond]):
            cwin = cand[L][cond][cl] & background
            n = len(cwin)
            if n < MIN_LOCI:
                continue
            for term, twin in term_windows.items():
                K = len(twin)
                hits = cwin & twin
                k = len(hits)
                if K < MIN_TERM or k < MIN_LOCI:
                    continue
                exp = n * K / N_BG
                rows.append(dict(
                    condition=cond, lipid_class=cl, go_term=term,
                    go_id=term_ids.get(term, ""), term_bg_windows=K,
                    candidate_loci=n, term_loci=k, expected_loci=round(exp, 4),
                    fold_enrichment=round(k / exp, 3),
                    p_value=hyper_sf(k, N_BG, K, n),
                    n_genes=len({g for g in cand_genes[L][cond][cl]
                                 if g in annotated and gene_window[g] in hits}),
                    loci_genes=";".join(sorted(
                        g for g in cand_genes[L][cond][cl]
                        if g in annotated and gene_window[g] in hits))))
    for r, qv in zip(rows, bh_adjust([r["p_value"] for r in rows])):
        r["p_adj_bh"] = qv

    # gene-density-matched permutation for BH-significant terms
    for r in rows:
        if r["p_adj_bh"] >= 0.05:
            r["p_perm"] = ""
            continue
        n = r["candidate_loci"]
        cwin = cand[L][r["condition"]][r["lipid_class"]] & background
        counts = defaultdict(int)
        for w in cwin:
            counts[stratum(w)] += 1
        twin = term_windows[r["go_term"]]
        ge = 0
        for _ in range(N_PERM):
            hit = 0
            for s, cnt in counts.items():
                pool = strata[s]
                hit += sum(1 for w in random.sample(pool, min(cnt, len(pool))) if w in twin)
            ge += hit >= r["term_loci"]
        r["p_perm"] = (ge + 1) / (N_PERM + 1)

    rows.sort(key=lambda r: (r["p_adj_bh"], r["p_value"], -r["fold_enrichment"]))
    fields = ["condition", "lipid_class", "go_term", "go_id", "term_bg_windows",
              "candidate_loci", "term_loci", "n_genes", "expected_loci",
              "fold_enrichment", "p_value", "p_adj_bh", "p_perm", "loci_genes"]
    path = f"{OUT}/GO_BP_enrichment_by_lipid_class_{L}_LOCUS.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    sig = [r for r in rows if r["p_adj_bh"] < 0.05]
    perm = [r for r in sig if r["p_perm"] != "" and r["p_perm"] < 0.05]
    print(f"\n=== {L}: {len(rows)} locus-level tests | {len(sig)} BH-significant "
          f"| {len(perm)} also permutation-significant ===")
    for r in sig:
        flag = "" if (r["p_perm"] != "" and r["p_perm"] < 0.05) else "   <- fails permutation"
        print(f"  {r['condition']:3} {r['lipid_class']:5} {r['go_term'][:42]:42} "
              f"loci={r['term_loci']:2}/{r['candidate_loci']:<4} fold={r['fold_enrichment']:5.1f} "
              f"q={r['p_adj_bh']:.2e} perm={r['p_perm']}{flag}")
    print(f"  -> {path}")
