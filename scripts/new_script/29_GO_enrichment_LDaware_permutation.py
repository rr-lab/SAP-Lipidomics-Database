#!/usr/bin/env python3
"""LD-aware permutation validation of class-stratified GO-BP enrichment.

Problem
-------
LD-based candidate mapping assigns every gene model in a linkage block to one
association, so genes in a block - especially tandem arrays - are not
independent observations. A gene-level hypergeometric test therefore treats one
locus as several draws and inflates significance.

Approach
--------
The observed statistic stays at gene level (number of candidate genes annotated
with a term), which preserves the resolution of the full GO background. Only the
NULL changes: instead of drawing independent genes, each permutation redraws the
same number of LINKAGE BLOCKS, each contributing the same number of genes as the
corresponding observed block, sampled from background blocks matched on
gene density. This asks the correct question - is the term overrepresented given
that we sampled this many LD blocks of this size from gene regions of this
density? - and a term carried by a single tandem array no longer survives.

Only terms already BH-significant at gene level are tested, so this is a
validation filter on published results rather than a new screen.
"""

import csv, os, re, numpy as np
from collections import defaultdict

REPO   = os.path.expanduser("~/mnt/SAP-Lipidomics-Database")
ANNOT  = f"{REPO}/data/annotation/gene_annotation.txt"
RANGES = f"{REPO}/data/LD_mapped/genes_ranges/genes.range"
MASTER = f"{REPO}/data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv"
GO     = f"{REPO}/data/LD_mapped/go_enrichment"
WIN, N_PERM, SEED = 250_000, 8000, 1
csv.field_size_limit(10**8)
rng = np.random.default_rng(SEED)

def parse_go_bp(v):
    out = []
    for t in (v or "").split(";"):
        t = t.strip()
        if t:
            m = re.match(r"^(.*)\((GO:[0-9]+)\)$", t)
            out.append(m.group(1).strip() if m else t)
    return out

gene_window, window_genes = {}, defaultdict(list)
for line in open(RANGES):
    p = line.split()
    if len(p) >= 4:
        w = (int(p[0]), int(p[1]) // WIN); gene_window[p[3]] = w; window_genes[w].append(p[3])

gene_terms = defaultdict(set)
for r in csv.DictReader(open(ANNOT), delimiter="\t"):
    gene_terms[r["GeneID"].strip()].update(parse_go_bp(r.get("GO_BP", "")))
annot = sorted(g for g, t in gene_terms.items() if t and g in gene_window)
idx = {g: i for i, g in enumerate(annot)}
NG = len(annot)

win_annot = defaultdict(list)
for g in annot: win_annot[gene_window[g]].append(idx[g])
bg_windows = [w for w in win_annot if win_annot[w]]
dens = np.array([len(win_annot[w]) for w in bg_windows])
qs = np.quantile(dens, [.2, .4, .6, .8])
strat = np.searchsorted(qs, dens, side="right")
by_strat = {s: [bg_windows[i] for i in np.where(strat == s)[0]] for s in range(5)}
w_strat = {w: int(strat[i]) for i, w in enumerate(bg_windows)}

term_vec = {}
for g in annot:
    for t in gene_terms[g]:
        term_vec.setdefault(t, np.zeros(NG, bool))[idx[g]] = True

layers = {"individual": lambda p: ([m.group(1)] if (m := re.match(r"^([A-Za-z0-9]+)\(", p.strip())) else []),
          "sumratio":   lambda p: [v for v in p.strip().replace("Sum_", "", 1).replace("_log10safe", "").split("_over_") if v]}
cand = {L: defaultdict(lambda: defaultdict(set)) for L in layers}
for r in csv.DictReader(open(MASTER), delimiter="\t"):
    L = r["layer"].strip().lower()
    if L in layers and r["GeneID"].strip() in idx:
        for ph in filter(None, map(str.strip, r.get("Phenotypes", "").split(";"))):
            for cl in layers[L](ph): cand[L][r["condition"].strip()][cl].add(r["GeneID"].strip())

def null_matrix(genes):
    """B x NG boolean matrix; each row redraws the same LD-block structure."""
    blocks = defaultdict(list)
    for g in genes: blocks[gene_window[g]].append(g)
    sizes_by_stratum = defaultdict(list)
    for w, gs in blocks.items(): sizes_by_stratum[w_strat[w]].append(len(gs))
    M = np.zeros((N_PERM, NG), bool)
    for s, sizes in sizes_by_stratum.items():
        pool = by_strat[s]
        picks = rng.integers(0, len(pool), size=(N_PERM, len(sizes)))
        for j, need in enumerate(sizes):
            for b in range(N_PERM):
                avail = win_annot[pool[picks[b, j]]]
                take = avail if len(avail) <= need else [avail[i] for i in rng.choice(len(avail), need, replace=False)]
                M[b, take] = True
    return M

print(f"annotated genes {NG} | background LD blocks {len(bg_windows)} | permutations {N_PERM}\n")
for L in ("individual", "sumratio"):
    src = f"{GO}/GO_BP_enrichment_by_lipid_class_{L}_LD.tsv"
    rows = [r for r in csv.DictReader(open(src), delimiter="\t") if float(r["p_adj_bh"]) < 0.05]
    out = []
    for (cond, cl), grp in {(c, k): [r for r in rows if r["condition"] == c and r["lipid_class"] == k]
                            for c in {r["condition"] for r in rows} for k in {r["lipid_class"] for r in rows}}.items():
        if not grp: continue
        genes = cand[L][cond][cl]
        M = null_matrix(genes)
        nblk = len({gene_window[g] for g in genes})
        for r in grp:
            obs = int(r["term_test_count"])
            null = M[:, term_vec[r["go_term"]]].sum(axis=1)   # bool matmul saturates; index+sum counts
            p = (int((null >= obs).sum()) + 1) / (N_PERM + 1)
            r.update(n_ld_blocks_in_set=nblk,
                     blocks_carrying_term=len({gene_window[g] for g in r["overlap_genes"].split(";") if g in gene_window}),
                     null_mean=round(float(null.mean()), 2), p_perm=round(p, 5))
            out.append(r)
    for r, q in zip(out, [None]*len(out)):
        pass
    ps = sorted(range(len(out)), key=lambda i: out[i]["p_perm"])
    run = 1.0; qv = [0.0]*len(out)
    for rr, i in enumerate(reversed(ps), 1):
        rank = len(out) - rr + 1
        run = min(run, out[i]["p_perm"] * len(out) / rank); qv[i] = min(run, 1.0)
    for r, q in zip(out, qv): r["q_perm_bh"] = round(q, 5)
    out.sort(key=lambda r: r["p_perm"])
    path = f"{GO}/GO_BP_enrichment_by_lipid_class_{L}_LDaware.tsv"
    fields = ["condition","lipid_class","go_term","go_id","term_test_count","n_ld_blocks_in_set",
              "blocks_carrying_term","fold_enrichment","expected","null_mean","p_value","p_adj_bh",
              "p_perm","q_perm_bh","overlap_genes"]
    with open(path,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",extrasaction="ignore"); w.writeheader(); w.writerows(out)
    keep=[r for r in out if r["q_perm_bh"]<0.05]
    print(f"=== {L}: {len(out)} gene-level BH-significant terms -> {len(keep)} survive the LD-aware null ===")
    for r in out:
        mark = "KEEP  " if r["q_perm_bh"]<0.05 else "DROP  "
        print(f"  {mark}{r['condition']:3} {r['lipid_class']:5} {r['go_term'][:40]:40} "
              f"genes={r['term_test_count']:>2} blocks={r['blocks_carrying_term']:>2} "
              f"obs/null={r['term_test_count']}/{r['null_mean']:<5} p_perm={r['p_perm']:.4f} q={r['q_perm_bh']:.3f}")
    print(f"  -> {path}\n")
