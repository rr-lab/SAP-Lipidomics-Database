#!/usr/bin/env python3
"""
GO_BP over-representation (Fisher exact, one-sided greater, BH-corrected) for the
CTL and LIN GWAS candidate-gene sets, combining individual + sum/ratio trait layers.

Reproduces run_combined_condition_go_enrichment.R and extends it with:
  (1) recurrence-filtered subsets (top ~10% / ~5% / ~1% by N_Phenotypes), and
  (2) a locus-clustering breakdown of every significant term, to flag hits that
      collapse onto a single genomic region (LD / physical-clustering artifact).

Background = all annotated genes carrying >=1 GO_BP term.
Thresholds match the R script: term_bg >= 5, term_overlap >= 3.
Recurrence per gene per condition = max N_Phenotypes across the two trait layers.
Loci: gene midpoints on the same chromosome are merged when within 2 x 25 kb
      (i.e. their +/-25 kb windows overlap), matching the manuscript's window logic.

Run:  python3 go_enrichment_recurrence.py
Requires: scipy, numpy
"""
import re, csv, os
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact

# ---- paths (edit if your checkout differs) --------------------------------
REPO   = "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database"
THESIS = "/Users/nirwantandukar/Documents/Github/Thesis_NirwanTandukar_GeneticsGenomics/ncsuthesis-0.6/Supplementary_tables/Chapter2"
ANNOT  = f"{REPO}/data/annotation/gene_annotation.txt"
OUT    = f"{REPO}/curated_gwas_themes/go_enrichment"
COND_FILES = {
 "CTL": [f"{THESIS}/SuppTable_S7_CTL_GWAS_candidate_genes_for_individual_lipid_traits.tsv",
         f"{THESIS}/SuppTable_S8_CTL_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv"],
 "LIN": [f"{THESIS}/SuppTable_S9_LIN_GWAS_candidate_genes_for_individual_lipid_traits.tsv",
         f"{THESIS}/SuppTable_S10_LIN_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv"],
}
CUTS = [("full", None), ("top10", 0.90), ("top5", 0.95), ("top1", 0.99)]
MIN_BG, MIN_OVERLAP, WINDOW = 5, 3, 25000
os.makedirs(OUT, exist_ok=True)

def parse_go(s):
    if not s: return []
    out = []
    for tok in s.split(";"):
        tok = tok.strip()
        if not tok: continue
        m = re.match(r"^(.*)\((GO:[0-9]+)\)$", tok)
        out.append((m.group(1).strip(), m.group(2)) if m else (tok, None))
    return out

# ---- background from annotation -------------------------------------------
gene_terms = defaultdict(set)
with open(ANNOT, newline='') as f:
    r = csv.DictReader(f, delimiter='\t'); fld = {k.strip(): k for k in r.fieldnames}
    for row in r:
        g = row[fld['GeneID']].strip()
        for nm, goid in parse_go(row.get(fld['GO_BP'], '')): gene_terms[g].add((nm, goid))
background = sorted(g for g, t in gene_terms.items() if t)
bg_set = set(background); n_bg = len(background)
term_genes = defaultdict(set); term_id = {}
for g in background:
    for nm, goid in gene_terms[g]:
        term_genes[nm].add(g)
        if goid: term_id[nm] = goid

# ---- per-condition recurrence + gene positions ----------------------------
cond_rec, cond_pos = {}, {}
for cond, files in COND_FILES.items():
    rec = defaultdict(int); pos = {}
    for fn in files:
        with open(fn, newline='') as f:
            r = csv.DictReader(f, delimiter='\t'); fld = {k.strip(): k for k in r.fieldnames}
            for row in r:
                g = row[fld['GeneID']].strip()
                try: n = int(float(row[fld['N_Phenotypes']]))
                except: n = 0
                if n > rec[g]: rec[g] = n
                try:
                    c = row[fld['Chromosome']].strip()
                    s = int(float(row[fld['Gene_Start']])); e = int(float(row[fld['Gene_End']]))
                    pos[g] = (c, (s + e) // 2)
                except: pass
    cond_rec[cond] = rec; cond_pos[cond] = pos

def bh(p):
    m = len(p); order = sorted(range(m), key=lambda i: p[i]); adj = [0]*m; prev = 1.0
    for rank, i in enumerate(reversed(order)):
        prev = min(prev, p[i]*m/(m-rank)); adj[i] = min(prev, 1.0)
    return adj

def enrich(genes):
    ts = set(genes) & bg_set; n = len(ts); rows = []
    for term, tg in term_genes.items():
        ov = ts & tg; a = len(ov); b = len(tg)
        if b < MIN_BG or a < MIN_OVERLAP: continue
        p = fisher_exact([[a, n-a], [b-a, n_bg-b-(n-a)]], alternative='greater')[1]
        exp = n*b/n_bg
        rows.append({'term': term, 'go_id': term_id.get(term, ''), 'bg': b, 'obs': a,
                     'exp': exp, 'fold': a/exp, 'p': p, 'genes': sorted(ov)})
    for r, q in zip(rows, bh([r['p'] for r in rows])): r['q'] = q
    rows.sort(key=lambda r: (r['q'], -r['fold']))
    return n, rows

def cluster(genes, pos):
    pts = defaultdict(list)
    for g in genes:
        if g in pos: pts[pos[g][0]].append((pos[g][1], g))
    loci = []
    for c, lst in pts.items():
        lst.sort(); cur = [lst[0]]
        for p, g in lst[1:]:
            if p - cur[-1][0] <= 2*WINDOW: cur.append((p, g))
            else: loci.append((c, cur)); cur = [(p, g)]
        loci.append((c, cur))
    npos = sum(len(v) for v in pts.values())
    return loci, npos

def wtsv(path, rows, cut=None):
    hdr = ['go_term','go_id','term_bg_count','term_test_count','expected','fold_enrichment','p_value','p_adj_bh','overlap_genes']
    if cut is not None: hdr = ['recurrence_cut'] + hdr
    with open(path, 'w') as o:
        o.write('\t'.join(hdr) + '\n')
        for r in rows:
            base = [r['term'], r['go_id'], r['bg'], r['obs'], f"{r['exp']:.3f}", f"{r['fold']:.3f}",
                    f"{r['p']:.3e}", f"{r['q']:.3e}", ';'.join(r['genes'])]
            o.write('\t'.join(map(str, ([cut]+base) if cut is not None else base)) + '\n')

locus_rows = []
for cond in ('CTL', 'LIN'):
    rec, pos = cond_rec[cond], cond_pos[cond]; vals = sorted(rec.values()); genes_all = list(rec)
    by_rec = []
    for label, qthr in CUTS:
        if qthr is None: sub, cutval = genes_all, 'all'
        else: cutval = int(np.quantile(vals, qthr)); sub = [g for g in genes_all if rec[g] >= cutval]
        n, rows = enrich(sub)
        if label == 'full':
            wtsv(f"{OUT}/{cond}_GO_BP_enrichment_full.tsv", rows)
        else:
            for r in rows: by_rec.append((f"{label}(rec>={cutval})", r))
        for r in [x for x in rows if x['q'] < 0.05]:
            loci, npos = cluster(r['genes'], pos)
            big = max((len(l[1]) for l in loci), default=0)
            desc = "; ".join(f"chr{c}:{min(m[0] for m in mm)}-{max(m[0] for m in mm)}({len(mm)}g)" for c, mm in loci)
            verdict = "SINGLE LOCUS" if len(loci) == 1 else ("mostly 1 locus" if big/max(npos,1) >= 0.6 else "multi-locus")
            locus_rows.append([cond, f"{label}(rec>={cutval})", r['term'], r['go_id'], r['obs'],
                               round(r['fold'],2), f"{r['q']:.2e}", npos, len(loci), big, verdict, desc, ';'.join(r['genes'])])
        print(f"{cond} {label} (rec>={cutval}): {len(sub)} genes, {n} w/GO, {len(rows)} tested, {sum(1 for x in rows if x['q']<0.05)} sig(q<0.05)")
    with open(f"{OUT}/{cond}_GO_BP_enrichment_by_recurrence.tsv", 'w') as o:
        o.write('recurrence_cut\tgo_term\tgo_id\tterm_bg_count\tterm_test_count\texpected\tfold_enrichment\tp_value\tp_adj_bh\toverlap_genes\n')
        for lab, r in by_rec:
            o.write('\t'.join(map(str, [lab, r['term'], r['go_id'], r['bg'], r['obs'], f"{r['exp']:.3f}",
                    f"{r['fold']:.3f}", f"{r['p']:.3e}", f"{r['q']:.3e}", ';'.join(r['genes'])])) + '\n')

with open(f"{OUT}/significant_terms_locus_breakdown.tsv", 'w') as o:
    o.write('condition\trecurrence_cut\tgo_term\tgo_id\tobs\tfold\tq\tn_genes_positioned\tn_independent_loci\tmax_genes_in_one_locus\tverdict\tloci\tgenes\n')
    for row in locus_rows: o.write('\t'.join(map(str, row)) + '\n')

print(f"\nBackground genes with GO_BP: {n_bg}")
print(f"Outputs written to: {OUT}")
