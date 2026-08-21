#!/usr/bin/env python3
"""
GO_BP over-representation for the LD-based candidate genes (r2 >= 0.4),
computed SEPARATELY for each condition x layer:
  CTL_individual, CTL_sumratio, LIN_individual, LIN_sumratio.

Fisher exact (one-sided greater), Benjamini-Hochberg FDR.
Background = all annotated genes with >=1 GO_BP term.
Term filters: term background >= 5 genes, overlap >= 3 candidate genes.

Run:  python3 go_enrichment_LD_by_condition_layer.py
Requires: scipy
"""
import re, csv, os
from collections import defaultdict
from scipy.stats import fisher_exact

REPO   = "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database"
ANNOT  = f"{REPO}/data/annotation/gene_annotation.txt"
MASTER = f"{REPO}/data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv"
OUT    = f"{REPO}/data/LD_mapped/go_enrichment"
os.makedirs(OUT, exist_ok=True)

def parse_go(s):
    if not s: return []
    out=[]
    for tok in s.split(";"):
        tok=tok.strip()
        if not tok: continue
        m=re.match(r"^(.*)\((GO:[0-9]+)\)$", tok)
        out.append((m.group(1).strip(), m.group(2)) if m else (tok,None))
    return out

# ---- background ----
gene_terms=defaultdict(set)
with open(ANNOT,newline='') as f:
    r=csv.DictReader(f,delimiter='\t'); fld={k.strip():k for k in r.fieldnames}
    for row in r:
        g=row[fld['GeneID']].strip()
        for nm,goid in parse_go(row.get(fld['GO_BP'],'')): gene_terms[g].add((nm,goid))
background=sorted(g for g,t in gene_terms.items() if t); bg=set(background); n_bg=len(background)
term_genes=defaultdict(set); term_id={}
for g in background:
    for nm,goid in gene_terms[g]:
        term_genes[nm].add(g)
        if goid: term_id[nm]=goid

# ---- candidate sets ----
sets=defaultdict(set)
with open(MASTER,newline='') as f:
    r=csv.DictReader(f,delimiter='\t')
    for row in r: sets[(row['condition'],row['layer'])].add(row['GeneID'].strip())

def bh(p):
    m=len(p); order=sorted(range(m),key=lambda i:p[i]); adj=[0]*m; prev=1.0
    for rank,i in enumerate(reversed(order)):
        prev=min(prev,p[i]*m/(m-rank)); adj[i]=min(prev,1.0)
    return adj

def enrich(genes):
    ts=set(genes)&bg; n=len(ts); rows=[]
    for term,tg in term_genes.items():
        ov=sorted(ts&tg); a=len(ov); b=len(tg)
        if b<5 or a<3: continue
        p=fisher_exact([[a,n-a],[b-a,n_bg-b-(n-a)]],alternative='greater')[1]
        exp=n*b/n_bg
        rows.append({'term':term,'go_id':term_id.get(term,''),'bg':b,'obs':a,'exp':exp,
                     'fold':a/exp,'p':p,'genes':';'.join(ov)})
    for r,q in zip(rows,bh([r['p'] for r in rows])): r['q']=q
    rows.sort(key=lambda r:(r['q'],-r['fold']))
    return n,rows

print(f"background genes with GO_BP: {n_bg}\n")
for cond,lay in [("CTL","individual"),("CTL","sumratio"),("LIN","individual"),("LIN","sumratio")]:
    n,rows=enrich(sets[(cond,lay)]); nsig=sum(1 for r in rows if r['q']<0.05)
    print(f"{cond}_{lay}: {len(sets[(cond,lay)])} genes, {n} with GO_BP, {len(rows)} terms tested, {nsig} sig (q<0.05)")
    with open(f"{OUT}/{cond}_{lay}_GO_BP_enrichment_LD.tsv","w") as o:
        o.write("go_term\tgo_id\tterm_bg_count\tterm_test_count\texpected\tfold_enrichment\tp_value\tp_adj_bh\toverlap_genes\n")
        for r in rows:
            o.write(f"{r['term']}\t{r['go_id']}\t{r['bg']}\t{r['obs']}\t{r['exp']:.3f}\t{r['fold']:.3f}\t{r['p']:.3e}\t{r['q']:.3e}\t{r['genes']}\n")
print(f"\noutputs -> {OUT}")
