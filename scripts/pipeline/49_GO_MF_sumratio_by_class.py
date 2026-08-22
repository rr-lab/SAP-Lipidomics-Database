#!/usr/bin/env python3
"""Class-stratified GO-MF enrichment for LD-mapped class-sum and ratio GWAS genes.
A ratio gene is assigned to both numerator and denominator lipid classes."""
import csv, math, re
from collections import defaultdict
from pathlib import Path
repo=Path(__file__).resolve().parents[2]
annot=repo/'data/annotation/gene_annotation.txt'
master=repo/'data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv'
out=repo/'data/LD_mapped/go_enrichment/GO_MF_enrichment_by_lipid_class_sumratio_LD.tsv'
summary=repo/'data/LD_mapped/go_enrichment/GO_MF_enrichment_by_lipid_class_sumratio_LD_summary.tsv'
def parse_terms(v):
    z=[]
    for x in (v or '').split(';'):
        x=x.strip()
        if x:
            m=re.match(r'^(.*)\((GO:\d+)\)$',x)
            z.append((m.group(1).strip(),m.group(2)) if m else (x,''))
    return z
def classes(ph):
    x=ph.strip().replace('Sum_','',1).replace('_log10safe','')
    return [v for v in x.split('_over_') if v]
def bh(p):
    order=sorted(range(len(p)),key=lambda i:p[i]); ans=[0.]*len(p); run=1.
    for rr,i in enumerate(reversed(order),1):
        rank=len(p)-rr+1; run=min(run,p[i]*len(p)/rank); ans[i]=min(run,1.)
    return ans
def fish(a,n,term,bg):
    den=math.comb(bg,n); return sum(math.comb(term,k)*math.comb(bg-term,n-k)/den for k in range(a,min(n,term)+1))
gene_terms=defaultdict(set)
with open(annot) as fh:
    for r in csv.DictReader(fh,delimiter='\t'): gene_terms[r['GeneID'].strip()].update(parse_terms(r.get('GO_MF','')))
bg={g for g,t in gene_terms.items() if t}; termgenes=defaultdict(set); termids={}
for g in bg:
    for term,go in gene_terms[g]: termgenes[term].add(g); termids[term]=go or termids.get(term,'')
classgenes=defaultdict(lambda:defaultdict(set)); classph=defaultdict(lambda:defaultdict(set))
with open(master) as fh:
    for r in csv.DictReader(fh,delimiter='\t'):
        if r['layer'].strip()!='sumratio': continue
        for ph in filter(None,map(str.strip,r.get('Phenotypes','').split(';'))):
            for cl in classes(ph): classgenes[r['condition'].strip()][cl].add(r['GeneID'].strip()); classph[r['condition'].strip()][cl].add(ph)
rows=[]; summ=[]
for cond in sorted(classgenes):
  for cl,genes in sorted(classgenes[cond].items()):
    test=genes&bg; summ.append({'condition':cond,'lipid_class':cl,'sumratio_phenotypes':len(classph[cond][cl]),'candidate_genes':len(genes),'GO_MF_testable_genes':len(test)})
    for term,tset in termgenes.items():
      ov=sorted(test&tset)
      if len(tset)<5 or len(ov)<3: continue
      exp=len(test)*len(tset)/len(bg); pv=fish(len(ov),len(test),len(tset),len(bg))
      rows.append({'condition':cond,'lipid_class':cl,'go_term':term,'go_id':termids.get(term,''),'term_bg_count':len(tset),'term_test_count':len(ov),'expected':exp,'fold_enrichment':len(ov)/exp,'p_value':pv,'overlap_genes':';'.join(ov)})
adj=bh([x['p_value'] for x in rows])
for r,q in zip(rows,adj): r['p_adj_bh']=q
rows.sort(key=lambda r:(r['p_adj_bh'],r['p_value'],-r['fold_enrichment']))
fields=['condition','lipid_class','go_term','go_id','term_bg_count','term_test_count','expected','fold_enrichment','p_value','p_adj_bh','overlap_genes']
with open(out,'w',newline='') as fh:
 w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
for r in summ:
 x=[v for v in rows if v['condition']==r['condition'] and v['lipid_class']==r['lipid_class']];r['GO_terms_tested']=len(x);r['BH_significant_terms']=sum(v['p_adj_bh']<.05 for v in x)
with open(summary,'w',newline='') as fh:
 w=csv.DictWriter(fh,fieldnames=list(summ[0]),delimiter='\t');w.writeheader();w.writerows(summ)
print('tests',len(rows),'nominal',sum(r['p_value']<.05 for r in rows),'BH',sum(r['p_adj_bh']<.05 for r in rows))
print(out)
