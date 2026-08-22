import sys,os,csv
sys.path.insert(0,os.path.expanduser('~/tmp/soldx')); csv.field_size_limit(10**8)
import pandas as pd, numpy as np
from stats_lite import overlap_stats
REPO=os.path.expanduser('~/mnt/SAP-Lipidomics-Database')
TAB=os.path.join(REPO,'table','overlap'); os.makedirs(TAB,exist_ok=True)

# chromosome lengths from the gene range file (max gene end per chr)
g=pd.read_csv(os.path.join(REPO,'data/LD_mapped/genes_ranges/genes.range'),sep='\t',header=None,
              names=['chr','start','end','gene'])
chrlen=g.groupby('chr')['end'].max().to_dict()
print('chrom lengths (Mb):',{k:round(v/1e6,1) for k,v in sorted(chrlen.items())},
      '| genome',round(sum(chrlen.values())/1e6,1),'Mb')

rows=list(csv.DictReader(open(os.path.join(REPO,'data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv')),delimiter='\t'))
d=pd.DataFrame(rows)
d['Chromosome']=d.Chromosome.astype(int); d['Gene_Start']=d.Gene_Start.astype(int)

out=[]
print()
for W in [100_000,250_000,500_000]:
    # genome-wide window universe = windows that contain >=1 annotated gene (testable space)
    g['win']=g['chr'].astype(str)+'_'+(g['start']//W).astype(str)
    NW=g['win'].nunique()
    d['win']=d.Chromosome.astype(str)+'_'+(d.Gene_Start//W).astype(str)
    print(f'--- window = {W//1000} kb : {NW} gene-containing windows genome-wide ---')
    for lab,sel in [('All layers',d),('Individual lipids',d[d.layer=='individual']),
                    ('Class sums / ratios',d[d.layer=='sumratio'])]:
        A=set(sel.loc[sel.condition=='CTL','win']); B=set(sel.loc[sel.condition=='LIN','win'])
        s=overlap_stats(A,B,NW); s.update(window_kb=W//1000,Layer=lab,universe=NW); out.append(s)
        print(f'  {lab:22s} CTL={len(A):4d} LIN={len(B):4d} shared={len(A&B):4d} '
              f'exp={s["expected"]:6.1f} fold={s["fold"]:5.2f} J={s["jaccard"]:.3f} P={s["p"]:.2e}')
    print()
r=pd.DataFrame(out)[['window_kb','Layer','universe','n_A','n_B','n_shared','expected','fold','jaccard','p']]
r.columns=['window_kb','Layer','n_windows_genome','n_CTL','n_LIN','n_shared','expected_shared','fold_enrichment','jaccard','p_hypergeom']
r.to_csv(os.path.join(TAB,'gwas_overlap_locus_level.csv'),index=False)

# inflation table
W=250_000
d['win']=d.Chromosome.astype(str)+'_'+(d.Gene_Start//W).astype(str)
inf=[]
print('=== Gene-count inflation (250 kb windows) ===')
for c in ['CTL','LIN']:
    for lay,lab in [('individual','individual lipids'),('sumratio','class sums/ratios')]:
        s=d[(d.condition==c)&(d.layer==lay)]
        inf.append(dict(Condition=c,Layer=lab,n_genes=s.GeneID.nunique(),n_loci=s.win.nunique(),
                        genes_per_locus=s.GeneID.nunique()/max(s.win.nunique(),1)))
        print(f'  {c} {lab:20s}: {s.GeneID.nunique():5d} genes -> {s.win.nunique():4d} loci '
              f'({inf[-1]["genes_per_locus"]:.1f} genes/locus)')
    s=d[d.condition==c]
    inf.append(dict(Condition=c,Layer='all layers',n_genes=s.GeneID.nunique(),n_loci=s.win.nunique(),
                    genes_per_locus=s.GeneID.nunique()/max(s.win.nunique(),1)))
    print(f'  {c} {"all layers":20s}: {s.GeneID.nunique():5d} genes -> {s.win.nunique():4d} loci '
          f'({inf[-1]["genes_per_locus"]:.1f} genes/locus)')
pd.DataFrame(inf).to_csv(os.path.join(TAB,'gwas_gene_to_locus_inflation.csv'),index=False)
