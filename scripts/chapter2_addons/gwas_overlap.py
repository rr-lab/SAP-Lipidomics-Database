import sys, os, re, csv, collections
sys.path.insert(0, os.path.expanduser('~/tmp/soldx'))
csv.field_size_limit(10**8)
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from stats_lite import overlap_stats, bh

REPO = os.path.expanduser('~/mnt/SAP-Lipidomics-Database')
FIG = os.path.join(REPO,'fig','overlap'); os.makedirs(FIG, exist_ok=True)
TAB = os.path.join(REPO,'table','overlap'); os.makedirs(TAB, exist_ok=True)

UNIV = {p.split()[3] for p in open(os.path.join(REPO,'data/LD_mapped/genes_ranges/genes.range')) if len(p.split())>=4}
N = len(UNIV)

ABBR = ['MGDG','DGDG','SQDG','GalCer','DGTS','AEG','LPC','LPE','Cer','SM','DG','TG','MG',
        'PC','PE','PG','PA','PS','FA']
FALLBACK = {'Fatty acid and conjugates':'FA','Fatty acid amide':'FA','Fatty amide':'FA',
            'Fatty acid ester':'FA','Oxidised fatty acid':'FA',
            'Triterpenoid':'Terpenoid','Tetraterpenoid':'Terpenoid','Diterpenoid':'Terpenoid',
            'Monoterpenoid':'Terpenoid','Sesquiterpenoid':'Terpenoid','Prenol':'Terpenoid',
            'Cyclic monoterpenoid':'Terpenoid','Sterol':'Sterol','Vitamin':'Other','Vitmain':'Other',
            'Coenzyme':'Other','Other sphingolipid':'SPB','Ceramide':'Cer','Galactosylceramide':'GalCer',
            'Cardiolipin':'CL','Glycerophosphocholine':'PC'}
CLSMAP = {'Fatty acyls':'FA','Terpenoid':'Terpenoid','Sterol':'Sterol','Prenol':'Terpenoid',
          'Sphingolipid':'SPB','Betaine lipid':'DGTS','Ether lipid':'AEG'}
meta = {}
for r in csv.DictReader(open(os.path.join(REPO,'data/metadata/final_lipid_classes.csv'))):
    meta[r['Lipids'].strip()] = (r['Class'].strip(), r['SubClass'].strip())

def pheno_classes(p, layer):
    p = p.strip()
    if layer == 'sumratio':
        if p in ABBR: return {p}
        m = re.match(r'^Sum_(.+?)_over_(.+?)_log10safe$', p)
        if m: return {m.group(1), m.group(2)}
        m = re.match(r'^Sum_(.+?)(_log10safe)?$', p)
        return {m.group(1)} if m else set()
    m = re.match(r'^([A-Za-z0-9]+)\(', p)
    if m and m.group(1) in ABBR: return {m.group(1)}
    cs = meta.get(p)
    if cs:
        cls, sub = cs
        if sub in FALLBACK: return {FALLBACK[sub]}
        if cls in CLSMAP:  return {CLSMAP[cls]}
        return {cls}
    return {'Other'}

rows = list(csv.DictReader(open(os.path.join(REPO,'data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv')), delimiter='\t'))
recs=[]
for x in rows:
    cls=set()
    for p in x['Phenotypes'].split(';'):
        if p.strip(): cls |= pheno_classes(p, x['layer'])
    recs.append(dict(condition=x['condition'], layer=x['layer'], gene=x['GeneID'],
                     name=x['GeneName'], n_phen=int(x['N_Phenotypes']), classes=cls))
df=pd.DataFrame(recs)

# ---------- 1. OVERALL ----------
out=[]
for lab, sel in [('All layers',df),('Individual lipids',df[df.layer=='individual']),
                 ('Class sums / ratios',df[df.layer=='sumratio'])]:
    A=set(sel.loc[sel.condition=='CTL','gene']); B=set(sel.loc[sel.condition=='LIN','gene'])
    s=overlap_stats(A,B,N); s.update(Layer=lab, CTL_only=len(A-B), LIN_only=len(B-A)); out.append(s)
ov=pd.DataFrame(out)[['Layer','n_A','n_B','n_shared','CTL_only','LIN_only','expected','fold','jaccard','p']]
ov.columns=['Layer','n_CTL','n_LIN','n_shared','CTL_only','LIN_only','expected_shared','fold_enrichment','jaccard','p_hypergeom']
ov.to_csv(os.path.join(TAB,'gwas_overlap_overall.csv'), index=False)
print('=== OVERALL (gene universe N = %d) ===' % N)
print(ov.to_string(index=False, float_format=lambda x:f'{x:.4g}'))

# ---------- 2. CLASS-CONCORDANCE OF SHARED GENES ----------
print('\n=== Of shared genes, how many are shared for the SAME lipid class? ===')
conc=[]
for lab, sel in [('All layers',df),('Individual lipids',df[df.layer=='individual'])]:
    cc=sel.groupby(['condition','gene']).classes.apply(lambda s:set().union(*s)).unstack(0)
    sh=cc.dropna()
    same=sum(1 for _,r in sh.iterrows() if r['CTL'] & r['LIN'])
    conc.append(dict(Layer=lab, n_shared=len(sh), same_class=same,
                     pct_same=100*same/len(sh) if len(sh) else np.nan))
conc=pd.DataFrame(conc); print(conc.to_string(index=False, float_format=lambda x:f'{x:.1f}'))
conc.to_csv(os.path.join(TAB,'gwas_overlap_class_concordance.csv'), index=False)

# ---------- 3. PER-CLASS, LAYER-STRATIFIED ----------
def gset(cond, cls, layer):
    s=df[(df.condition==cond)&(df.layer==layer)]
    return set(s.loc[s.classes.map(lambda c: cls in c),'gene'])

per=[]
for layer in ['individual','sumratio']:
    cl=sorted({c for s in df.loc[df.layer==layer,'classes'] for c in s})
    for c in cl:
        A,B=gset('CTL',c,layer), gset('LIN',c,layer)
        if len(A)<3 or len(B)<3: continue
        s=overlap_stats(A,B,N); s.update(Class=c, Layer=layer, CTL_only=len(A-B), LIN_only=len(B-A))
        per.append(s)
per=pd.DataFrame(per)
for L in per.Layer.unique():
    per.loc[per.Layer==L,'q_BH']=bh(per.loc[per.Layer==L,'p'].values)
per=per[['Layer','Class','n_A','n_B','n_shared','CTL_only','LIN_only','expected','fold','jaccard','p','q_BH']]
per.columns=['Layer','Class','n_CTL','n_LIN','n_shared','CTL_only','LIN_only','expected_shared','fold_enrichment','jaccard','p_hypergeom','q_BH']
per=per.sort_values(['Layer','n_CTL'],ascending=[True,False])
per.to_csv(os.path.join(TAB,'gwas_overlap_by_class.csv'), index=False)
print('\n=== PER-CLASS, individual-lipid layer ===')
print(per[per.Layer=='individual'].drop(columns='Layer').to_string(index=False, float_format=lambda x:f'{x:.4g}'))
print('\n=== PER-CLASS, class-sum/ratio layer ===')
print(per[per.Layer=='sumratio'].drop(columns='Layer').to_string(index=False, float_format=lambda x:f'{x:.4g}'))

# ---------- FIGURE ----------
def theme(ax,xlab=None,ylab=None,title=None):
    ax.set_facecolor('white')
    for s in ('top','right'): ax.spines[s].set_visible(False)
    for s in ('left','bottom'): ax.spines[s].set_color('black'); ax.spines[s].set_linewidth(.8)
    ax.grid(False); ax.tick_params(colors='black',labelsize=10,length=3,width=.8)
    if xlab: ax.set_xlabel(xlab,fontsize=12,fontweight='bold')
    if ylab: ax.set_ylabel(ylab,fontsize=12,fontweight='bold')
    if title: ax.set_title(title,fontsize=12.5,fontweight='bold',loc='left')

pi=per[per.Layer=='individual'].set_index('Class')
MAIN=[c for c in ['TG','PC','DG','PE','MGDG','DGDG','SQDG','MG','PG','PA','PS','LPC','LPE','FA','Cer','SM','GalCer','AEG','SPB','Terpenoid','Sterol'] if c in pi.index]
pi=pi.reindex(MAIN)

fig=plt.figure(figsize=(15.5,10.5))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.35],hspace=.30,wspace=.26)

ax=fig.add_subplot(gs[0,0]); ax.set_aspect('equal'); ax.axis('off')
r0=ov.iloc[0]; a,b,s=int(r0.n_CTL),int(r0.n_LIN),int(r0.n_shared)
ax.add_patch(Circle((-.30,0),.40,color='#3B0F70',alpha=.60,ec='black',lw=1.2))
ax.add_patch(Circle(( .32,0),.60,color='#FDE725',alpha=.65,ec='black',lw=1.2))
ax.text(-.62,0,f'{a-s}',ha='center',va='center',fontsize=18,fontweight='bold',color='white')
ax.text( .03,0,f'{s}',  ha='center',va='center',fontsize=18,fontweight='bold')
ax.text( .60,0,f'{b-s}',ha='center',va='center',fontsize=18,fontweight='bold')
ax.text(-.62,.50,'CTL',ha='center',fontsize=14,fontweight='bold',color='#3B0F70')
ax.text( .60,.68,'LIN',ha='center',fontsize=14,fontweight='bold',color='#9A8500')
ax.set_xlim(-1.15,1.15); ax.set_ylim(-.9,.9)
ax.set_title(f'A  Candidate genes, CTL vs LIN (all layers)\n'
             f'{s} shared vs {r0.expected_shared:.0f} expected by chance  '
             f'({r0.fold_enrichment:.1f}$\\times$, $P$ = {r0.p_hypergeom:.0e})\n'
             f'Jaccard = {r0.jaccard:.3f}  —  only {100*r0.jaccard:.1f}% of the union is shared',
             fontsize=11.5,fontweight='bold')

ax=fig.add_subplot(gs[0,1])
y=np.arange(len(ov))[::-1]
ax.barh(y,ov.fold_enrichment,color='#0072B2',edgecolor='black',lw=.6,height=.5)
ax.axvline(1,color='#D55E00',ls='--',lw=1.5)
for yi,(_,r) in zip(y,ov.iterrows()):
    ax.text(r.fold_enrichment+.15,yi,f'{int(r.n_shared)} shared · Jaccard {r.jaccard:.3f} · $P$={r.p_hypergeom:.0e}',
            va='center',fontsize=9.5)
ax.set_yticks(y); ax.set_yticklabels(ov.Layer,fontsize=11)
theme(ax,xlab='Fold enrichment of overlap over chance',
      title='B  Overlap exceeds chance but stays small')
ax.set_xlim(0,ov.fold_enrichment.max()*2.4)
ax.text(1.05,-0.55,'chance',color='#D55E00',fontsize=9,rotation=90,va='bottom')

ax=fig.add_subplot(gs[1,0])
y=np.arange(len(pi))[::-1]
ax.barh(y,pi.CTL_only,color='#3B0F70',edgecolor='black',lw=.4,height=.6,label='CTL only')
ax.barh(y,pi.n_shared,left=pi.CTL_only,color='#1B7837',edgecolor='black',lw=.4,height=.6,label='Shared')
ax.barh(y,pi.LIN_only,left=pi.CTL_only+pi.n_shared,color='#FDE725',edgecolor='black',lw=.4,height=.6,label='LIN only')
tot=(pi.CTL_only+pi.n_shared+pi.LIN_only)
for yi,(cl,r) in zip(y,pi.iterrows()):
    st='***' if r.q_BH<.001 else '**' if r.q_BH<.01 else '*' if r.q_BH<.05 else ''
    ax.text(tot[cl]+22,yi,f'{int(r.n_shared)} shared ({r.fold_enrichment:.0f}$\\times${st})',va='center',fontsize=9)
ax.set_yticks(y); ax.set_yticklabels(pi.index,fontsize=10.5)
theme(ax,xlab='LD-mapped candidate genes (individual-lipid GWAS)',
      title='C  Per lipid class, the two environments barely share loci')
ax.legend(frameon=False,fontsize=10,loc='lower right'); ax.set_xlim(0,tot.max()*1.42)

ax=fig.add_subplot(gs[1,1])
J=np.full((len(MAIN),len(MAIN)),np.nan)
for i,a_ in enumerate(MAIN):
    A=gset('CTL',a_,'individual')
    for j,b_ in enumerate(MAIN):
        B=gset('LIN',b_,'individual')
        J[i,j]=len(A&B)/len(A|B) if (A or B) else np.nan
im=ax.imshow(J,cmap='magma_r',vmin=0,vmax=np.nanmax(J))
ax.set_xticks(range(len(MAIN))); ax.set_xticklabels(MAIN,rotation=90,fontsize=9)
ax.set_yticks(range(len(MAIN))); ax.set_yticklabels(MAIN,fontsize=9)
for i in range(len(MAIN)):
    ax.add_patch(plt.Rectangle((i-.5,i-.5),1,1,fill=False,edgecolor='#00A0FF',lw=1.8))
theme(ax,xlab='LIN candidate set',ylab='CTL candidate set',
      title='D  Cross-class similarity (blue boxes = same class)')
cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.03); cb.set_label('Jaccard index',fontsize=10)
pd.DataFrame(J,index=[f'CTL_{c}' for c in MAIN],columns=[f'LIN_{c}' for c in MAIN]).to_csv(
    os.path.join(TAB,'gwas_overlap_crossclass_jaccard.csv'))

fig.suptitle('CTL and LIN GWAS identify largely non-overlapping candidate genes',
             fontsize=17,fontweight='bold',y=.985)
fig.savefig(os.path.join(FIG,'FigOverlap_CTL_vs_LIN.png'),dpi=300,facecolor='white',bbox_inches='tight')
plt.close(fig)
print('\nWrote figure + 4 tables.')
