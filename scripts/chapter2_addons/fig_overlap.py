import sys,os,csv
sys.path.insert(0,os.path.expanduser('~/tmp/soldx')); csv.field_size_limit(10**8)
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
REPO=os.path.expanduser('~/mnt/SAP-Lipidomics-Database')
FIG=os.path.join(REPO,'fig','overlap'); TAB=os.path.join(REPO,'table','overlap')

ovg=pd.read_csv(os.path.join(TAB,'gwas_overlap_overall.csv'))
ovl=pd.read_csv(os.path.join(TAB,'gwas_overlap_locus_level.csv'))
per=pd.read_csv(os.path.join(TAB,'gwas_overlap_by_class.csv'))
inf=pd.read_csv(os.path.join(TAB,'gwas_gene_to_locus_inflation.csv'))

def theme(ax,xlab=None,ylab=None,title=None):
    ax.set_facecolor('white')
    for s in ('top','right'): ax.spines[s].set_visible(False)
    for s in ('left','bottom'): ax.spines[s].set_color('black'); ax.spines[s].set_linewidth(.8)
    ax.grid(False); ax.tick_params(colors='black',labelsize=10,length=3,width=.8)
    if xlab: ax.set_xlabel(xlab,fontsize=12,fontweight='bold')
    if ylab: ax.set_ylabel(ylab,fontsize=12,fontweight='bold')
    if title: ax.set_title(title,fontsize=12.5,fontweight='bold',loc='left')

CT,LN,SH='#3B0F70','#FDE725','#1B7837'
fig=plt.figure(figsize=(16,11.5))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.3],hspace=.33,wspace=.25)

# ---- A: two Venns ----
axA=fig.add_subplot(gs[0,0]); axA.axis('off'); axA.set_xlim(0,2); axA.set_ylim(0,1.15); axA.set_aspect('equal')
def venn(ax,cx,a,b,s,ttl,sub):
    ax.add_patch(Circle((cx-.145,.55),.20,color=CT,alpha=.60,ec='black',lw=1.1))
    ax.add_patch(Circle((cx+.165,.55),.30,color=LN,alpha=.65,ec='black',lw=1.1))
    ax.text(cx-.28,.55,f'{a-s}',ha='center',va='center',fontsize=13,fontweight='bold',color='white')
    ax.text(cx+.005,.55,f'{s}', ha='center',va='center',fontsize=13,fontweight='bold')
    ax.text(cx+.30,.55,f'{b-s}',ha='center',va='center',fontsize=13,fontweight='bold')
    ax.text(cx-.28,.81,'CTL',ha='center',fontsize=11,fontweight='bold',color=CT)
    ax.text(cx+.30,.90,'LIN',ha='center',fontsize=11,fontweight='bold',color='#9A8500')
    ax.text(cx,1.03,ttl,ha='center',fontsize=12,fontweight='bold')
    ax.text(cx,.20,sub,ha='center',va='top',fontsize=9.8,linespacing=1.5)
rg=ovg[ovg.Layer=='All layers'].iloc[0]
rl=ovl[(ovl.Layer=='All layers')&(ovl.window_kb==250)].iloc[0]
venn(axA,.52,int(rg.n_CTL),int(rg.n_LIN),int(rg.n_shared),'Gene level',
     f'{rg.fold_enrichment:.2f}$\\times$ chance\n$P$ = {rg.p_hypergeom:.0e}')
venn(axA,1.50,int(rl.n_CTL),int(rl.n_LIN),int(rl.n_shared),'Locus level (250 kb)',
     f'{rl.fold_enrichment:.2f}$\\times$ chance\n$P$ = {rl.p_hypergeom:.2f}')
axA.set_title('A  Gene-level overlap is an LD artefact',fontsize=13,fontweight='bold',loc='left')

# ---- B: fold enrichment across resolutions ----
axB=fig.add_subplot(gs[0,1])
lays=['Individual lipids','Class sums / ratios']
res=[('Gene',None)]+[('%d kb'%w,w) for w in [100,250,500]]
x=np.arange(len(res)); w=.36
for k,(lay,col) in enumerate(zip(lays,['#0072B2','#D55E00'])):
    vals,ps=[],[]
    for lab,win in res:
        if win is None:
            r=ovg[ovg.Layer==lay].iloc[0]
        else:
            r=ovl[(ovl.Layer==lay)&(ovl.window_kb==win)].iloc[0]
        vals.append(r.fold_enrichment); ps.append(r.p_hypergeom)
    bars=axB.bar(x+(k-.5)*w,vals,width=w,color=col,edgecolor='black',lw=.6,label=lay)
    for xi,v,p in zip(x+(k-.5)*w,vals,ps):
        st='***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'n.s.'
        axB.text(xi,v+.08,st,ha='center',fontsize=9.5,fontweight='bold')
axB.axhline(1,color='black',ls='--',lw=1.3)
axB.set_xticks(x); axB.set_xticklabels([r[0] for r in res],fontsize=11)
theme(axB,xlab='Resolution at which overlap is counted',ylab='Fold enrichment over chance',
      title='B  Collapsing LD blocks removes the signal')
axB.legend(frameon=False,fontsize=10.5); axB.set_ylim(0,max(3.2,ovg.fold_enrichment.max()*1.35))
axB.text(3.42,1.06,'chance',fontsize=9,color='black')

# ---- C: per-class (individual layer) ----
axC=fig.add_subplot(gs[1,0])
pi=per[per.Layer=='individual'].set_index('Class')
order=[c for c in ['TG','PC','FA','Terpenoid','DG','SQDG','AEG'] if c in pi.index]
pi=pi.reindex(order)
y=np.arange(len(pi))[::-1]
axC.barh(y,pi.CTL_only,color=CT,edgecolor='black',lw=.4,height=.6,label='CTL only')
axC.barh(y,pi.n_shared,left=pi.CTL_only,color=SH,edgecolor='black',lw=.4,height=.6,label='Shared')
axC.barh(y,pi.LIN_only,left=pi.CTL_only+pi.n_shared,color=LN,edgecolor='black',lw=.4,height=.6,label='LIN only')
tot=pi.CTL_only+pi.n_shared+pi.LIN_only
for yi,(cl,r) in zip(y,pi.iterrows()):
    st='***' if r.q_BH<.001 else '**' if r.q_BH<.01 else '*' if r.q_BH<.05 else 'n.s.'
    axC.text(tot[cl]+45,yi,f'{int(r.n_shared)} shared · {r.fold_enrichment:.1f}$\\times$ {st}',va='center',fontsize=9.5)
axC.set_yticks(y); axC.set_yticklabels(pi.index,fontsize=11)
theme(axC,xlab='LD-mapped candidate genes (individual-lipid GWAS)',
      title='C  Per lipid class: SQDG, DG and AEG share zero genes')
axC.legend(frameon=False,fontsize=10,loc='lower right'); axC.set_xlim(0,tot.max()*1.45)

# ---- D: gene -> locus inflation ----
axD=fig.add_subplot(gs[1,1])
sub=inf[inf.Layer!='all layers'].copy()
sub['lab']=sub.Condition+'\n'+sub.Layer.str.replace('class sums/ratios','sums/ratios')
x=np.arange(len(sub)); w=.38
axD.bar(x-w/2,sub.n_genes,width=w,color='#999999',edgecolor='black',lw=.6,label='Candidate genes')
axD.bar(x+w/2,sub.n_loci,width=w,color='#009E73',edgecolor='black',lw=.6,label='Independent loci (250 kb)')
for xi,r in zip(x,sub.itertuples()):
    axD.text(xi-w/2,r.n_genes*1.10,f'{r.n_genes}',ha='center',fontsize=10,fontweight='bold')
    axD.text(xi+w/2,r.n_loci*1.10,f'{r.n_loci}',ha='center',fontsize=10,fontweight='bold')
    axD.text(xi,sub.n_genes.max()*3.0,f'{r.genes_per_locus:.1f}\ngenes/locus',ha='center',
             fontsize=9.5,style='italic',linespacing=1.3)
axD.set_xticks(x); axD.set_xticklabels(sub.lab,fontsize=10)
theme(axD,ylab='Count',title='D  Candidate-gene counts are inflated by LD')
axD.set_yscale('log'); axD.set_ylim(10,sub.n_genes.max()*6)
axD.legend(frameon=False,fontsize=10.5,loc='upper left')

fig.suptitle('CTL and LIN GWAS share no more loci than expected by chance',
             fontsize=18,fontweight='bold',y=.985)
fig.savefig(os.path.join(FIG,'FigOverlap_CTL_vs_LIN.png'),dpi=300,facecolor='white',bbox_inches='tight')
print('saved', os.path.join(FIG,'FigOverlap_CTL_vs_LIN.png'))
