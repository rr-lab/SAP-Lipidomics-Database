import pandas as pd, re, os
BASE=os.path.expanduser('~/mnt/SAP-Lipidomics-Database')
CLASSES=['DG','DGDG','LPC','LPE','MG','MGDG','PA','PC','PE','PG','PS','SQDG','TG']
def cls(n):
    m=re.match(r'^([A-Za-z0-9]+)\(', n.strip())
    return m.group(1) if m else None
out=[]
for cond,f in [('CTL','Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv'),
               ('LIN','Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv')]:
    d=pd.read_csv(os.path.join(BASE,'data/raw',f))
    meta=['LineRaw','PlotID','row','col']
    lip=[c for c in d.columns if c not in meta]
    X=d[lip].apply(pd.to_numeric,errors='coerce').fillna(0)
    tic=X.sum(axis=1)
    rec=pd.DataFrame({'Condition':cond,'Taxa':d['LineRaw'],'Class':'TotalLipid','sum_raw':tic,'pct_TIC':100.0})
    out.append(rec)
    for c in CLASSES:
        cols=[l for l in lip if cls(l)==c]
        if not cols: continue
        s=X[cols].sum(axis=1)
        out.append(pd.DataFrame({'Condition':cond,'Taxa':d['LineRaw'],'Class':c,
                                 'sum_raw':s,'pct_TIC':100*s/tic}))
    print(f'# {cond}: n={len(d)} lipid_cols={len(lip)} assigned={sum(1 for l in lip if cls(l) in CLASSES)}',flush=True)
res=pd.concat(out,ignore_index=True)
res.to_csv(os.path.expanduser('~/tmp/soldx/class_sums.csv'),index=False,float_format='%.6g')
print('# rows',len(res))
