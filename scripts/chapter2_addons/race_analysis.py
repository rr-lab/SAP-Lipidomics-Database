import sys, os, math
sys.path.insert(0, os.path.expanduser('~/tmp/soldx'))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from stats_lite import kruskal, bh

REPO = os.path.expanduser('~/mnt/SAP-Lipidomics-Database')
FIG  = os.path.join(REPO, 'fig', 'race'); os.makedirs(FIG, exist_ok=True)
TAB  = os.path.join(REPO, 'table', 'race'); os.makedirs(TAB, exist_ok=True)

CLASSES = ['TotalLipid','MGDG','PC','DG','DGDG','TG','PE','SQDG','MG','PG','PA','LPC','PS','LPE']
PURE = ['Bicolor','Caudatum','Durra','Guinea','Kafir']
RACE_ORDER = PURE + ['Mixed']
RACE_COL = {'Bicolor':'#E69F00','Caudatum':'#0072B2','Durra':'#009E73',
            'Guinea':'#CC79A7','Kafir':'#D55E00','Mixed':'#999999'}
KCOL = {1:'#E69F00',2:'#56B4E9',3:'#009E73',4:'#F0E442',5:'#0072B2',6:'#CC79A7'}

# ---------- plot_theme (ggplot theme_minimal(base_size=14) equivalent) ----------
def apply_theme(ax, xlab=None, ylab=None, title=None):
    ax.set_facecolor('white')
    for s in ('top','right'): ax.spines[s].set_visible(False)
    for s in ('left','bottom'):
        ax.spines[s].set_visible(True); ax.spines[s].set_color('black'); ax.spines[s].set_linewidth(0.8)
    ax.grid(False)
    ax.tick_params(colors='black', labelsize=10, length=3, width=0.8)
    if xlab is not None: ax.set_xlabel(xlab, fontsize=12, fontweight='bold')
    if ylab is not None: ax.set_ylabel(ylab, fontsize=12, fontweight='bold')
    if title is not None: ax.set_title(title, fontsize=12, fontweight='bold')

def race_group(x):
    if pd.isna(x) or str(x) == 'NA': return np.nan
    x = str(x)
    if 'verticilliflorum' in x.lower(): return np.nan
    return x if x in PURE else 'Mixed'

d = pd.read_csv(os.path.expanduser('~/tmp/soldx/joined.csv'))
d['RaceGroup'] = d['Original_Race'].map(race_group)
d['TotalLipid'] = pd.to_numeric(d['TotalLipid'], errors='coerce')
d['TotalLipid_log10'] = np.log10(d['TotalLipid'])

def run_tests(df, groupcol, order):
    rows = []
    for cond in ['CTL','LIN']:
        sub = df[(df.Condition == cond) & df[groupcol].notna()]
        for tr in CLASSES:
            col = 'TotalLipid_log10' if tr == 'TotalLipid' else tr
            grps, labs = [], []
            for g in order:
                v = sub.loc[sub[groupcol] == g, col].dropna().values
                if len(v) >= 5: grps.append(v); labs.append(g)
            if len(grps) < 3: continue
            H, p, k, eps2 = kruskal(grps)
            meds = {g: float(np.median(v)) for g, v in zip(labs, grps)}
            hi = max(meds, key=meds.get); lo = min(meds, key=meds.get)
            rows.append(dict(Condition=cond, Grouping=groupcol, Trait=tr, n=int(sum(len(g) for g in grps)),
                             k_groups=k, H=H, p=p, epsilon2=eps2, highest=hi, lowest=lo,
                             **{f'median_{g}': meds.get(g, np.nan) for g in order}))
    r = pd.DataFrame(rows)
    for cond in r.Condition.unique():
        m = r.Condition == cond
        r.loc[m, 'q_BH'] = bh(r.loc[m, 'p'].values)
    return r

res_race = run_tests(d, 'RaceGroup', RACE_ORDER)
res_k    = run_tests(d, 'K.Cluster', [1,2,3,4,5,6])
res = pd.concat([res_race, res_k], ignore_index=True)
res.to_csv(os.path.join(TAB, 'race_structure_lipid_tests.csv'), index=False)

def stars(q):
    return '***' if q < .001 else '**' if q < .01 else '*' if q < .05 else 'ns'

# ================= FIGURE: boxplot grid per condition =================
def boxgrid(cond, outfile):
    sub = d[(d.Condition == cond) & d.RaceGroup.notna()]
    rr = res_race[res_race.Condition == cond].set_index('Trait')
    fig, axes = plt.subplots(4, 4, figsize=(15.5, 13.5))
    axes = axes.ravel()
    for i, tr in enumerate(CLASSES):
        ax = axes[i]
        col = 'TotalLipid_log10' if tr == 'TotalLipid' else tr
        data, labs, cols = [], [], []
        for g in RACE_ORDER:
            v = sub.loc[sub.RaceGroup == g, col].dropna().values
            if len(v) >= 5: data.append(v); labs.append(g); cols.append(RACE_COL[g])
        bp = ax.boxplot(data, patch_artist=True, widths=.62, showfliers=False,
                        medianprops=dict(color='black', lw=1.4),
                        whiskerprops=dict(color='black', lw=.8),
                        capprops=dict(color='black', lw=.8),
                        boxprops=dict(lw=.8, edgecolor='black'))
        for patch, c in zip(bp['boxes'], cols): patch.set_facecolor(c); patch.set_alpha(.75)
        rng = np.random.default_rng(42)
        for j, v in enumerate(data):
            ax.scatter(j + 1 + rng.uniform(-.16, .16, len(v)), v, s=5, color='black', alpha=.35, zorder=3, linewidths=0)
        ylab = 'log10 total lipid signal' if tr == 'TotalLipid' else '% of TIC'
        apply_theme(ax, ylab=ylab)
        ax.set_xticks(range(1, len(labs) + 1))
        ax.set_xticklabels([l[:4] for l in labs], fontsize=9)
        if tr in rr.index:
            q = rr.loc[tr, 'q_BH']; e = rr.loc[tr, 'epsilon2']
            ax.set_title(f"{tr}   $\\varepsilon^2$={e:.3f}  {stars(q)}", fontsize=12, fontweight='bold')
        else:
            ax.set_title(tr, fontsize=12, fontweight='bold')
    for j in range(len(CLASSES), len(axes)): axes[j].axis('off')
    handles = [Line2D([0],[0], marker='s', color='none', markerfacecolor=RACE_COL[g],
                      markeredgecolor='black', markersize=11, label=g) for g in RACE_ORDER]
    fig.legend(handles=handles, loc='lower right', bbox_to_anchor=(.97,.06), ncol=2,
               frameon=False, fontsize=12)
    fig.suptitle(f'Lipid class composition by botanical race — {"Control (CTL)" if cond=="CTL" else "Low input (LIN)"}',
                 fontsize=16, fontweight='bold', y=.995)
    fig.tight_layout(rect=[0,0,1,.985])
    fig.savefig(outfile, dpi=300, facecolor='white'); plt.close(fig)

boxgrid('CTL', os.path.join(FIG, 'FigRace_A_CTL_boxplots.png'))
boxgrid('LIN', os.path.join(FIG, 'FigRace_B_LIN_boxplots.png'))

# ================= FIGURE: variance explained =================
fig, axes = plt.subplots(1, 2, figsize=(12.5, 7), sharey=True)
for ax, (gc, tab, ttl) in zip(axes, [('RaceGroup', res_race, 'Botanical race (6 groups)'),
                                     ('K.Cluster', res_k, 'Genetic cluster (K = 6)')]):
    piv = tab.pivot_table(index='Trait', columns='Condition', values='epsilon2')
    qiv = tab.pivot_table(index='Trait', columns='Condition', values='q_BH')
    piv = piv.reindex([t for t in CLASSES if t in piv.index])
    y = np.arange(len(piv))[::-1]
    ax.barh(y + .19, piv['CTL'], height=.36, color='#3B0F70', edgecolor='black', lw=.5, label='CTL')
    ax.barh(y - .19, piv['LIN'], height=.36, color='#FDE725', edgecolor='black', lw=.5, label='LIN')
    for yi, tr in zip(y, piv.index):
        for off, c in [(.19,'CTL'), (-.19,'LIN')]:
            v = piv.loc[tr, c]; q = qiv.loc[tr, c]
            if np.isfinite(v) and q < .05:
                ax.text(v + .004, yi + off, stars(q), va='center', fontsize=9)
    ax.set_yticks(y); ax.set_yticklabels(piv.index, fontsize=11)
    apply_theme(ax, xlab='Variance explained  ($\\varepsilon^2$)', title=ttl)
    ax.set_xlim(0, max(.02, np.nanmax(piv.values) * 1.32))
    nsig = int((qiv < .05).sum().sum())
    ax.text(.98, .02, ('no trait significant after BH correction' if nsig == 0
                       else f'{nsig} trait(s) significant after BH correction'),
            transform=ax.transAxes, ha='right', va='bottom', fontsize=10.5, style='italic')
axes[0].legend(frameon=False, fontsize=12, loc='center right')
axes[1].annotate('PS is the single exception:\n11% of variance, $q$ = 4e-7',
                 xy=(.110, 1.81), xytext=(.055, 4.6), fontsize=10.5,
                 arrowprops=dict(arrowstyle='->', lw=1.1, color='black'))
fig.suptitle('Botanical race and genetic structure explain almost none of the lipid-class variation',
             fontsize=15, fontweight='bold')
fig.tight_layout(rect=[0,0,1,.95])
fig.savefig(os.path.join(FIG, 'FigRace_C_variance_explained.png'), dpi=300, facecolor='white')
plt.close(fig)

print(res[['Condition','Grouping','Trait','n','H','p','q_BH','epsilon2','highest','lowest']]
      .to_string(index=False, float_format=lambda x: f'{x:.4g}'))
