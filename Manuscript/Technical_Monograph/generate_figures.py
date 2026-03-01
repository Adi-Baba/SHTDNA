"""
Publication-quality figure generator for SHT DNA Technical Monograph.
Clean academic style: minimal ink, clear typography, Nature/Science aesthetics.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import os
import warnings
warnings.filterwarnings('ignore')

# ── Output directory ──────────────────────────────────────────────────────────
OUTDIR = "/home/aditya/Documents/SHTDNA/Manuscript/Technical_Monograph/images"
os.makedirs(OUTDIR, exist_ok=True)
os.chdir("/home/aditya/Documents/SHTDNA")

# ── Academic style ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':         'serif',
    'font.serif':          ['DejaVu Serif', 'Times New Roman', 'Georgia'],
    'mathtext.fontset':    'dejavuserif',
    'font.size':           9,
    'axes.titlesize':      10,
    'axes.titleweight':    'bold',
    'axes.labelsize':      9,
    'xtick.labelsize':     8,
    'ytick.labelsize':     8,
    'legend.fontsize':     8,
    'figure.dpi':          300,
    'savefig.dpi':         300,
    'savefig.bbox':        'tight',
    'savefig.pad_inches':  0.06,
    'savefig.facecolor':   'white',
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    'axes.linewidth':      0.7,
    'axes.labelpad':       4,
    'xtick.major.width':   0.7,
    'ytick.major.width':   0.7,
    'xtick.major.size':    3,
    'ytick.major.size':    3,
    'xtick.direction':     'out',
    'ytick.direction':     'out',
    'axes.grid':           True,
    'grid.color':          '#E5E5E5',
    'grid.linewidth':      0.5,
    'lines.linewidth':     1.4,
    'patch.linewidth':     0.6,
    'legend.framealpha':   0.95,
    'legend.edgecolor':    '#CCCCCC',
    'legend.frameon':      True,
    'legend.handlelength': 1.2,
    'figure.facecolor':    'white',
    'axes.facecolor':      '#FAFAFA',
})

# ── Colour palette (desaturated, print-safe) ──────────────────────────────────
C = {
    'nucleosome': '#B2182B',
    'a_tract':    '#2166AC',
    'damaged':    '#D6604D',
    'other':      '#878787',
    'protein':    '#4D9221',
    'highlight':  '#762A83',
    'regression': '#454545',
    'band':       '#BBBBBB',
    'null':       '#9ECAE1',
    'perm_line':  '#B2182B',
    'bar1':       '#2166AC',
    'bar2':       '#4393C3',
    'bar3':       '#92C5DE',
}

MKR = {'Nucleosome': '*', 'A-tract': 's', 'Damaged': 'D', 'Other': 'o'}
SZ  = {'Nucleosome': 110, 'A-tract': 42,  'Damaged': 65,  'Other': 22}
COL = {
    'Nucleosome': C['nucleosome'],
    'A-tract':    C['a_tract'],
    'Damaged':    C['damaged'],
    'Other':      C['other'],
}
NUCLEOSOMES = {'1AOI', '1F66'}
A_TRACTS    = {'10MH','1BF4','1BY4','1DE9','1D02','1DH3','1DNK','1BNK','1BNZ','1CA5'}
DAMAGED     = {'1AHD'}

# ── Load & merge data ─────────────────────────────────────────────────────────
vij = (pd.read_csv("CurvatureMode/curvature_results.csv")[['pdb_id','mean_kappa']]
         .rename(columns={'mean_kappa': 'vij'}))
sha = (pd.read_csv("EccentricityMode/eccentricity_results.csv")[['pdb_id','mean_epsilon']]
         .rename(columns={'mean_epsilon': 'sha'}))
mam = (pd.read_csv("RoughnessMode/roughness_results.csv")[['pdb_id','mean_omega']]
         .rename(columns={'mean_omega': 'mam'}))
df  = vij.merge(sha, on='pdb_id').merge(mam, on='pdb_id').dropna().reset_index(drop=True)

def bio_class(pid):
    if pid in NUCLEOSOMES: return 'Nucleosome'
    if pid in A_TRACTS:    return 'A-tract'
    if pid in DAMAGED:     return 'Damaged'
    return 'Other'

df['cls'] = df['pdb_id'].apply(bio_class)
N = len(df)
print(f"Loaded {N} structures. Generating figures ...\n")


def _legend_h(classes, fmkr=None, fcol=None):
    fmkr = fmkr or MKR
    fcol = fcol or COL
    return [
        Line2D([0],[0], marker=fmkr[c], color='w',
               markerfacecolor=fcol[c], markersize=6.5,
               markeredgewidth=0.4, markeredgecolor='#444', label=c)
        for c in classes
    ]

def _annotate(ax, xs, ys, pids, offset=(5, -13)):
    for i, pid in enumerate(pids):
        if pid in NUCLEOSOMES | DAMAGED:
            col = C['nucleosome'] if pid in NUCLEOSOMES else C['damaged']
            ax.annotate(pid, xy=(xs[i], ys[i]),
                        xytext=offset, textcoords='offset points',
                        fontsize=7, color=col, fontweight='bold',
                        arrowprops=dict(arrowstyle='-', color=col, lw=0.5))


# ─────────────────────────────────────────────────────────────────────────────
# FIG 1  Pairwise scatter matrix
# ─────────────────────────────────────────────────────────────────────────────
def fig_pairwise_scatter():
    pairs = [
        ('vij', 'mam',  r'$\bar{v}_0$  (L$^{-1}$)',    r'$\bar{\Omega}$  (L$^3$)'),
        ('vij', 'sha',  r'$\bar{v}_0$  (L$^{-1}$)',    r'$\bar{\eta}_r$'),
        ('sha', 'mam',  r'$\bar{\eta}_r$',              r'$\bar{\Omega}$  (L$^3$)'),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.6),
                              gridspec_kw={'wspace': 0.40})

    for ax, (xc, yc, xl, yl) in zip(axes, pairs):
        for cls in ['Other', 'A-tract', 'Damaged', 'Nucleosome']:
            sub = df[df['cls'] == cls]
            ax.scatter(sub[xc], sub[yc],
                       c=COL[cls], marker=MKR[cls], s=SZ[cls],
                       alpha=0.85, linewidths=0.3, edgecolors='white', zorder=3)

        r, p   = stats.pearsonr(df[xc], df[yc])
        m, b   = np.polyfit(df[xc], df[yc], 1)
        xl_arr = np.linspace(df[xc].min(), df[xc].max(), 200)
        lc     = C['highlight'] if abs(r) > 0.5 else '#AAAAAA'
        ax.plot(xl_arr, m * xl_arr + b, color=lc, lw=1.5, ls='--', zorder=2)

        pstr = r'$p \approx 0$' if p < 1e-10 else f'$p = {p:.2f}$'
        fc   = '#FFF5F0' if abs(r) > 0.3 else 'white'
        ax.text(0.96, 0.96, f'$r = {r:.3f}$\n{pstr}',
                transform=ax.transAxes, ha='right', va='top', fontsize=7.5,
                bbox=dict(boxstyle='round,pad=0.28', fc=fc, ec='#CCCCCC', lw=0.6))
        ax.set_xlabel(xl); ax.set_ylabel(yl)

    _annotate(axes[0], df['vij'].values, df['mam'].values, df['pdb_id'].values)

    axes[2].legend(handles=_legend_h(['Nucleosome','A-tract','Damaged','Other']),
                   loc='upper right', handletextpad=0.3)
    fig.suptitle(r'Pairwise Relationships Among SHT Geometric Units ($N=65$)',
                 y=1.01, fontsize=10)
    fig.savefig(f'{OUTDIR}/pairwise_scatter.png')
    plt.close()
    print('  ✓ pairwise_scatter.png')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 2  Correlation matrix heatmap
# ─────────────────────────────────────────────────────────────────────────────
def fig_correlation_matrix():
    sub  = df[['vij', 'sha', 'mam']]
    corr = sub.corr()
    pmat = np.ones((3, 3))
    cols = list(corr.columns)
    for i, ci in enumerate(cols):
        for j, cj in enumerate(cols):
            if i != j:
                _, pmat[i, j] = stats.pearsonr(sub[ci], sub[cj])

    tick_labels = [r'Vij ($v_0$)', r'Shambhian ($\eta_r$)', r'Mamton ($\Omega$)']

    fig, ax = plt.subplots(figsize=(3.9, 3.5))
    ax.set_facecolor('white')
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im   = ax.imshow(corr.values, cmap='RdBu_r', norm=norm, aspect='equal')

    ax.set_xticks(range(3)); ax.set_yticks(range(3))
    ax.set_xticklabels(tick_labels, rotation=32, ha='right', fontsize=8)
    ax.set_yticklabels(tick_labels, fontsize=8)

    for i in range(3):
        for j in range(3):
            rv = corr.values[i, j]
            pv = pmat[i, j]
            star = '***' if pv < 0.001 else ('**' if pv < 0.01 else ('*' if pv < 0.05 else 'ns'))
            tc   = 'white' if abs(rv) > 0.55 else '#222222'
            ax.text(j, i, f'{rv:.3f}\n({star})',
                    ha='center', va='center', fontsize=8.5, color=tc, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, fraction=0.044, pad=0.04, shrink=0.88)
    cbar.set_label('Pearson $r$', fontsize=8)
    cbar.ax.tick_params(labelsize=7.5)
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(r'Inter-Unit Correlation Matrix ($N = 65$)', fontsize=10, pad=10)

    plt.tight_layout()
    fig.savefig(f'{OUTDIR}/correlation_matrix.png')
    plt.close()
    print('  ✓ correlation_matrix.png')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 3  PCA manifold + scree
# ─────────────────────────────────────────────────────────────────────────────
def fig_pca_manifold():
    X   = StandardScaler().fit_transform(df[['vij', 'sha', 'mam']].values)
    pca = PCA(n_components=3)
    P   = pca.fit_transform(X)
    ev  = pca.explained_variance_ratio_ * 100

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.1),
                              gridspec_kw={'width_ratios': [1.6, 1], 'wspace': 0.38})

    # Panel A – scatter
    ax = axes[0]
    for cls in ['Other', 'A-tract', 'Damaged', 'Nucleosome']:
        idx = (df['cls'] == cls).values
        ax.scatter(P[idx, 0], P[idx, 1],
                   c=COL[cls], marker=MKR[cls], s=SZ[cls],
                   alpha=0.88, linewidths=0.3, edgecolors='white', zorder=3, label=cls)

    ax.axhline(0, color='#CCCCCC', lw=0.5, ls=':')
    ax.axvline(0, color='#CCCCCC', lw=0.5, ls=':')
    _annotate(ax, P[:, 0], P[:, 1], df['pdb_id'].values, offset=(5, -14))
    ax.set_xlabel(f'PC1 — Integrity Mode ({ev[0]:.1f}%)', fontsize=9)
    ax.set_ylabel(f'PC2 — Shape Mode ({ev[1]:.1f}%)', fontsize=9)
    ax.set_title(r'PCA Geometric Manifold ($N = 65$)', fontsize=10)
    ax.legend(handles=_legend_h(['Nucleosome','A-tract','Damaged','Other']),
              loc='upper left', handletextpad=0.3)

    # Panel B – scree
    ax2 = axes[1]
    labels_scree = ['PC1\nIntegrity', 'PC2\nShape', 'PC3\nResidual']
    bars = ax2.bar(labels_scree, ev,
                   color=[C['bar1'], C['bar2'], C['bar3']],
                   width=0.52, edgecolor='white', linewidth=0.6, alpha=0.92, zorder=3)
    for bar, v in zip(bars, ev):
        ax2.text(bar.get_x() + bar.get_width() / 2, v + 0.5,
                 f'{v:.1f}%', ha='center', va='bottom', fontsize=8.5, fontweight='bold')

    ax2r = ax2.twinx()
    ax2r.plot(labels_scree, np.cumsum(ev), 'o-', color='#333333',
              lw=1.4, ms=5, zorder=5, clip_on=False)
    ax2r.set_ylabel('Cumulative (%)', fontsize=8)
    ax2r.set_ylim(0, 115)
    ax2r.spines['top'].set_visible(False)
    ax2r.tick_params(labelsize=7.5)

    ax2.set_ylabel('Variance Explained (%)', fontsize=9)
    ax2.set_title(f'Scree Plot (PC1+PC2 = {ev[0]+ev[1]:.1f}%)', fontsize=10)
    ax2.set_ylim(0, max(ev) * 1.22)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(f'{OUTDIR}/pca_manifold.png')
    plt.close()
    print('  ✓ pca_manifold.png')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 4  Feature importance
# ─────────────────────────────────────────────────────────────────────────────
def fig_feature_importance():
    vm, vs = df['vij'].mean(), df['vij'].std()
    mm, ms = df['mam'].mean(), df['mam'].std()
    sm, ss = df['sha'].mean(), df['sha'].std()

    def auto_label(row):
        vz = (row['vij'] - vm) / vs
        mz = (row['mam'] - mm) / ms
        sz = (row['sha'] - sm) / ss
        if vz > 3.5 and mz > 3.5:    return 'Nucleosome'
        if sz < -2.0:                  return 'Damaged'
        if vz < -0.5 and mz < -0.3:   return 'A-tract'
        if vz > 0.4 or mz > 0.4:      return 'Protein-bound'
        return 'Free duplex'

    df['fl'] = df.apply(auto_label, axis=1)
    X  = df[['vij', 'sha', 'mam']].values
    y  = df['fl'].values
    rf = RandomForestClassifier(n_estimators=400, max_depth=9,
                                 random_state=42, class_weight='balanced')
    rf.fit(X, y)
    imp = rf.feature_importances_

    names  = [r'Vij ($v_0$)', r'Shambhian ($\eta_r$)', r'Mamton ($\Omega$)']
    order  = np.argsort(imp)          # ascending for horizontal bars
    snames = [names[i] for i in order]
    svals  = imp[order] * 100
    sbars  = [C['bar3'], C['bar2'], C['bar1']]

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7),
                              gridspec_kw={'width_ratios': [1, 0.9], 'wspace': 0.40})

    # Panel A – horizontal bar
    ax = axes[0]
    bars = ax.barh(range(3), svals, color=sbars,
                   height=0.52, edgecolor='white', linewidth=0.6, zorder=3)
    for bar, v in zip(bars, svals):
        ax.text(v + 0.3, bar.get_y() + bar.get_height() / 2,
                f'{v:.1f}%', va='center', fontsize=8.5, fontweight='bold')
    ax.set_yticks(range(3))
    ax.set_yticklabels(snames, fontsize=8.5)
    ax.set_xlabel('Gini Importance (%)', fontsize=9)
    ax.set_title('Random Forest Feature Importance\n(400 trees, $N = 65$)', fontsize=10)
    ax.set_xlim(0, max(svals) * 1.28)
    ax.yaxis.grid(False)

    # Panel B – donut chart
    ax2 = axes[1]
    pie_order  = np.argsort(imp)[::-1]
    wedge_cols = [C['bar1'], C['bar2'], C['bar3']]
    wedges, texts, autotexts = ax2.pie(
        imp[pie_order] * 100,
        labels=[names[i] for i in pie_order],
        colors=[wedge_cols[k] for k in range(3)],
        autopct='%1.1f%%',
        startangle=90,
        pctdistance=0.70,
        wedgeprops=dict(linewidth=0.9, edgecolor='white', width=0.52),
    )
    for at in autotexts: at.set_fontsize(8.5); at.set_fontweight('bold')
    for t in texts:      t.set_fontsize(8.5)
    ax2.set_title('Importance Distribution', fontsize=10)

    plt.tight_layout()
    fig.savefig(f'{OUTDIR}/feature_importance.png')
    plt.close()
    print('  ✓ feature_importance.png')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 5  Functional classification space
# ─────────────────────────────────────────────────────────────────────────────
def fig_functional_classification_space():
    X   = StandardScaler().fit_transform(df[['vij', 'sha', 'mam']].values)
    pca = PCA(n_components=2)
    P   = pca.fit_transform(X)
    ev  = pca.explained_variance_ratio_ * 100

    vm, vs = df['vij'].mean(), df['vij'].std()
    mm, ms = df['mam'].mean(), df['mam'].std()
    sm, ss = df['sha'].mean(), df['sha'].std()

    def func_label(row):
        vz = (row['vij'] - vm) / vs
        mz = (row['mam'] - mm) / ms
        sz = (row['sha'] - sm) / ss
        if vz > 3.5 and mz > 3.5:    return 'Nucleosome'
        if sz < -2.0:                  return 'Damaged'
        if vz < -0.5 and mz < -0.3:   return 'A-tract'
        if vz > 0.4 or mz > 0.4:      return 'Protein-bound'
        return 'Free duplex'

    df['fl'] = df.apply(func_label, axis=1)

    FCOL = {'Nucleosome': C['nucleosome'], 'A-tract': C['a_tract'],
            'Damaged': C['damaged'],       'Protein-bound': C['protein'],
            'Free duplex': '#666666'}
    FMKR = {'Nucleosome':'*','A-tract':'s','Damaged':'D','Protein-bound':'v','Free duplex':'o'}
    FSZ  = {'Nucleosome':120,'A-tract':42,'Damaged':70,'Protein-bound':36,'Free duplex':22}

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2),
                              gridspec_kw={'wspace': 0.38})

    cls_order = ['Free duplex','A-tract','Protein-bound','Damaged','Nucleosome']

    # Panel A – PCA space
    ax = axes[0]
    for cls in cls_order:
        idx = (df['fl'] == cls).values
        ax.scatter(P[idx, 0], P[idx, 1],
                   c=FCOL[cls], marker=FMKR[cls], s=FSZ[cls],
                   alpha=0.88, linewidths=0.3, edgecolors='white', zorder=3, label=cls)
    ax.axhline(0, color='#CCCCCC', lw=0.5, ls=':')
    ax.axvline(0, color='#CCCCCC', lw=0.5, ls=':')
    _annotate(ax, P[:, 0], P[:, 1], df['pdb_id'].values)
    ax.set_xlabel(f'PC1 — Integrity Mode ({ev[0]:.1f}%)', fontsize=9)
    ax.set_ylabel(f'PC2 — Shape Mode ({ev[1]:.1f}%)', fontsize=9)
    ax.set_title('Classification in PCA Space', fontsize=10)
    leg_h = [Line2D([0],[0],marker=FMKR[c],color='w',markerfacecolor=FCOL[c],
                    markersize=6.5,markeredgewidth=0.4,markeredgecolor='#444',label=c)
             for c in cls_order]
    ax.legend(handles=leg_h, loc='upper left', handletextpad=0.3, fontsize=7.5)

    # Panel B – Vij vs log Mamton
    ax2 = axes[1]
    xr = df['vij'].values
    yr = np.log10(df['mam'].values + 1)
    for cls in cls_order:
        idx = (df['fl'] == cls).values
        ax2.scatter(xr[idx], yr[idx],
                    c=FCOL[cls], marker=FMKR[cls], s=FSZ[cls],
                    alpha=0.88, linewidths=0.3, edgecolors='white', zorder=3)
    _annotate(ax2, xr, yr, df['pdb_id'].values)
    ax2.set_xlabel(r'$\bar{v}_0$  (L$^{-1}$)', fontsize=9)
    ax2.set_ylabel(r'$\log_{10}(\bar{\Omega}+1)$', fontsize=9)
    ax2.set_title('Classification: Vij vs. Mamton Plane', fontsize=10)

    plt.tight_layout()
    fig.savefig(f'{OUTDIR}/functional_classification_space.png')
    plt.close()
    print('  ✓ functional_classification_space.png')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 6  Unit validation (stiffness / flatness scaling)
# ─────────────────────────────────────────────────────────────────────────────
def fig_stiffness_flatness_scaling():
    rest = df[~df['pdb_id'].isin(A_TRACTS) & ~df['pdb_id'].isin(NUCLEOSOMES)]
    grps_v   = [df[df['pdb_id'].isin(A_TRACTS)]['vij'].values,
                rest['vij'].values,
                df[df['pdb_id'].isin(NUCLEOSOMES)]['vij'].values]
    grps_m   = [df[df['pdb_id'].isin(A_TRACTS)]['mam'].values,
                rest['mam'].values,
                df[df['pdb_id'].isin(NUCLEOSOMES)]['mam'].values]
    pos   = [1, 2, 3]
    glabs = [f'A-tract\n($n={len(grps_v[0])}$)',
             f'Other\n($n={len(grps_v[1])}$)',
             f'Nucleosome\n($n={len(grps_v[2])}$)']
    gcols = [C['a_tract'], C['other'], C['nucleosome']]
    rng   = np.random.default_rng(42)

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0),
                              gridspec_kw={'wspace': 0.38})

    for (ax, grps, ylabel, title, logscale) in [
        (axes[0], grps_v,
         r'$\bar{v}_0$  (L$^{-1}$)',
         r'Curvature Vij  ($p = 1.79\times10^{-8}$)', False),
        (axes[1], grps_m,
         r'$\bar{\Omega}$  (L$^3$, log scale)',
         r'Roughness Mamton  ($p \approx 0$)', True),
    ]:
        vp = ax.violinplot(grps, positions=pos, widths=0.55,
                           showmedians=True, showextrema=False)
        for pc, col in zip(vp['bodies'], gcols):
            pc.set_facecolor(col); pc.set_alpha(0.25)
            pc.set_edgecolor(col); pc.set_linewidth(0.9)
        vp['cmedians'].set_color('#111111')
        vp['cmedians'].set_linewidth(2.0)

        for p2, grp, col in zip(pos, grps, gcols):
            jit = rng.uniform(-0.11, 0.11, len(grp))
            ax.scatter(p2 + jit, grp, c=col, s=18, alpha=0.78,
                       edgecolors='white', linewidths=0.3, zorder=4)

        if logscale:
            ax.set_yscale('log')
            ax.yaxis.grid(True, which='both')
        ax.set_xticks(pos)
        ax.set_xticklabels(glabs, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=9.5)

        if not logscale:
            ymax = max(np.concatenate(grps))
            yo   = ymax * 1.08
            ax.annotate('', xy=(3, yo), xytext=(1, yo),
                        arrowprops=dict(arrowstyle='-', color='#333', lw=0.9))
            ax.text(2, yo * 1.005, '***', ha='center', va='bottom', fontsize=9)

    fig.suptitle('Unit Validation: Group Separation by SHT Geometric Units',
                 y=1.01, fontsize=10)
    plt.tight_layout()
    fig.savefig(f'{OUTDIR}/stiffness_flatness_scaling.png')
    plt.close()
    print('  ✓ stiffness_flatness_scaling.png')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 7  Geometric state map
# ─────────────────────────────────────────────────────────────────────────────
def fig_geometric_state_map():
    x = df['vij'].values
    y = np.log10(df['mam'].values + 1)
    r, _ = stats.pearsonr(x, y)
    m, b  = np.polyfit(x, y, 1)

    # 95 % confidence band via bootstrap
    rng    = np.random.default_rng(42)
    xl_fit = np.linspace(x.min(), x.max(), 200)
    preds  = np.array([
        np.polyval(np.polyfit(x[idx := rng.integers(0, len(x), len(x))], y[idx], 1), xl_fit)
        for _ in range(2000)
    ])
    lo, hi = np.percentile(preds, [2.5, 97.5], axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3),
                              gridspec_kw={'wspace': 0.40})

    # Panel A – state map with regression
    ax = axes[0]
    ax.fill_between(xl_fit, lo, hi, color=C['band'], alpha=0.40, zorder=1,
                    label='95% CI band')
    ax.plot(xl_fit, m * xl_fit + b, color=C['regression'],
            lw=1.6, ls='--', zorder=2, label=f'$r = {r:.3f}$,  $p \\approx 0$')

    for cls in ['Other', 'A-tract', 'Damaged', 'Nucleosome']:
        idx = (df['cls'] == cls).values
        ax.scatter(x[idx], y[idx],
                   c=COL[cls], marker=MKR[cls], s=SZ[cls],
                   alpha=0.90, linewidths=0.35, edgecolors='white', zorder=3, label=cls)

    _annotate(ax, x, y, df['pdb_id'].values, offset=(5, -13))
    ax.set_xlabel(r'$\bar{v}_0$  (L$^{-1}$)', fontsize=9)
    ax.set_ylabel(r'$\log_{10}(\bar{\Omega}+1)$  (L$^3$)', fontsize=9)
    ax.set_title(f'DNA Geometric State Map ($N = {N}$)', fontsize=10)

    leg_h = _legend_h(['Nucleosome','A-tract','Damaged','Other'])
    leg_h += [
        Line2D([0],[0], color=C['band'],      linewidth=5, alpha=0.5, label='95% CI'),
        Line2D([0],[0], ls='--', color=C['regression'], lw=1.4,
               label=f'$r = {r:.3f}$,  $p \\approx 0$'),
    ]
    ax.legend(handles=leg_h, loc='upper left', handletextpad=0.3, fontsize=7.5)

    # Panel B – Shambhian colour-coded
    ax2 = axes[1]
    sc = ax2.scatter(x, df['sha'].values,
                     c=y, cmap='viridis', s=32,
                     alpha=0.90, linewidths=0.3, edgecolors='white', zorder=3)
    cbar = fig.colorbar(sc, ax=ax2, fraction=0.045, pad=0.03, shrink=0.90)
    cbar.set_label(r'$\log_{10}(\bar{\Omega}+1)$', fontsize=8)
    cbar.ax.tick_params(labelsize=7.5)
    _annotate(ax2, x, df['sha'].values, df['pdb_id'].values)
    ax2.set_xlabel(r'$\bar{v}_0$  (L$^{-1}$)', fontsize=9)
    ax2.set_ylabel(r'$\bar{\eta}_r$  (Shambhian)', fontsize=9)
    ax2.set_title('Mamton Magnitude in Vij–Shambhian Plane', fontsize=10)

    plt.tight_layout()
    fig.savefig(f'{OUTDIR}/geometric_state_map.png')
    plt.close()
    print('  ✓ geometric_state_map.png')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 8  Null hypothesis test
# ─────────────────────────────────────────────────────────────────────────────
def fig_null_test():
    rng   = np.random.default_rng(99)
    xobs  = df['vij'].values
    yobs  = df['mam'].values
    r_obs, _ = stats.pearsonr(xobs, yobs)

    null_r = np.array([
        stats.pearsonr(xobs, rng.permutation(yobs))[0]
        for _ in range(10_000)
    ])
    mu, sd = null_r.mean(), null_r.std()
    Z      = (r_obs - mu) / sd

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1),
                              gridspec_kw={'wspace': 0.40})

    # Panel A – null distribution
    ax = axes[0]
    cnts, _, _ = ax.hist(null_r, bins=55, color=C['null'],
                          edgecolor='white', linewidth=0.4,
                          density=True, zorder=2,
                          label=r'Permuted ($n=10{,}000$)')
    xn = np.linspace(null_r.min(), null_r.max(), 300)
    ax.plot(xn, stats.norm.pdf(xn, mu, sd), color='#2166AC',
            lw=1.6, zorder=3, label=r'$\mathcal{N}(\mu,\sigma)$ fit')
    ax.axvline(r_obs, color=C['perm_line'], lw=2.2, zorder=5,
               label=f'Observed $r = {r_obs:.3f}$')

    ymax = cnts.max() * 1.22
    ax.annotate(
        f'$r_{{\\rm obs}} = {r_obs:.3f}$\n$Z = {Z:.1f}$\n$p \\approx 0$',
        xy=(r_obs, ymax * 0.55),
        xytext=(r_obs - (null_r.max()-null_r.min()) * 0.38, ymax * 0.70),
        fontsize=8, color=C['perm_line'], fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=C['perm_line'], lw=1.2))
    ax.set_ylim(0, ymax)
    ax.set_xlabel('Permuted Pearson $r$', fontsize=9)
    ax.set_ylabel('Density', fontsize=9)
    ax.set_title('Null Distribution for Vij–Mamton\nCorrelation', fontsize=10)
    ax.legend(loc='upper left', handletextpad=0.3, fontsize=7.5)

    # Panel B – Z-score on standard normal
    ax2 = axes[1]
    zx  = np.linspace(-5, 6.5, 500)
    ax2.fill_between(zx, stats.norm.pdf(zx), where=(zx < -1.96),
                     color=C['null'], alpha=0.40)
    ax2.fill_between(zx, stats.norm.pdf(zx), where=(zx > 1.96),
                     color=C['null'], alpha=0.40, label='$p < 0.05$ tails')
    ax2.plot(zx, stats.norm.pdf(zx), color='#2166AC', lw=1.8,
             label=r'$\mathcal{N}(0,1)$')
    zplot = min(abs(Z), 5.9) * np.sign(Z)
    ax2.axvline(zplot, color=C['perm_line'], lw=2.2,
                label=f'$Z = {Z:.1f}$ (off-scale)')
    ax2.text(4.1, 0.19, f'$Z = {Z:.1f}\\sigma$', fontsize=9.5,
             color=C['perm_line'], fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.28', fc='#FFF5F0',
                       ec=C['perm_line'], lw=0.7))
    ax2.set_xlabel('$Z$-score', fontsize=9)
    ax2.set_ylabel('Density', fontsize=9)
    ax2.set_title(f'Observed $Z = {Z:.1f}$ on $\\mathcal{{N}}(0,1)$', fontsize=10)
    ax2.legend(loc='upper left', handletextpad=0.3, fontsize=7.5)
    ax2.set_xlim(-5, 7)

    fig.suptitle('Permutation Null Hypothesis Test: Vij–Mamton Coupling',
                 y=1.01, fontsize=10)
    plt.tight_layout()
    fig.savefig(f'{OUTDIR}/null_test_proof.png')
    plt.close()
    print('  ✓ null_test_proof.png')


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    fig_pairwise_scatter()
    fig_correlation_matrix()
    fig_pca_manifold()
    fig_feature_importance()
    fig_functional_classification_space()
    fig_stiffness_flatness_scaling()
    fig_geometric_state_map()
    fig_null_test()
    print(f'\nAll 8 figures saved to {OUTDIR}')
