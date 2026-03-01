"""
arc_length_curvature_comparison.py  (v2)

Uses EXACTLY the same centroid trajectory C(h) as the official Vij code,
then computes BOTH parametrizations from the same curve:

  Vij (h-parametrized) :  kappa_h = ||d2C/dh2||           (existing SHT definition)
  Arc-length curvature :  kappa_s = ||C'' x C'|| / ||C'||^3  (intrinsic, coord-free)

All centroid curves come from compute_sigma_sht_from_pdb() - same code path
as the production curvature_vij_results.csv.
"""

import sys, os, warnings, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr, spearmanr
from pathlib import Path

warnings.filterwarnings('ignore')

BASE       = '/home/aditya/Documents/SHTDNA'
STRUCT_DIR = f'{BASE}/PDB/dataset/hard_dataset/structures'
CODE_DIR   = f'{BASE}/PDB/code'
sys.path.insert(0, CODE_DIR)
import sht_dna_analysis


def kappa_from_C(C, h):
    """
    Compute both h-parametrized and arc-length curvature from centroid curve C(h).
    Returns (mean_kappa_h, mean_kappa_s, mean_dsdh)
    """
    if len(h) < 4:
        return None, None, None

    # Height-parametrized curvature (official Vij)
    dC_dh   = np.gradient(C, h, axis=0)
    d2C_dh2 = np.gradient(dC_dh, h, axis=0)
    kappa_h = np.linalg.norm(d2C_dh2, axis=1)

    # Arc-length parameter
    diffs     = np.diff(C, axis=0)
    seg_lens  = np.linalg.norm(diffs, axis=1)
    s         = np.concatenate([[0.], np.cumsum(seg_lens)])
    if s[-1] < 1e-6:
        return float(np.mean(kappa_h)), None, None

    ds_dh  = np.abs(np.gradient(s, h))

    dC_ds   = np.gradient(C, s, axis=0)
    d2C_ds2 = np.gradient(dC_ds, s, axis=0)
    T  = dC_ds
    T2 = d2C_ds2

    if C.shape[1] >= 3:
        cross   = np.cross(T2[:,:3], T[:,:3])
        kappa_s = np.linalg.norm(cross, axis=1) / (np.linalg.norm(T[:,:3], axis=1)**3 + 1e-12)
    else:
        numer   = np.abs(T[:,0]*T2[:,1] - T[:,1]*T2[:,0])
        denom   = (T[:,0]**2 + T[:,1]**2)**1.5 + 1e-12
        kappa_s = numer / denom

    return float(np.mean(kappa_h)), float(np.nanmean(kappa_s)), float(np.mean(ds_dh))


# Load structure list
vij_df = pd.read_csv(f'{BASE}/CurvatureMode/curvature_vij_results.csv')
pids   = vij_df['pdb_id'].str.replace('.cif','',regex=False).str.replace('.pdb','',regex=False).str.upper().str.strip()

print(f"Processing {len(pids)} structures using SHT centroid trajectories...")
results = []
failed  = 0

for i, pid in enumerate(pids):
    cif_path = f'{STRUCT_DIR}/{pid}.cif'
    if not os.path.exists(cif_path):
        cif_path = f'{STRUCT_DIR}/{pid.lower()}.cif'
    if not os.path.exists(cif_path):
        failed += 1; continue

    try:
        res    = sht_dna_analysis.compute_sigma_sht_from_pdb(str(cif_path), num_slices=100)
        h      = np.asarray(res['h'])
        C_axis = np.asarray(res['C_axis'])
        if len(h) < 4:
            failed += 1; continue
    except Exception:
        failed += 1; continue

    kh, ks, dsdh = kappa_from_C(C_axis, h)
    mean_vij_official = float(np.mean(np.asarray(res['kappa_h']))) if len(res['kappa_h']) else np.nan

    results.append({
        'pdb_id':        pid,
        'kappa_h':       kh,
        'kappa_s':       ks,
        'dsdh':          dsdh,
        'vij_official':  mean_vij_official,
    })

    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{len(pids)} done, {failed} failed")

print(f"Finished. {len(results)} OK, {failed} failed.")

df = pd.DataFrame(results).dropna(subset=['kappa_h', 'kappa_s'])
df.to_csv(f'{BASE}/DiscoveryMode/arc_length_comparison.csv', index=False)
print(f"Saved {len(df)} rows to arc_length_comparison.csv")

# Sanity check
r_sanity, _ = pearsonr(df['kappa_h'], df['vij_official'])
print(f"\nSanity check: r(recomputed kappa_h, official mean_vij) = {r_sanity:.4f}")

r_main,  p_main  = pearsonr(df['kappa_h'], df['kappa_s'])
rho_main, p_rho  = spearmanr(df['kappa_h'], df['kappa_s'])
print(f"\n=== Arc-Length vs Height-Parametrized Curvature ===")
print(f"N = {len(df)}")
print(f"Pearson  r(kappa_h, kappa_s)  = {r_main:.4f}, p = {p_main:.2e}")
print(f"Spearman r(kappa_h, kappa_s)  = {rho_main:.4f}, p = {p_rho:.2e}")
print(f"Mean ds/dh                    = {df['dsdh'].mean():.4f} +/- {df['dsdh'].std():.4f}")

# Predictive power vs Mamton
mamton = pd.read_csv(f'{BASE}/RoughnessMode/roughness_mamton_results.csv')
mamton['pid'] = mamton['pdb_id'].str.replace('.cif','',regex=False).str.upper()
df['pid']     = df['pdb_id'].str.upper()
m = df.merge(mamton[['pid','mean_mamton']], on='pid').dropna()
print(f"\n=== Predictive Power vs Mamton roughness (N = {len(m)}) ===")
r_ks_mamton, p_ks_mamton = pearsonr(m['kappa_s'], m['mean_mamton'])
r_kh_mamton, p_kh_mamton = pearsonr(m['kappa_h'], m['mean_mamton'])
r_vo_mamton, p_vo_mamton = pearsonr(m['vij_official'], m['mean_mamton'])
print(f"r(kappa_s arc-len, Mamton)  = {r_ks_mamton:.4f}, p = {p_ks_mamton:.2e}")
print(f"r(kappa_h h-param, Mamton)  = {r_kh_mamton:.4f}, p = {p_kh_mamton:.2e}")
print(f"r(Vij official,    Mamton)  = {r_vo_mamton:.4f}, p = {p_vo_mamton:.2e}")

# Figure
os.makedirs(f'{BASE}/DiscoveryMode/robustness_figures', exist_ok=True)
fig = plt.figure(figsize=(16, 5))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38)

ax1 = fig.add_subplot(gs[0])
ax1.scatter(df['kappa_h'], df['kappa_s'], alpha=0.6, s=40, color='#1976D2', edgecolors='white', lw=0.3)
xr = np.linspace(df['kappa_h'].min(), df['kappa_h'].max(), 200)
m_, b_ = np.polyfit(df['kappa_h'], df['kappa_s'], 1)
ax1.plot(xr, m_*xr+b_, 'r-', lw=2)
ax1.set_xlabel('kappa_h — Vij (height-param)', fontsize=11)
ax1.set_ylabel('kappa_s — Arc-length curvature', fontsize=11)
ax1.set_title(f'Two Curvature Parametrizations\nr = {r_main:.3f}, N = {len(df)}', fontsize=11)
ax1.spines[['top','right']].set_visible(False)

ax2 = fig.add_subplot(gs[1])
ax2.hist(df['dsdh'], bins=30, color='#43A047', edgecolor='white', alpha=0.85)
ax2.axvline(1.0, color='red', ls='--', lw=2, label='ds/dh=1 (equivalent)')
ax2.axvline(df['dsdh'].mean(), color='navy', ls='--', lw=1.5, label=f"Mean={df['dsdh'].mean():.3f}")
ax2.set_xlabel('ds/dh ratio', fontsize=11)
ax2.set_ylabel('Count', fontsize=11)
ax2.set_title(f"Parametrization Conversion Factor\nMean={df['dsdh'].mean():.3f}+/-{df['dsdh'].std():.3f}", fontsize=11)
ax2.legend(fontsize=9)
ax2.spines[['top','right']].set_visible(False)

ax3 = fig.add_subplot(gs[2])
bars = ax3.bar(['kappa_s\n(arc-len)', 'Vij (h-param)'],
               [abs(r_ks_mamton), abs(r_vo_mamton)],
               color=['#EF5350','#1976D2'], width=0.45, edgecolor='white', lw=1.5)
ax3.set_ylabel('|r| with Mamton roughness', fontsize=11)
ax3.set_title('Predictive Power\nvs Structural Roughness (Mamton)', fontsize=11)
ax3.set_ylim(0, 1.0)
for bar, v in zip(bars, [abs(r_ks_mamton), abs(r_vo_mamton)]):
    ax3.text(bar.get_x()+bar.get_width()/2, v+0.02, f'{v:.3f}',
             ha='center', fontsize=12, fontweight='bold')
ax3.spines[['top','right']].set_visible(False)

fig.suptitle('Arc-Length vs Height-Parametrized Curvature: Compatibility Analysis',
             fontsize=12, fontweight='bold')
fig.savefig(f'{BASE}/DiscoveryMode/robustness_figures/arc_length_comparison.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nFigure saved.")

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"r(kappa_h, kappa_s) = {r_main:.4f} — independence of the two metrics")
print(f"r(Vij, Mamton)      = {r_vo_mamton:.4f} — Vij's coupling strength")
print(f"r(kappa_s, Mamton)  = {r_ks_mamton:.4f} — arc-length's coupling strength")
print(f"Vij signal advantage= {abs(r_vo_mamton)-abs(r_ks_mamton):.4f} over arc-length")
