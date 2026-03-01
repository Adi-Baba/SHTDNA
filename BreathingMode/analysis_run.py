import numpy as np
import pandas as pd
from pathlib import Path
import json
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import pearsonr
import argparse


def clean_metadata(meta_csv):
    df = pd.read_csv(meta_csv)
    df = df.rename(columns={df.columns[0]: 'PDB'})
    # Clean resolution like [2.3]
    def parse_res(x):
        try:
            if pd.isna(x):
                return np.nan
            s = str(x).strip()
            s = s.strip('[]')
            return float(s)
        except Exception:
            return np.nan
    if 'Resolution' in df.columns:
        df['Resolution_val'] = df['Resolution'].apply(parse_res)
    else:
        df['Resolution_val'] = np.nan
    return df


def bootstrap_struct_corr(df, n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(df)
    stats = []
    x = df['at_content'].values
    y = df['std_beta'].values
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        r, _ = pearsonr(x[idx], y[idx])
        stats.append(r)
    arr = np.array(stats)
    return np.percentile(arr, [2.5, 50, 97.5]), arr


def bootstrap_variance_ratio(site_files, n_boot=1000, seed=0):
    # Preload site arrays into memory (faster than repeated np.load)
    at_arrays = []
    gc_arrays = []
    for p in site_files:
        arr = np.load(p)
        at = arr.get('beta_at', np.array([]))
        gc = arr.get('beta_gc', np.array([]))
        if at.size > 0:
            at_arrays.append(at)
        else:
            at_arrays.append(np.array([]))
        if gc.size > 0:
            gc_arrays.append(gc)
        else:
            gc_arrays.append(np.array([]))

    rng = np.random.default_rng(seed)
    n = len(at_arrays)
    stats = np.empty(n_boot)
    for k in range(n_boot):
        idx = rng.integers(0, n, n)
        # gather concatenated samples
        at_list = [at_arrays[i] for i in idx if at_arrays[i].size > 0]
        gc_list = [gc_arrays[i] for i in idx if gc_arrays[i].size > 0]
        if len(at_list) == 0 or len(gc_list) == 0:
            stats[k] = np.nan
            continue
        at_all = np.concatenate(at_list)
        gc_all = np.concatenate(gc_list)
        if at_all.size == 0 or gc_all.size == 0:
            stats[k] = np.nan
            continue
        stats[k] = np.std(at_all, ddof=1) / np.std(gc_all, ddof=1)
    return np.nanpercentile(stats, [2.5, 50, 97.5]), stats


def compute_effect_sizes(all_at, all_gc):
    # variance ratio
    vr = np.std(all_at, ddof=1) / np.std(all_gc, ddof=1)
    # Cohen's d for means
    n1, n2 = len(all_at), len(all_gc)
    m1, m2 = np.mean(all_at), np.mean(all_gc)
    s1, s2 = np.var(all_at, ddof=1), np.var(all_gc, ddof=1)
    pooled_sd = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    cohen_d = (m1 - m2) / pooled_sd if pooled_sd > 0 else np.nan
    return {'variance_ratio': float(vr), 'cohen_d': float(cohen_d),
            'mean_at': float(m1), 'mean_gc': float(m2),
            'std_at': float(np.std(all_at, ddof=1)), 'std_gc': float(np.std(all_gc, ddof=1))}


def run_mixed_effects(site_files, out_dir):
    # Build site-level dataframe
    rows = []
    for p in site_files:
        pdb = Path(p).stem.replace('_sites', '')
        arr = np.load(p)
        if 'beta_at' in arr:
            for v in arr['beta_at']:
                rows.append({'pdb_id': pdb, 'base': 'AT', 'beta': float(v)})
        if 'beta_gc' in arr:
            for v in arr['beta_gc']:
                rows.append({'pdb_id': pdb, 'base': 'GC', 'beta': float(v)})
    df_site = pd.DataFrame(rows)
    df_site = df_site.dropna()
    df_site['base_type'] = df_site['base'].map({'AT': 1, 'GC': 0})

    # Mixed effects: beta ~ base_type + (1|pdb_id)
    md = smf.mixedlm('beta ~ base_type', df_site, groups=df_site['pdb_id'])
    mdf = md.fit(reml=True)
    out = mdf.summary().as_text()
    (out_dir / 'mixed_effects_summary.txt').write_text(out)
    return mdf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--breathing-csv', default='BreathingMode/breathing_results.csv')
    parser.add_argument('--meta', default='PDB/dataset/hard_dataset/metadata.csv')
    parser.add_argument('--site-dir', default='BreathingMode/site_data')
    parser.add_argument('--out-dir', default='BreathingMode/analysis_results')
    parser.add_argument('--nboot', type=int, default=1000)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.breathing_csv)
    meta = clean_metadata(args.meta)
    # Merge metadata (upper-case PDB codes)
    df = df.merge(meta[['PDB', 'Resolution_val', 'Method']], left_on='pdb_id', right_on='PDB', how='left')

    # Bootstrap Pearson r
    ci_r, r_arr = bootstrap_struct_corr(df, n_boot=args.nboot, seed=1)
    summary = {'pearson_r_ci': ci_r.tolist(), 'pearson_r_median': float(np.median(r_arr))}

    # Gather site files
    site_dir = Path(args.site_dir)
    site_files = sorted([str(p) for p in site_dir.glob('*_sites.npz')])

    # Bootstrap variance ratio
    ci_vr, vr_arr = bootstrap_variance_ratio(site_files, n_boot=args.nboot, seed=2)
    summary['variance_ratio_ci'] = ci_vr.tolist()
    summary['variance_ratio_median'] = float(np.nanmedian(vr_arr))

    # Compute pooled effect sizes
    at_all = []
    gc_all = []
    for p in site_files:
        arr = np.load(p)
        if 'beta_at' in arr and arr['beta_at'].size > 0:
            at_all.append(arr['beta_at'])
        if 'beta_gc' in arr and arr['beta_gc'].size > 0:
            gc_all.append(arr['beta_gc'])
    if len(at_all) == 0 or len(gc_all) == 0:
        print('No site data available for effect sizes')
    else:
        at_all = np.concatenate(at_all)
        gc_all = np.concatenate(gc_all)
        es = compute_effect_sizes(at_all, gc_all)
        summary['effect_sizes'] = es

    # Stratifications
    df['seq_len_bin'] = pd.qcut(df['seq_len'].fillna(0), q=4, duplicates='drop')
    df['at_content_bin'] = pd.qcut(df['at_content'].fillna(0), q=4, duplicates='drop')
    # Resolution category
    def res_cat(x):
        try:
            if np.isnan(x):
                return 'unknown'
            x = float(x)
            if x <= 2.0:
                return 'high'
            elif x <= 3.0:
                return 'mid'
            else:
                return 'low'
        except Exception:
            return 'unknown'
    df['res_cat'] = df['Resolution_val'].apply(res_cat)

    strat_results = {}
    for col in ['seq_len_bin', 'res_cat', 'at_content_bin', 'Method']:
        groups = df.groupby(col)
        gstats = {}
        for name, g in groups:
            if g.shape[0] < 5:
                continue
            r = g['at_content'].corr(g['std_beta'])
            gstats[str(name)] = {'n_structures': int(g.shape[0]), 'pearson_r': float(r)}
        strat_results[col] = gstats
    summary['stratifications'] = strat_results

    # Save summary
    (out_dir / 'analysis_summary.json').write_text(json.dumps(summary, indent=2))

    # Mixed-effects
    try:
        mdf = run_mixed_effects(site_files, out_dir)
        summary['mixed_effects_coef'] = mdf.params.to_dict()
        (out_dir / 'analysis_summary.json').write_text(json.dumps(summary, indent=2))
    except Exception as e:
        (out_dir / 'mixed_effects_error.txt').write_text(str(e))

    print('Analysis complete. Results written to', out_dir)


if __name__ == '__main__':
    main()
