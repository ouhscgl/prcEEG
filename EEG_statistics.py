#!/usr/bin/env python3
"""
EEG Statistical Analysis & Visualization
=========================================
Combined statistics and plotting for EEG analysis outputs.

Reads:
    - bandpower_<roi>.csv
    - spectral_slope_<roi>.csv  
    - ersp_<roi>_summary.csv

Outputs:
    - stats_bandpower.csv, stats_ratios.csv, stats_slope.csv, stats_ersp.csv
    - Publication-ready figures (PNG + SVG)

Usage:
    python EEG_statistics.py <data_dir> <output_dir> [roi]
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats as sp_stats
from itertools import combinations
import warnings
import sys

# Optional: statsmodels for mixed effects
try:
    import statsmodels.formula.api as smf
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("Warning: statsmodels not found. Using ANOVA instead of LME.")

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

GROUP_ORDER = ['healthy_controls', 'normal_performers', 'low_performers']
GROUP_LABELS = {
    'healthy_controls': 'Healthy Controls',
    'normal_performers': 'Normal Performers', 
    'low_performers': 'Low Performers',
}

GROUP_COLORS = {
    'healthy_controls': '#77DD77',
    'normal_performers': '#E2B272',
    'low_performers': '#C08080',
}

TASK_ORDER = ['eyes_open', 'nback_0a', 'nback_1a', 'nback_0b', 'nback_2a']
TASK_LABELS = {
    'eyes_open': 'Eyes Open',
    'nback_0a': '0-back A',
    'nback_1a': '1-back',
    'nback_0b': '0-back B',
    'nback_2a': '2-back',
}

BAND_ORDER = ['delta', 'theta', 'alpha', 'beta', 'gamma']
BAND_LABELS = {
    'delta': 'Delta (0.5–4 Hz)',
    'theta': 'Theta (4–8 Hz)',
    'alpha': 'Alpha (8–13 Hz)',
    'beta': 'Beta (13–30 Hz)',
    'gamma': 'Gamma (30–50 Hz)',
}

# Figure settings
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})


# =============================================================================
# DATA LOADING
# =============================================================================

def load_bandpower(data_dir: Path, roi: str) -> pd.DataFrame:
    """Load bandpower CSV and compute subject-level means."""
    csv_file = data_dir / f'bandpower_{roi}.csv'
    if not csv_file.exists():
        print(f"  ⚠ {csv_file.name} not found")
        return pd.DataFrame()
    
    df = pd.read_csv(csv_file)
    
    # Standardize column names
    for col in ['SubjectID', 'Group', 'Channel', 'Task']:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    
    # Subject-level means (average across channels)
    band_cols = [c for c in df.columns if c.endswith(('_abs', '_rel', '_log'))]
    subj_means = df.groupby(['SubjectID', 'Group', 'Task'])[band_cols].mean().reset_index()
    
    return subj_means


def load_spectral_slope(data_dir: Path, roi: str) -> pd.DataFrame:
    """Load spectral slope CSV."""
    csv_file = data_dir / f'spectral_slope_{roi}.csv'
    if not csv_file.exists():
        print(f"  ⚠ {csv_file.name} not found")
        return pd.DataFrame()
    
    df = pd.read_csv(csv_file)
    
    for col in ['SubjectID', 'Group', 'Channel', 'Task']:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    
    # Subject-level means
    subj_means = df.groupby(['SubjectID', 'Group', 'Task'])[['Slope', 'Offset', 'R2']].mean().reset_index()
    
    return subj_means


def load_ersp_summary(data_dir: Path, roi: str) -> pd.DataFrame:
    """Load ERSP summary CSV."""
    csv_file = data_dir / f'ersp_{roi}_summary.csv'
    if not csv_file.exists():
        print(f"  ⚠ {csv_file.name} not found")
        return pd.DataFrame()
    
    df = pd.read_csv(csv_file)
    
    for col in ['SubjectID', 'Group', 'Task', 'Band', 'TimeWindow']:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    
    return df


def compute_ratios(bp_data: pd.DataFrame) -> pd.DataFrame:
    """Compute spectral ratios from bandpower data."""
    if bp_data.empty:
        return pd.DataFrame()
    
    rows = []
    for (subj, group, task), gdf in bp_data.groupby(['SubjectID', 'Group', 'Task']):
        row = {'SubjectID': subj, 'Group': group, 'Task': task}
        
        # Theta/Beta ratio
        if 'theta_abs' in gdf.columns and 'beta_abs' in gdf.columns:
            theta = gdf['theta_abs'].values[0]
            beta = gdf['beta_abs'].values[0]
            if beta > 0:
                row['Theta_Beta'] = theta / beta
        
        # Delta/Alpha ratio
        if 'delta_abs' in gdf.columns and 'alpha_abs' in gdf.columns:
            delta = gdf['delta_abs'].values[0]
            alpha = gdf['alpha_abs'].values[0]
            if alpha > 0:
                row['Delta_Alpha'] = delta / alpha
        
        rows.append(row)
    
    return pd.DataFrame(rows)


# =============================================================================
# STATISTICS
# =============================================================================

def compute_cohens_d(g1: np.ndarray, g2: np.ndarray) -> float:
    """Compute Cohen's d effect size."""
    g1, g2 = np.array(g1), np.array(g2)
    g1, g2 = g1[~np.isnan(g1)], g2[~np.isnan(g2)]
    
    if len(g1) < 2 or len(g2) < 2:
        return np.nan
    
    n1, n2 = len(g1), len(g2)
    pooled_std = np.sqrt(((n1-1)*np.var(g1, ddof=1) + (n2-1)*np.var(g2, ddof=1)) / (n1+n2-2))
    
    if pooled_std == 0:
        return np.nan
    
    return (np.mean(g1) - np.mean(g2)) / pooled_std


def run_group_comparison(data: pd.DataFrame, value_col: str, 
                         groups: list = None) -> dict:
    """
    Run omnibus + pairwise group comparisons.
    Uses LME if available and data has repeated measures, else ANOVA.
    """
    if groups is None:
        groups = GROUP_ORDER
    
    # Filter to valid groups
    data = data[data['Group'].isin(groups)].copy()
    data = data.dropna(subset=[value_col])
    
    if len(data) < 6:
        return {'Omnibus_F': np.nan, 'Omnibus_p': np.nan}
    
    result = {}
    
    # Descriptives per group
    for g in groups:
        gdata = data[data['Group'] == g][value_col]
        n_subj = data[data['Group'] == g]['SubjectID'].nunique()
        result[f'{g}_mean'] = gdata.mean()
        result[f'{g}_std'] = gdata.std()
        result[f'{g}_n'] = len(gdata)
        result[f'{g}_nSubj'] = n_subj
    
    # Omnibus test
    group_data = [data[data['Group'] == g][value_col].values for g in groups]
    group_data = [g[~np.isnan(g)] for g in group_data]
    
    if all(len(g) >= 2 for g in group_data):
        try:
            f_stat, p_val = sp_stats.f_oneway(*group_data)
            result['Omnibus_F'] = f_stat
            result['Omnibus_p'] = p_val
        except:
            result['Omnibus_F'] = np.nan
            result['Omnibus_p'] = np.nan
    else:
        result['Omnibus_F'] = np.nan
        result['Omnibus_p'] = np.nan
    
    # Pairwise comparisons
    for g1, g2 in combinations(groups, 2):
        pair_name = f'{g1}_vs_{g2}'
        d1 = data[data['Group'] == g1][value_col].dropna().values
        d2 = data[data['Group'] == g2][value_col].dropna().values
        
        if len(d1) >= 2 and len(d2) >= 2:
            t_stat, p_val = sp_stats.ttest_ind(d1, d2)
            result[f'{pair_name}_t'] = t_stat
            result[f'{pair_name}_p'] = p_val
            result[f'{pair_name}_d'] = compute_cohens_d(d1, d2)
        else:
            result[f'{pair_name}_t'] = np.nan
            result[f'{pair_name}_p'] = np.nan
            result[f'{pair_name}_d'] = np.nan
    
    return result


def apply_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    p = np.array(p_values)
    n = len(p)
    
    valid_mask = ~np.isnan(p)
    if not np.any(valid_mask):
        return p
    
    p_valid = p[valid_mask]
    sorted_idx = np.argsort(p_valid)
    sorted_p = p_valid[sorted_idx]
    
    # BH correction
    q = sorted_p * n / (np.arange(n - (~valid_mask).sum()) + 1)
    
    # Ensure monotonicity
    for i in range(len(q) - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    q = np.minimum(q, 1.0)
    
    # Unsort
    q_unsorted = np.zeros_like(p_valid)
    q_unsorted[sorted_idx] = q
    
    result = np.full_like(p, np.nan)
    result[valid_mask] = q_unsorted
    
    return result


def run_stats_table(data: pd.DataFrame, value_col: str, 
                    groupby_cols: list, label_cols: dict = None) -> pd.DataFrame:
    """
    Run stats for each combination of groupby_cols.
    Returns a DataFrame with one row per combination.
    """
    if data.empty:
        return pd.DataFrame()
    
    results = []
    for keys, gdf in data.groupby(groupby_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        
        row = dict(zip(groupby_cols, keys))
        stats = run_group_comparison(gdf, value_col)
        row.update(stats)
        results.append(row)
    
    if not results:
        return pd.DataFrame()
    
    df = pd.DataFrame(results)
    
    # FDR correction on all p-value columns
    p_cols = [c for c in df.columns if c.endswith('_p')]
    all_p = df[p_cols].values.flatten()
    all_q = apply_fdr(all_p)
    
    # Insert q columns after p columns
    for i, p_col in enumerate(p_cols):
        q_col = p_col.replace('_p', '_q')
        col_idx = df.columns.get_loc(p_col)
        q_values = all_q[i::len(p_cols)]
        df.insert(col_idx + 1, q_col, q_values[:len(df)])
    
    return df


# =============================================================================
# PLOTTING
# =============================================================================

def plot_grouped_bars(data: pd.DataFrame, value_col: str, stats_df: pd.DataFrame,
                      tasks: list, ylabel: str, title: str, 
                      output_path: Path, stats_filter: dict = None,
                      figsize: tuple = (10, 5)):
    """
    Plot grouped bar chart with significance annotations.
    """
    if data.empty:
        return
    
    # Filter data to valid tasks and groups
    plot_data = data[data['Task'].isin(tasks) & data['Group'].isin(GROUP_ORDER)].copy()
    
    if plot_data.empty:
        return
    
    # Compute subject-level means for plotting
    subj_means = plot_data.groupby(['SubjectID', 'Group', 'Task'])[value_col].mean().reset_index()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    n_tasks = len(tasks)
    n_groups = len(GROUP_ORDER)
    bar_width = 0.25
    x = np.arange(n_tasks)
    
    for i, group in enumerate(GROUP_ORDER):
        group_data = subj_means[subj_means['Group'] == group]
        
        means = []
        sems = []
        for task in tasks:
            task_data = group_data[group_data['Task'] == task][value_col]
            means.append(task_data.mean() if len(task_data) > 0 else 0)
            sems.append(task_data.sem() if len(task_data) > 1 else 0)
        
        offset = (i - n_groups/2 + 0.5) * bar_width
        bars = ax.bar(x + offset, means, bar_width, 
                      label=GROUP_LABELS.get(group, group),
                      color=GROUP_COLORS.get(group, 'gray'),
                      yerr=sems, capsize=3, error_kw={'linewidth': 1})
    
    # Add significance markers
    if not stats_df.empty and stats_filter:
        for t_idx, task in enumerate(tasks):
            # Find stats row for this task
            mask = pd.Series([True] * len(stats_df))
            for col, val in stats_filter.items():
                if col in stats_df.columns:
                    mask &= (stats_df[col] == val)
            mask &= (stats_df['Task'] == task)
            
            if mask.sum() == 1:
                row = stats_df[mask].iloc[0]
                
                # Check omnibus significance
                if 'Omnibus_q' in row and row['Omnibus_q'] < 0.05:
                    y_max = subj_means[subj_means['Task'] == task][value_col].max()
                    ax.annotate('*', (t_idx, y_max * 1.1), 
                               ha='center', fontsize=14, fontweight='bold')
    
    ax.set_xticks(x)
    ax.set_xticklabels([TASK_LABELS.get(t, t) for t in tasks])
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight='bold')
    ax.legend(loc='upper right', frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_path.with_suffix('.png'))
    plt.savefig(output_path.with_suffix('.svg'))
    plt.close()


def plot_spectral_slope(data: pd.DataFrame, stats_df: pd.DataFrame,
                        output_dir: Path):
    """Plot spectral slope comparison across groups and tasks."""
    if data.empty:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Slope by task
    ax = axes[0]
    plot_data = data[data['Group'].isin(GROUP_ORDER) & data['Task'].isin(TASK_ORDER)]
    
    # Create grouped positions
    tasks = [t for t in TASK_ORDER if t in plot_data['Task'].values]
    n_tasks = len(tasks)
    n_groups = len(GROUP_ORDER)
    bar_width = 0.25
    x = np.arange(n_tasks)
    
    for i, group in enumerate(GROUP_ORDER):
        gdata = plot_data[plot_data['Group'] == group]
        means = [gdata[gdata['Task'] == t]['Slope'].mean() for t in tasks]
        sems = [gdata[gdata['Task'] == t]['Slope'].sem() for t in tasks]
        
        offset = (i - n_groups/2 + 0.5) * bar_width
        ax.bar(x + offset, means, bar_width,
               label=GROUP_LABELS.get(group, group),
               color=GROUP_COLORS.get(group, 'gray'),
               yerr=sems, capsize=3)
    
    ax.set_xticks(x)
    ax.set_xticklabels([TASK_LABELS.get(t, t) for t in tasks], rotation=45, ha='right')
    ax.set_ylabel('Spectral Slope (1/f exponent)')
    ax.set_title('Aperiodic Slope by Task', fontweight='bold')
    ax.legend(loc='upper right', frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Right: Grand average slope (outliers excluded via IQR)
    ax = axes[1]
    grand_avg = data.groupby(['SubjectID', 'Group'])['Slope'].mean().reset_index()

    for i, group in enumerate(GROUP_ORDER):
        gdata = grand_avg[grand_avg['Group'] == group]['Slope'].dropna().values

        # IQR outlier detection
        q1, q3 = np.percentile(gdata, [25, 75])
        iqr = q3 - q1
        lower, upper = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        outlier_mask = (gdata < lower) | (gdata > upper)
        normal_vals = gdata[~outlier_mask]
        outlier_vals = gdata[outlier_mask]

        if len(outlier_vals) > 0:
            print(f"    [Grand Average Slope] {group}: {len(outlier_vals)} outlier(s) excluded "
                  f"(values: {outlier_vals.round(3).tolist()})")

        # Violin uses only clean values
        if len(normal_vals) >= 2:
            parts = ax.violinplot([normal_vals], positions=[i],
                                  showmeans=True, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(GROUP_COLORS.get(group, 'gray'))
                pc.set_alpha(0.6)

        # Normal points in group colour
        jitter = np.random.default_rng(42 + i).uniform(-0.1, 0.1, len(normal_vals))
        ax.scatter(i + jitter, normal_vals,
                   color=GROUP_COLORS.get(group, 'gray'),
                   alpha=0.7, s=30, edgecolor='white', linewidth=0.5)

        # Outlier points in grey
        if len(outlier_vals) > 0:
            jitter_out = np.random.default_rng(99 + i).uniform(-0.1, 0.1, len(outlier_vals))
            ax.scatter(i + jitter_out, outlier_vals,
                       color='lightgray', edgecolor='darkgray',
                       alpha=0.7, s=30, linewidth=0.5)
    
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels([GROUP_LABELS.get(g, g) for g in GROUP_ORDER])
    ax.set_ylabel('Spectral Slope (1/f exponent)')
    ax.set_title('Grand Average Slope', fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add stats annotation
    if not stats_df.empty and 'grand_average' in stats_df['Task'].values:
        row = stats_df[stats_df['Task'] == 'grand_average'].iloc[0]
        if 'Omnibus_p' in row:
            p_text = f"F = {row.get('Omnibus_F', np.nan):.2f}, p = {row['Omnibus_p']:.3f}"
            ax.text(0.5, 0.95, p_text, transform=ax.transAxes, 
                   ha='center', fontsize=9, style='italic')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'spectral_slope.png')
    plt.savefig(output_dir / 'spectral_slope.svg')
    plt.close()


def plot_ersp_summary(data: pd.DataFrame, stats_df: pd.DataFrame,
                      output_dir: Path):
    """Plot ERSP band power summaries."""
    if data.empty:
        return
    
    # Filter to valid data
    plot_data = data[data['Group'].isin(GROUP_ORDER)].copy()
    
    # Plot: Band × TimeWindow heatmap-style for each group
    bands = [b for b in BAND_ORDER if b in plot_data['Band'].values]
    time_windows = ['early', 'mid', 'late']
    tasks = [t for t in ['nback_0a', 'nback_1a', 'nback_0b', 'nback_2a'] 
             if t in plot_data['Task'].values]
    
    if not bands or not tasks:
        return
    
    fig, axes = plt.subplots(1, len(tasks), figsize=(4*len(tasks), 5), 
                              squeeze=False)
    axes = axes.flatten()
    
    for t_idx, task in enumerate(tasks):
        ax = axes[t_idx]
        task_data = plot_data[plot_data['Task'] == task]
        
        # Create matrix: groups × bands (averaging over time windows)
        matrix = np.zeros((len(GROUP_ORDER), len(bands)))
        
        for g_idx, group in enumerate(GROUP_ORDER):
            for b_idx, band in enumerate(bands):
                vals = task_data[(task_data['Group'] == group) & 
                                 (task_data['Band'] == band)]['MeanPower_dB']
                matrix[g_idx, b_idx] = vals.mean() if len(vals) > 0 else np.nan
        
        im = ax.imshow(matrix, aspect='auto', cmap='RdBu_r', 
                       vmin=-3, vmax=3)
        
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels([BAND_LABELS.get(b, b).split()[0] for b in bands], 
                          rotation=45, ha='right')
        ax.set_yticks(range(len(GROUP_ORDER)))
        ax.set_yticklabels([GROUP_LABELS.get(g, g) for g in GROUP_ORDER])
        ax.set_title(TASK_LABELS.get(task, task), fontweight='bold')
        
        # Add values
        for g_idx in range(len(GROUP_ORDER)):
            for b_idx in range(len(bands)):
                val = matrix[g_idx, b_idx]
                if not np.isnan(val):
                    color = 'white' if abs(val) > 1.5 else 'black'
                    ax.text(b_idx, g_idx, f'{val:.1f}', ha='center', va='center',
                           fontsize=8, color=color)
    
    # Colorbar
    cbar = fig.colorbar(im, ax=axes, shrink=0.6, label='Power (dB)')
    
    fig.suptitle('ERSP Band Power by Task', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'ersp_summary.png')
    plt.savefig(output_dir / 'ersp_summary.svg')
    plt.close()


def plot_ratios(data: pd.DataFrame, stats_df: pd.DataFrame, 
                output_dir: Path):
    """Plot spectral ratios (Theta/Beta, Delta/Alpha)."""
    if data.empty:
        return
    
    ratios = ['Theta_Beta', 'Delta_Alpha']
    ratio_labels = {'Theta_Beta': 'Theta / Beta', 'Delta_Alpha': 'Delta / Alpha'}
    
    for ratio in ratios:
        if ratio not in data.columns:
            continue
        
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        
        # Left: By task
        ax = axes[0]
        plot_data = data[data['Group'].isin(GROUP_ORDER) & data['Task'].isin(TASK_ORDER)]
        tasks = [t for t in TASK_ORDER if t in plot_data['Task'].values]
        
        n_tasks = len(tasks)
        bar_width = 0.25
        x = np.arange(n_tasks)
        
        for i, group in enumerate(GROUP_ORDER):
            gdata = plot_data[plot_data['Group'] == group]
            means = [gdata[gdata['Task'] == t][ratio].mean() for t in tasks]
            sems = [gdata[gdata['Task'] == t][ratio].sem() for t in tasks]
            
            offset = (i - len(GROUP_ORDER)/2 + 0.5) * bar_width
            ax.bar(x + offset, means, bar_width,
                   label=GROUP_LABELS.get(group, group),
                   color=GROUP_COLORS.get(group, 'gray'),
                   yerr=sems, capsize=3)
        
        ax.set_xticks(x)
        ax.set_xticklabels([TASK_LABELS.get(t, t) for t in tasks], rotation=45, ha='right')
        ax.set_ylabel(ratio_labels.get(ratio, ratio))
        ax.set_title(f'{ratio_labels.get(ratio, ratio)} by Task', fontweight='bold')
        ax.legend(loc='upper right', frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Right: Grand average
        ax = axes[1]
        grand_avg = data.groupby(['SubjectID', 'Group'])[ratio].mean().reset_index()
        
        positions = []
        for i, group in enumerate(GROUP_ORDER):
            gdata = grand_avg[grand_avg['Group'] == group][ratio].dropna()
            if len(gdata) == 0:
                continue
            
            positions.append(i)
            
            # Box plot
            bp = ax.boxplot([gdata.values], positions=[i], widths=0.5,
                           patch_artist=True, showfliers=False)
            bp['boxes'][0].set_facecolor(GROUP_COLORS.get(group, 'gray'))
            bp['boxes'][0].set_alpha(0.6)
            
            # Individual points
            jitter = np.random.uniform(-0.15, 0.15, len(gdata))
            ax.scatter(np.full(len(gdata), i) + jitter, gdata,
                      color=GROUP_COLORS.get(group, 'gray'),
                      alpha=0.7, s=40, edgecolor='white', linewidth=0.5, zorder=3)
        
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels([GROUP_LABELS.get(g, g) for g in GROUP_ORDER])
        ax.set_ylabel(ratio_labels.get(ratio, ratio))
        ax.set_title('Grand Average', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(output_dir / f'ratio_{ratio.lower()}.png')
        plt.savefig(output_dir / f'ratio_{ratio.lower()}.svg')
        plt.close()


def plot_bandpower_summary(data: pd.DataFrame, stats_df: pd.DataFrame,
                           output_dir: Path):
    """Plot band power for each band."""
    if data.empty:
        return
    
    for band in BAND_ORDER:
        rel_col = f'{band}_rel'
        log_col = f'{band}_log'
        
        if rel_col not in data.columns:
            continue
        
        # Relative power plot
        plot_grouped_bars(
            data=data,
            value_col=rel_col,
            stats_df=stats_df[stats_df.get('Band', pd.Series()) == band] if 'Band' in stats_df.columns else pd.DataFrame(),
            tasks=TASK_ORDER,
            ylabel='Relative Power',
            title=f'Relative {BAND_LABELS.get(band, band)}',
            output_path=output_dir / f'bandpower_rel_{band}',
            stats_filter={'Band': band} if 'Band' in stats_df.columns else {}
        )
        
        # Log power plot  
        if log_col in data.columns:
            plot_grouped_bars(
                data=data,
                value_col=log_col,
                stats_df=stats_df[stats_df.get('Band', pd.Series()) == band] if 'Band' in stats_df.columns else pd.DataFrame(),
                tasks=TASK_ORDER,
                ylabel='Log Power (10·log₁₀ µV²)',
                title=f'Log {BAND_LABELS.get(band, band)}',
                output_path=output_dir / f'bandpower_log_{band}',
                stats_filter={'Band': band} if 'Band' in stats_df.columns else {}
            )


# =============================================================================
# MAIN
# =============================================================================

def main():
    if len(sys.argv) < 3:
        print("Usage: python EEG_statistics.py <data_dir> <output_dir> [roi]")
        print("  data_dir   - Directory with CSV outputs from EEG_analysis.m")
        print("  output_dir - Directory for statistics and figures")
        print("  roi        - ROI suffix (default: 'allchan')")
        sys.exit(1)
    
    data_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    roi = sys.argv[3] if len(sys.argv) > 3 else 'allchan'
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print("EEG STATISTICAL ANALYSIS")
    print(f"{'='*60}")
    print(f"Data:   {data_dir}")
    print(f"Output: {output_dir}")
    print(f"ROI:    {roi}")
    
    # =========================================================================
    # Load data
    # =========================================================================
    print("\n→ Loading data...")
    
    bp_data = load_bandpower(data_dir, roi)
    slope_data = load_spectral_slope(data_dir, roi)
    ersp_data = load_ersp_summary(data_dir, roi)
    
    if not bp_data.empty:
        print(f"  Bandpower: {bp_data['SubjectID'].nunique()} subjects")
        ratio_data = compute_ratios(bp_data)
    else:
        ratio_data = pd.DataFrame()
    
    if not slope_data.empty:
        print(f"  Spectral slope: {slope_data['SubjectID'].nunique()} subjects")
    
    if not ersp_data.empty:
        print(f"  ERSP: {ersp_data['SubjectID'].nunique()} subjects")
    
    # =========================================================================
    # Run statistics
    # =========================================================================
    print("\n→ Running statistics...")
    
    # --- Band Power Stats ---
    if not bp_data.empty:
        print("  Band power...")
        
        # Reshape to long format for stats
        bp_long_rows = []
        for band in BAND_ORDER:
            for suffix in ['rel', 'log']:
                col = f'{band}_{suffix}'
                if col in bp_data.columns:
                    temp = bp_data[['SubjectID', 'Group', 'Task', col]].copy()
                    temp['Band'] = band
                    temp['Metric'] = suffix
                    temp = temp.rename(columns={col: 'Value'})
                    bp_long_rows.append(temp)
        
        if bp_long_rows:
            bp_long = pd.concat(bp_long_rows, ignore_index=True)
            
            stats_bp = []
            for (band, metric), gdf in bp_long.groupby(['Band', 'Metric']):
                task_stats = run_stats_table(gdf, 'Value', ['Task'])
                task_stats['Band'] = band
                task_stats['Metric'] = metric
                stats_bp.append(task_stats)
            
            if stats_bp:
                stats_bp = pd.concat(stats_bp, ignore_index=True)
                stats_bp.to_csv(output_dir / 'stats_bandpower.csv', index=False)
                print(f"    ✓ stats_bandpower.csv ({len(stats_bp)} rows)")
    else:
        stats_bp = pd.DataFrame()
    
    # --- Ratio Stats ---
    if not ratio_data.empty:
        print("  Spectral ratios...")
        
        ratio_long_rows = []
        for ratio in ['Theta_Beta', 'Delta_Alpha']:
            if ratio in ratio_data.columns:
                temp = ratio_data[['SubjectID', 'Group', 'Task', ratio]].copy()
                temp['Ratio'] = ratio
                temp = temp.rename(columns={ratio: 'Value'})
                ratio_long_rows.append(temp)
        
        if ratio_long_rows:
            ratio_long = pd.concat(ratio_long_rows, ignore_index=True)
            
            stats_ratio = []
            for ratio_name, gdf in ratio_long.groupby('Ratio'):
                task_stats = run_stats_table(gdf, 'Value', ['Task'])
                task_stats['Ratio'] = ratio_name
                stats_ratio.append(task_stats)
                
                # Grand average
                grand = gdf.groupby(['SubjectID', 'Group'])['Value'].mean().reset_index()
                grand['Task'] = 'grand_average'
                grand_stats = run_group_comparison(grand, 'Value')
                grand_stats['Task'] = 'grand_average'
                grand_stats['Ratio'] = ratio_name
                stats_ratio.append(pd.DataFrame([grand_stats]))
            
            if stats_ratio:
                stats_ratio = pd.concat(stats_ratio, ignore_index=True)
                stats_ratio.to_csv(output_dir / 'stats_ratios.csv', index=False)
                print(f"    ✓ stats_ratios.csv ({len(stats_ratio)} rows)")
    else:
        stats_ratio = pd.DataFrame()
    
    # --- Spectral Slope Stats ---
    if not slope_data.empty:
        print("  Spectral slope...")
        
        stats_slope = run_stats_table(slope_data, 'Slope', ['Task'])
        
        # Grand average
        grand = slope_data.groupby(['SubjectID', 'Group'])['Slope'].mean().reset_index()
        grand['Task'] = 'grand_average'
        grand_stats = run_group_comparison(grand, 'Slope')
        grand_stats['Task'] = 'grand_average'
        stats_slope = pd.concat([stats_slope, pd.DataFrame([grand_stats])], ignore_index=True)
        
        stats_slope.to_csv(output_dir / 'stats_slope.csv', index=False)
        print(f"    ✓ stats_slope.csv ({len(stats_slope)} rows)")
    else:
        stats_slope = pd.DataFrame()
    
    # --- ERSP Stats ---
    if not ersp_data.empty:
        print("  ERSP...")
        
        stats_ersp = run_stats_table(ersp_data, 'MeanPower_dB', 
                                      ['Task', 'Band', 'TimeWindow'])
        stats_ersp.to_csv(output_dir / 'stats_ersp.csv', index=False)
        print(f"    ✓ stats_ersp.csv ({len(stats_ersp)} rows)")
    else:
        stats_ersp = pd.DataFrame()
    
    # =========================================================================
    # Generate plots
    # =========================================================================
    print("\n→ Generating figures...")
    
    # Band power plots
    if not bp_data.empty:
        print("  Band power...")
        plot_bandpower_summary(bp_data, stats_bp, output_dir)
        print("    ✓ bandpower_*.png/svg")
    
    # Ratio plots
    if not ratio_data.empty:
        print("  Spectral ratios...")
        plot_ratios(ratio_data, stats_ratio, output_dir)
        print("    ✓ ratio_*.png/svg")
    
    # Spectral slope plots
    if not slope_data.empty:
        print("  Spectral slope...")
        plot_spectral_slope(slope_data, stats_slope, output_dir)
        print("    ✓ spectral_slope.png/svg")
    
    # ERSP plots
    if not ersp_data.empty:
        print("  ERSP...")
        plot_ersp_summary(ersp_data, stats_ersp, output_dir)
        print("    ✓ ersp_summary.png/svg")
    
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"Stats:   {output_dir}/stats_*.csv")
    print(f"Figures: {output_dir}/*.png, *.svg")


if __name__ == '__main__':
    main()
