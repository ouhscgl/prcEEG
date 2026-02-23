#!/usr/bin/env python3
"""
EEG Statistical Visualization
==============================
Plots group comparisons for band power, spectral ratios, and ERP metrics.
Reads .mat files from EEG_analysis and stats CSVs from EEG_statistics.

Usage:
    python EEG_plotStatistics.py <data_dir> <stats_dir> <output_dir> [roi]

    data_dir  - Directory with bandpower_<roi>.mat and erp_<roi>.mat
    stats_dir - Directory with stats_bandpower_*.csv, stats_ratios.csv, stats_erp.csv
    output_dir- Directory for output figures
    roi       - ROI suffix (default: 'pfc')
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from itertools import combinations
import warnings
import sys

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

TASK_ORDER = ['eyes_open', 'nback_0a', 'nback_1a', 'nback_0b', 'nback_2a']
TASK_LABELS = {
    'eyes_open': 'Eyes Open',
    'nback_0a': '0-back A',
    'nback_1a': '1-back',
    'nback_0b': '0-back B',
    'nback_2a': '2-back',
}

GROUP_COLORS = {
    'healthy_controls': (0.467, 0.867, 0.467),
    'normal_performers': (0.886, 0.698, 0.447),
    'low_performers': (0.753, 0.502, 0.494),
}

BAND_ORDER = ['delta', 'theta', 'alpha', 'beta', 'gamma']
BAND_LABELS = {
    'delta': 'Delta (0.5–4 Hz)',
    'theta': 'Theta (4–8 Hz)',
    'alpha': 'Alpha (8–13 Hz)',
    'beta': 'Beta (13–30 Hz)',
    'gamma': 'Gamma (30–50 Hz)',
}

RATIO_ORDER = ['Theta/Beta', 'Delta/Alpha']
RATIO_LABELS = {
    'Theta/Beta': 'Theta / Beta Ratio',
    'Delta/Alpha': 'Delta / Alpha Ratio',
}

ERP_METRIC_ORDER = ['PeakAmp', 'MeanAmp', 'PeakLat', 'CentLat']
ERP_METRIC_LABELS = {
    'PeakAmp': 'P300 Peak Amplitude (µV)',
    'MeanAmp': 'P300 Mean Amplitude (µV)',
    'PeakLat': 'P300 Peak Latency (ms)',
    'CentLat': 'P300 Centroid Latency (ms)',
}

# ERP tasks (no eyes_open)
ERP_TASK_ORDER = ['nback_0a', 'nback_1a', 'nback_0b', 'nback_2a']


# =============================================================================
# DATA LOADING
# =============================================================================

def load_bandpower_subj_means(data_dir: Path, roi: str) -> pd.DataFrame:
    """
    Load bandpower CSV → compute subject-level means (mean across channels).
    Returns long-format: SubjectID, Group, Task, Band, Value_abs, Value_rel, Value_log
    """
    csv_file = data_dir / f'bandpower_{roi}.csv'
    if not csv_file.exists():
        print(f"  Warning: {csv_file} not found")
        return pd.DataFrame()

    bp = pd.read_csv(csv_file)

    for col in ['SubjectID', 'Group', 'Channel', 'Task']:
        if col in bp.columns:
            bp[col] = bp[col].astype(str).str.strip()

    rows = []
    for band in BAND_ORDER:
        abs_col = f'{band}_abs'
        rel_col = f'{band}_rel'
        log_col = f'{band}_log'

        if abs_col not in bp.columns:
            continue

        # Subject mean across channels
        grouped = bp.groupby(['SubjectID', 'Group', 'Task']).agg(
            Value_abs=(abs_col, 'mean'),
            Value_rel=(rel_col, 'mean'),
            Value_log=(log_col, 'mean'),
        ).reset_index()

        grouped['Band'] = band
        rows.append(grouped)

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, ignore_index=True)


def load_erp_subj_means(data_dir: Path, roi: str) -> pd.DataFrame:
    """
    Load ERP CSV → compute subject-level means (correct trials only).
    Returns: SubjectID, Group, Task, PeakAmp, PeakLat, MeanAmp, CentLat
    """
    csv_file = data_dir / f'erp_{roi}.csv'
    if not csv_file.exists():
        print(f"  Warning: {csv_file} not found")
        return pd.DataFrame()

    erp = pd.read_csv(csv_file)

    for col in ['SubjectID', 'Group', 'Task', 'Accuracy']:
        if col in erp.columns:
            erp[col] = erp[col].astype(str).str.strip()

    # Filter to correct trials only
    #erp = erp[erp['Accuracy'] == 'correct']

    # Convert amplitude from V to µV if values are tiny
    for col in ['PeakAmp', 'MeanAmp']:
        if col in erp.columns:
            erp[col] = pd.to_numeric(erp[col], errors='coerce')
            if erp[col].abs().median() < 0.001:
                erp[col] = erp[col] * 1e6

    for col in ['PeakLat', 'CentLat']:
        if col in erp.columns:
            erp[col] = pd.to_numeric(erp[col], errors='coerce')

    # Subject means
    metric_cols = [c for c in ['PeakAmp', 'PeakLat', 'MeanAmp', 'CentLat'] if c in erp.columns]
    subj_means = erp.groupby(['SubjectID', 'Group', 'Task'])[metric_cols].mean().reset_index()

    return subj_means


def load_erp_waveforms(data_dir: Path, roi: str):
    """Load ERP waveform CSV for grand-average waveform plots."""
    csv_file = data_dir / f'waveforms_{roi}.csv'
    if not csv_file.exists():
        print(f"  ⚠ {csv_file} not found, skipping waveforms")
        return None, None, None

    df = pd.read_csv(csv_file)

    # Time-point columns are like t-200.0, t0.0, t800.0 (t followed by number or minus)
    import re
    time_cols = [c for c in df.columns if re.match(r'^t-?\d', c)]
    info_cols = [c for c in df.columns if c not in time_cols]

    info = df[info_cols].copy()
    for col in ['SubjectID', 'Group', 'Task', 'Accuracy']:
        if col in info.columns:
            info[col] = info[col].astype(str).str.strip()

    wf_matrix = df[time_cols].values
    times = np.array([float(c[1:]) for c in time_cols])  # strip leading 't'

    return wf_matrix, info, times


def load_stats_csv(stats_dir: Path, filename: str) -> pd.DataFrame:
    """Load a stats CSV, return empty DataFrame if not found."""
    filepath = stats_dir / filename
    if filepath.exists():
        return pd.read_csv(filepath)
    return pd.DataFrame()


def compute_ratios_subj_means(data_dir: Path, roi: str) -> pd.DataFrame:
    """Compute subject-level ratio means from per-channel bandpower."""
    csv_file = data_dir / f'bandpower_{roi}.csv'
    if not csv_file.exists():
        return pd.DataFrame()

    bp = pd.read_csv(csv_file)
    for col in ['SubjectID', 'Group', 'Channel', 'Task']:
        if col in bp.columns:
            bp[col] = bp[col].astype(str).str.strip()

    rows = []

    # Theta/Beta
    if 'theta_abs' in bp.columns and 'beta_abs' in bp.columns:
        bp['_tb_ratio'] = pd.to_numeric(bp['theta_abs'], errors='coerce') / \
                          pd.to_numeric(bp['beta_abs'], errors='coerce')
        grouped = bp.groupby(['SubjectID', 'Group', 'Task'])['_tb_ratio'].mean().reset_index()
        grouped.columns = ['SubjectID', 'Group', 'Task', 'Value']
        grouped['Ratio'] = 'Theta/Beta'
        rows.append(grouped)

    # Delta/Alpha
    if 'delta_abs' in bp.columns and 'alpha_abs' in bp.columns:
        bp['_da_ratio'] = pd.to_numeric(bp['delta_abs'], errors='coerce') / \
                          pd.to_numeric(bp['alpha_abs'], errors='coerce')
        grouped = bp.groupby(['SubjectID', 'Group', 'Task'])['_da_ratio'].mean().reset_index()
        grouped.columns = ['SubjectID', 'Group', 'Task', 'Value']
        grouped['Ratio'] = 'Delta/Alpha'
        rows.append(grouped)

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, ignore_index=True)


# =============================================================================
# STATISTICS HELPERS
# =============================================================================

def get_pairwise_p(stats_df: pd.DataFrame, task: str, group1: str, group2: str,
                   task_col='Task') -> float:
    """Extract pairwise p-value from stats CSV for a given task and group pair."""
    if stats_df.empty:
        return np.nan

    row = stats_df[stats_df[task_col] == task]
    if row.empty:
        return np.nan

    pair_name = f'{group1}_vs_{group2}'
    p_col = f'{pair_name}_p'

    if p_col in row.columns:
        return row[p_col].values[0]

    # Try reversed order
    pair_name_rev = f'{group2}_vs_{group1}'
    p_col_rev = f'{pair_name_rev}_p'
    if p_col_rev in row.columns:
        return row[p_col_rev].values[0]

    return np.nan


def sig_stars(p: float) -> str:
    if np.isnan(p):
        return ''
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''


# =============================================================================
# PLOTTING
# =============================================================================

def add_significance_bracket(ax, x1, x2, y, p_val, height=0.02):
    """Add significance bracket between two x positions."""
    stars = sig_stars(p_val)
    if not stars:
        return y  # no bracket drawn, return same y

    y_bracket = y
    bracket_height = height * (ax.get_ylim()[1] - ax.get_ylim()[0])

    ax.plot([x1, x1, x2, x2],
            [y_bracket, y_bracket + bracket_height, y_bracket + bracket_height, y_bracket],
            'k-', linewidth=1)
    ax.text((x1 + x2) / 2, y_bracket + bracket_height, stars,
            ha='center', va='bottom', fontsize=10, fontweight='bold')

    return y_bracket + bracket_height * 2.5


def plot_grouped_bars(subj_data: pd.DataFrame, stats_df: pd.DataFrame,
                      value_col: str, filter_col: str, filter_val: str,
                      tasks: list, ylabel: str, title: str,
                      output_path: Path, stats_task_col='Task',
                      stats_filter_col=None, stats_filter_val=None,
                      show_individual_markers=True):
    """
    Generic grouped bar plot: tasks on x-axis, groups as bars.

    subj_data: DataFrame with SubjectID, Group, Task, and value_col
    stats_df:  Stats CSV for p-values
    filter_col/val: Optional column to filter subj_data (e.g., Band=='delta')
    show_individual_markers: If True, overlay individual subject data points
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    # Filter data
    df = subj_data.copy()
    if filter_col and filter_val:
        df = df[df[filter_col] == filter_val]

    # Filter stats
    sdf = stats_df.copy()
    if stats_filter_col and stats_filter_val and not sdf.empty:
        sdf = sdf[sdf[stats_filter_col] == stats_filter_val]

    # Filter to requested tasks
    df = df[df['Task'].isin(tasks)]

    n_tasks = len(tasks)
    n_groups = len(GROUP_ORDER)
    bar_width = 0.25
    x_base = np.arange(n_tasks)

    for i, group in enumerate(GROUP_ORDER):
        gdf = df[df['Group'] == group]

        means = []
        sems = []
        for task in tasks:
            vals = gdf[gdf['Task'] == task][value_col].dropna().values
            if len(vals) > 0:
                means.append(np.mean(vals))
                sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
            else:
                means.append(0)
                sems.append(0)

        x_pos = x_base + (i - 1) * bar_width

        ax.bar(x_pos, means, bar_width,
               label=GROUP_LABELS[group],
               color=GROUP_COLORS[group],
               edgecolor='white', linewidth=1)
        ax.errorbar(x_pos, means, yerr=sems, fmt='none',
                    color='black', capsize=3, linewidth=1)

        # Individual subject points with jitter (optional)
        if show_individual_markers:
            for t_idx, task in enumerate(tasks):
                vals = gdf[gdf['Task'] == task][value_col].dropna().values
                if len(vals) > 0:
                    jitter = np.random.default_rng(42 + i + t_idx).uniform(-0.06, 0.06, len(vals))
                    ax.scatter(x_pos[t_idx] + jitter, vals,
                               color=np.array(GROUP_COLORS[group]) * 0.7,
                               edgecolor='black', s=30, alpha=0.6, zorder=3, linewidths=0.5)

    # Significance brackets
    if not sdf.empty:
        for t_idx, task in enumerate(tasks):
            # Get max bar height for this task cluster
            all_vals = df[df['Task'] == task][value_col].dropna()
            if all_vals.empty:
                continue
            max_y = all_vals.max() * 1.05 if all_vals.max() > 0 else all_vals.mean() + all_vals.std() * 2

            bracket_y = max_y
            for g1, g2 in combinations(range(n_groups), 2):
                p = get_pairwise_p(sdf, task, GROUP_ORDER[g1], GROUP_ORDER[g2],
                                   task_col=stats_task_col)
                if not np.isnan(p) and p < 0.05:
                    x1 = t_idx + (g1 - 1) * bar_width
                    x2 = t_idx + (g2 - 1) * bar_width
                    bracket_y = add_significance_bracket(ax, x1, x2, bracket_y, p)

    # Formatting
    ax.set_xticks(x_base)
    ax.set_xticklabels([TASK_LABELS.get(t, t) for t in tasks], fontsize=11)
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best', frameon=False, fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.15)

    plt.tight_layout()
    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close()
    print(f"  ✓ {output_path.stem}")


def plot_ratio_grand_average(subj_data: pd.DataFrame, stats_df: pd.DataFrame,
                              ratio_name: str, output_path: Path,
                              show_individual_markers=False):
    """
    Plot grand average (across all tasks) for a spectral ratio.
    One bar per group, comparing overall ratio levels.
    
    subj_data: DataFrame with SubjectID, Group, Task, Value, Ratio columns
    stats_df: Stats CSV (unused for now, but could add grand avg stats)
    ratio_name: 'Theta/Beta' or 'Delta/Alpha'
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Filter to this ratio
    df = subj_data[subj_data['Ratio'] == ratio_name].copy()
    
    if df.empty:
        plt.close()
        return
    
    # Compute subject-level grand average (mean across all tasks)
    grand_avg = df.groupby(['SubjectID', 'Group'])['Value'].mean().reset_index()
    
    # Calculate group statistics with IQR-based outlier exclusion
    group_stats = []
    for group in GROUP_ORDER:
        values = grand_avg[grand_avg['Group'] == group]['Value'].dropna().values
        if len(values) == 0:
            continue

        # IQR outlier detection (same method as plot_delta_dar scatter)
        q1, q3 = np.percentile(values, [25, 75])
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        outlier_mask = (values < lower_bound) | (values > upper_bound)
        normal_vals = values[~outlier_mask]
        outlier_vals = values[outlier_mask]

        if len(outlier_vals) > 0:
            print(f"    [{ratio_name} grand avg] {group}: {len(outlier_vals)} outlier(s) excluded "
                  f"(values: {outlier_vals.round(3).tolist()})")

        # Mean and SEM computed WITHOUT outliers
        clean_vals = normal_vals if len(normal_vals) > 0 else values
        group_stats.append({
            'Group': group,
            'Mean': clean_vals.mean(),
            'SEM': clean_vals.std(ddof=1) / np.sqrt(len(clean_vals)) if len(clean_vals) > 1 else 0,
            'values': normal_vals,      # non-outlier points (for markers)
            'outlier_vals': outlier_vals,
            'all_values': values,       # kept for pairwise stats
        })

    x_pos = np.arange(len(group_stats))

    # Plot bars
    for i, gs in enumerate(group_stats):
        ax.bar(i, gs['Mean'], 0.6,
               color=GROUP_COLORS[gs['Group']],
               edgecolor='white', linewidth=1)
        ax.errorbar(i, gs['Mean'], yerr=gs['SEM'], fmt='none',
                    color='black', capsize=4, linewidth=1.5)

        # Individual data points (optional)
        if show_individual_markers:
            # Normal points in group colour
            if len(gs['values']) > 0:
                jitter = np.random.default_rng(42 + i).uniform(-0.15, 0.15, len(gs['values']))
                ax.scatter(i + jitter, gs['values'],
                           color=np.array(GROUP_COLORS[gs['Group']]) * 0.7,
                           edgecolor='black', s=50, alpha=0.7, zorder=3, linewidths=0.5)
            # Outlier points in grey
            if len(gs['outlier_vals']) > 0:
                jitter_out = np.random.default_rng(99 + i).uniform(-0.15, 0.15, len(gs['outlier_vals']))
                ax.scatter(i + jitter_out, gs['outlier_vals'],
                           color='lightgray', edgecolor='darkgray',
                           s=50, alpha=0.7, zorder=2, linewidths=0.5,
                           label='_nolegend_')
    
    # Pairwise comparisons — intentionally uses all_values (including outliers)
    # so the p-values remain consistent with EEG_statistics.py.
    # If you exclude outliers in EEG_statistics.py, switch to gs['values'] here.
    from scipy import stats as sp_stats
    p_values = []
    comparisons = list(combinations(range(len(group_stats)), 2))
    
    for i, j in comparisons:
        vals1 = group_stats[i]['all_values']
        vals2 = group_stats[j]['all_values']
        if len(vals1) >= 2 and len(vals2) >= 2:
            _, p = sp_stats.ttest_ind(vals1, vals2)
            p_values.append((i, j, p))

    # Add significance brackets
    if p_values:
        max_y = max([gs['Mean'] + gs['SEM'] for gs in group_stats])
        bracket_y = max_y * 1.05
        for i, j, p in p_values:
            if p < 0.05:
                bracket_y = add_significance_bracket(ax, i, j, bracket_y, p)
    
    # Formatting
    ax.set_xticks(x_pos)
    ax.set_xticklabels([GROUP_LABELS[gs['Group']] for gs in group_stats], fontsize=11)
    ax.set_ylabel('Ratio', fontsize=12, fontweight='bold')
    ax.set_title(f'{RATIO_LABELS[ratio_name]} (Grand Average)', fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.2)
    
    plt.tight_layout()
    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close()
    print(f"  ✓ {output_path.stem}")


def plot_delta_dar(stats_dir: Path, output_dir: Path):
    """
    Plot Delta DAR (change in Delta/Alpha Ratio from Eyes Open to 2-back).
    
    Reads subject_delta_dar.csv from stats_dir.
    Creates two types of plots for each metric:
    - Bar plot (mean ± SEM, no individual markers)
    - Scatter plot (individual data points only, no bars)
    
    Metrics: absolute difference, percent change, log ratio
    
    Negative values = good (reduced slow/fast ratio during cognitive task)
    Positive values = poor cortical activation
    """
    from scipy import stats as sp_stats
    
    csv_file = stats_dir / 'subject_delta_dar.csv'
    if not csv_file.exists():
        print(f"  ⚠ {csv_file} not found, skipping Delta DAR plot")
        return
    
    df = pd.read_csv(csv_file)
    
    # Ensure Group column is clean
    df['Group'] = df['Group'].astype(str).str.strip()
    
    # Compute metrics if not present
    if 'DeltaDAR_pct' not in df.columns and 'DAR_2back' in df.columns and 'DAR_EO' in df.columns:
        df['DeltaDAR_pct'] = (df['DAR_2back'] - df['DAR_EO']) / df['DAR_EO'] * 100
    if 'DeltaDAR_log' not in df.columns and 'DAR_2back' in df.columns and 'DAR_EO' in df.columns:
        df['DeltaDAR_log'] = np.log10(df['DAR_2back'] / df['DAR_EO'])
    
    # Define the metrics to plot
    metrics = [
        {
            'col': 'DeltaDAR',
            'ylabel': 'ΔDAR (2-back − Eyes Open)',
            'title': 'Delta/Alpha Ratio Change\n(Absolute Difference)',
            'filename': 'delta_dar_absolute',
        },
        {
            'col': 'DeltaDAR_pct',
            'ylabel': 'ΔDAR (% change)',
            'title': 'Delta/Alpha Ratio Change\n(Percent Change)',
            'filename': 'delta_dar_percent',
        },
        {
            'col': 'DeltaDAR_log',
            'ylabel': 'ΔDAR (log₁₀ ratio)',
            'title': 'Delta/Alpha Ratio Change\n(Log Ratio)',
            'filename': 'delta_dar_log',
        },
    ]
    
    for metric in metrics:
        col = metric['col']
        
        if col not in df.columns:
            print(f"  ⚠ Column {col} not found, skipping")
            continue
        
        # Calculate group statistics
        group_stats = []
        for group in GROUP_ORDER:
            values = df[df['Group'] == group][col].dropna()
            values = values[np.isfinite(values)]
            if len(values) > 0:
                group_stats.append({
                    'Group': group,
                    'Mean': values.mean(),
                    'SEM': values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 else 0,
                    'values': values.values,
                })
        
        if not group_stats:
            print(f"  ⚠ No valid {col} data found")
            continue
        
        # Pairwise comparisons (computed once, used for both plots)
        p_values = []
        comparisons = list(combinations(range(len(group_stats)), 2))
        for i, j in comparisons:
            vals1 = group_stats[i]['values']
            vals2 = group_stats[j]['values']
            if len(vals1) >= 2 and len(vals2) >= 2:
                _, p = sp_stats.ttest_ind(vals1, vals2)
                p_values.append((i, j, p))
        
        x_pos = np.arange(len(group_stats))
        
        # =====================================================================
        # PLOT 1: Bar plot (no individual markers)
        # =====================================================================
        fig, ax = plt.subplots(figsize=(8, 6))
        
        for i, gs in enumerate(group_stats):
            ax.bar(i, gs['Mean'], 0.6,
                   color=GROUP_COLORS[gs['Group']],
                   edgecolor='white', linewidth=1)
            ax.errorbar(i, gs['Mean'], yerr=gs['SEM'], fmt='none',
                        color='black', capsize=4, linewidth=1.5)
        
        # Add significance brackets
        if p_values:
            all_means = [gs['Mean'] for gs in group_stats]
            all_sems = [gs['SEM'] for gs in group_stats]
            max_y = max([m + s for m, s in zip(all_means, all_sems)])
            min_y = min([m - s for m, s in zip(all_means, all_sems)])
            
            bracket_y = max_y * 1.1 if max_y > 0 else max_y + 0.05 * abs(min_y - max_y)
            
            for i, j, p in p_values:
                if p < 0.05:
                    bracket_y = add_significance_bracket(ax, i, j, bracket_y, p)
        
        ax.axhline(0, color='gray', linewidth=1, linestyle='--', alpha=0.5)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([GROUP_LABELS[gs['Group']] for gs in group_stats], fontsize=11)
        ax.set_ylabel(metric['ylabel'], fontsize=12, fontweight='bold')
        ax.set_title(metric['title'], fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        ymin, ymax = ax.get_ylim()
        y_range = ymax - ymin
        ax.set_ylim(ymin - 0.1 * y_range, ymax + 0.15 * y_range)
        
        plt.tight_layout()
        output_path = output_dir / metric['filename']
        plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
        plt.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
        plt.close()
        print(f"  ✓ {metric['filename']}")
        
        # =====================================================================
        # PLOT 2: Scatter plot (individual data points only, no bars)
        #         Outliers shown in grey, mean calculated without outliers
        # =====================================================================
        fig, ax = plt.subplots(figsize=(8, 6))
        
        for i, gs in enumerate(group_stats):
            vals = gs['values']
            
            # Detect outliers using IQR method
            q1, q3 = np.percentile(vals, [25, 75])
            iqr = q3 - q1
            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr
            
            outlier_mask = (vals < lower_bound) | (vals > upper_bound)
            normal_vals = vals[~outlier_mask]
            outlier_vals = vals[outlier_mask]
            
            # Jitter for normal points
            jitter_normal = np.random.default_rng(42 + i).uniform(-0.2, 0.2, len(normal_vals))
            ax.scatter(i + jitter_normal, normal_vals,
                       color=GROUP_COLORS[gs['Group']],
                       edgecolor='black', s=80, alpha=0.8, zorder=3, linewidths=0.5)
            
            # Jitter for outlier points (grey)
            if len(outlier_vals) > 0:
                jitter_outlier = np.random.default_rng(99 + i).uniform(-0.2, 0.2, len(outlier_vals))
                ax.scatter(i + jitter_outlier, outlier_vals,
                           color='lightgray',
                           edgecolor='darkgray', s=80, alpha=0.7, zorder=2, linewidths=0.5)
            
            # Add horizontal line for group mean (calculated WITHOUT outliers)
            if len(normal_vals) > 0:
                mean_no_outliers = np.mean(normal_vals)
                ax.hlines(mean_no_outliers, i - 0.3, i + 0.3, color='black', linewidth=2, zorder=4)
        
        # Add significance brackets based on individual data range
        if p_values:
            all_vals = np.concatenate([gs['values'] for gs in group_stats])
            max_y = np.max(all_vals)
            min_y = np.min(all_vals)
            
            bracket_y = max_y + 0.05 * (max_y - min_y)
            
            for i, j, p in p_values:
                if p < 0.05:
                    bracket_y = add_significance_bracket(ax, i, j, bracket_y, p)
        
        ax.axhline(0, color='gray', linewidth=1, linestyle='--', alpha=0.5)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([GROUP_LABELS[gs['Group']] for gs in group_stats], fontsize=11)
        ax.set_ylabel(metric['ylabel'], fontsize=12, fontweight='bold')
        ax.set_title(metric['title'] + '\n(Individual Subjects)', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        ymin, ymax = ax.get_ylim()
        y_range = ymax - ymin
        ax.set_ylim(ymin - 0.1 * y_range, ymax + 0.15 * y_range)
        
        plt.tight_layout()
        output_path = output_dir / f"{metric['filename']}_scatter"
        plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
        plt.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
        plt.close()
        print(f"  ✓ {metric['filename']}_scatter")


def plot_erp_waveforms(data_dir: Path, roi: str, output_dir: Path):
    """
    Plot grand-average ERP waveforms per group per task.
    One subplot per task, overlaid group waveforms with SEM shading.
    """
    wf_matrix, wf_info, times = load_erp_waveforms(data_dir, roi)
    if wf_matrix is None or wf_info is None:
        print("  ⚠ Waveform data not available, skipping waveform plots")
        return

    # Filter correct trials
    correct_mask = wf_info['Accuracy'] == 'correct'
    wf_matrix = wf_matrix[correct_mask.values]
    wf_info = wf_info[correct_mask].reset_index(drop=True)

    # Convert V → µV if needed
    if np.nanmedian(np.abs(wf_matrix)) < 0.001:
        wf_matrix = wf_matrix * 1e6

    erp_tasks = [t for t in ERP_TASK_ORDER if t in wf_info['Task'].unique()]
    n_tasks = len(erp_tasks)

    if n_tasks == 0:
        print("  ⚠ No ERP task data found")
        return

    fig, axes = plt.subplots(1, n_tasks, figsize=(4 * n_tasks, 5), sharey=True)
    if n_tasks == 1:
        axes = [axes]

    for ax, task in zip(axes, erp_tasks):
        task_mask = wf_info['Task'] == task

        for group in GROUP_ORDER:
            group_mask = task_mask & (wf_info['Group'] == group)
            trials = wf_matrix[group_mask.values]

            if len(trials) == 0:
                continue

            # Compute subject-level averages first, then grand average
            subj_ids = wf_info.loc[group_mask, 'SubjectID'].values
            unique_subjs = np.unique(subj_ids)
            subj_avgs = []
            for sid in unique_subjs:
                s_mask = group_mask & (wf_info['SubjectID'] == sid)
                s_trials = wf_matrix[s_mask.values]
                if len(s_trials) > 0:
                    subj_avgs.append(np.mean(s_trials, axis=0))

            if not subj_avgs:
                continue

            subj_avgs = np.array(subj_avgs)
            grand_mean = np.mean(subj_avgs, axis=0)
            grand_sem = np.std(subj_avgs, axis=0, ddof=1) / np.sqrt(len(subj_avgs))

            color = GROUP_COLORS[group]
            ax.plot(times, grand_mean, color=color, linewidth=1.5,
                    label=GROUP_LABELS[group])
            ax.fill_between(times, grand_mean - grand_sem, grand_mean + grand_sem,
                            color=color, alpha=0.2)

        # P300 window shading
        ax.axvspan(250, 500, color='gray', alpha=0.08)
        ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
        ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')

        ax.set_xlabel('Time (ms)', fontsize=10)
        ax.set_title(TASK_LABELS.get(task, task), fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].set_ylabel('Amplitude (µV)', fontsize=11, fontweight='bold')
    axes[-1].legend(loc='upper right', frameon=False, fontsize=9)

    fig.suptitle(f'Grand Average ERP Waveforms ({roi.upper()})',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()
    out = output_dir / f'erp_waveforms_{roi}'
    plt.savefig(out.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.savefig(out.with_suffix('.svg'), bbox_inches='tight')
    plt.close()
    print(f"  ✓ erp_waveforms_{roi}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    if len(sys.argv) >= 4:
        data_dir = Path(sys.argv[1])
        stats_dir = Path(sys.argv[2])
        output_dir = Path(sys.argv[3])
        roi = sys.argv[4] if len(sys.argv) >= 5 else 'pfc'
    else:
        print("Usage: python EEG_plotStatistics.py <data_dir> <stats_dir> <output_dir> [roi]")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n===== EEG PLOTTING =====")
    print(f"Data:   {data_dir}")
    print(f"Stats:  {stats_dir}")
    print(f"Output: {output_dir}")
    print(f"ROI:    {roi}")

    # --- Load subject-level data ---
    print("\nLoading data...")
    bp_subj = load_bandpower_subj_means(data_dir, roi)
    erp_subj = load_erp_subj_means(data_dir, roi)
    ratio_subj = compute_ratios_subj_means(data_dir, roi)

    if not bp_subj.empty:
        print(f"  Bandpower: {bp_subj['SubjectID'].nunique()} subjects")
    if not erp_subj.empty:
        print(f"  ERP: {erp_subj['SubjectID'].nunique()} subjects")
    if not ratio_subj.empty:
        print(f"  Ratios: {ratio_subj['SubjectID'].nunique()} subjects")

    # --- Load stats CSVs ---
    stats_rel = load_stats_csv(stats_dir, 'stats_bandpower_rel.csv')
    stats_log = load_stats_csv(stats_dir, 'stats_bandpower_log.csv')
    stats_abs = load_stats_csv(stats_dir, 'stats_bandpower_abs.csv')
    stats_ratios = load_stats_csv(stats_dir, 'stats_ratios.csv')
    stats_erp = load_stats_csv(stats_dir, 'stats_erp.csv')

    # =====================================================================
    #  BANDPOWER PLOTS — one figure per band, all 5 tasks (no markers)
    # =====================================================================
    if not bp_subj.empty:
        print("\n→ Plotting relative band power...")
        for band in BAND_ORDER:
            if band not in bp_subj['Band'].values:
                continue
            plot_grouped_bars(
                subj_data=bp_subj,
                stats_df=stats_rel,
                value_col='Value_rel',
                filter_col='Band', filter_val=band,
                tasks=TASK_ORDER,
                ylabel='Relative Power',
                title=f'Relative {BAND_LABELS[band]}',
                output_path=output_dir / f'bandpower_rel_{band}',
                stats_task_col='Task',
                stats_filter_col='Band', stats_filter_val=band,
                show_individual_markers=False,
            )

        print("\n→ Plotting log band power...")
        for band in BAND_ORDER:
            if band not in bp_subj['Band'].values:
                continue
            plot_grouped_bars(
                subj_data=bp_subj,
                stats_df=stats_log,
                value_col='Value_log',
                filter_col='Band', filter_val=band,
                tasks=TASK_ORDER,
                ylabel='Log Power (10·log₁₀ µV²)',
                title=f'Log {BAND_LABELS[band]}',
                output_path=output_dir / f'bandpower_log_{band}',
                stats_task_col='Task',
                stats_filter_col='Band', stats_filter_val=band,
                show_individual_markers=False,
            )

    # =====================================================================
    #  RATIO PLOTS — one figure per ratio, all 5 tasks + grand average
    # =====================================================================
    if not ratio_subj.empty:
        print("\n→ Plotting spectral ratios (per task)...")
        for ratio_name in RATIO_ORDER:
            if ratio_name not in ratio_subj['Ratio'].values:
                continue

            # stats CSV uses 'Ratio' column
            plot_grouped_bars(
                subj_data=ratio_subj,
                stats_df=stats_ratios,
                value_col='Value',
                filter_col='Ratio', filter_val=ratio_name,
                tasks=TASK_ORDER,
                ylabel='Ratio',
                title=RATIO_LABELS[ratio_name],
                output_path=output_dir / f'ratio_{ratio_name.replace("/", "_").lower()}',
                stats_task_col='Task',
                stats_filter_col='Ratio', stats_filter_val=ratio_name,
                show_individual_markers=False,
            )
        
        print("\n→ Plotting spectral ratios (grand average)...")
        for ratio_name in RATIO_ORDER:
            if ratio_name not in ratio_subj['Ratio'].values:
                continue
            
            plot_ratio_grand_average(
                subj_data=ratio_subj,
                stats_df=stats_ratios,
                ratio_name=ratio_name,
                output_path=output_dir / f'ratio_{ratio_name.replace("/", "_").lower()}_grand_avg',
                show_individual_markers=False,
            )

    # =====================================================================
    #  ERP METRIC PLOTS — bar plots for amplitude and latency (no markers)
    # =====================================================================
    if not erp_subj.empty:
        print("\n→ Plotting ERP metrics...")
        for metric in ERP_METRIC_ORDER:
            if metric not in erp_subj.columns:
                continue

            # Filter stats to this metric
            plot_grouped_bars(
                subj_data=erp_subj,
                stats_df=stats_erp,
                value_col=metric,
                filter_col=None, filter_val=None,
                tasks=ERP_TASK_ORDER,
                ylabel=ERP_METRIC_LABELS[metric],
                title=ERP_METRIC_LABELS[metric],
                output_path=output_dir / f'erp_{metric.lower()}',
                stats_task_col='Task',
                stats_filter_col='ERPMetric', stats_filter_val=metric,
                show_individual_markers=False,
            )

    # =====================================================================
    #  DELTA DAR PLOT — Change in Delta/Alpha Ratio (2-back - Eyes Open)
    # =====================================================================
    print("\n→ Plotting Delta DAR...")
    try:
        plot_delta_dar(stats_dir, output_dir)
    except Exception as e:
        print(f"  ⚠ Delta DAR plotting failed: {e}")

    # =====================================================================
    #  ERP WAVEFORM PLOTS — grand-average with SEM shading
    # =====================================================================
    print("\n→ Plotting ERP waveforms...")
    try:
        plot_erp_waveforms(data_dir, roi, output_dir)
    except Exception as e:
        print(f"  ⚠ Waveform plotting failed: {e}")

    print(f"\n===== PLOTTING COMPLETE =====")
    print(f"Figures saved to {output_dir}")


if __name__ == '__main__':
    main()
