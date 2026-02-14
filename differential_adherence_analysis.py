#!/usr/bin/env python3
"""
Differential Adherence Framework: Orthogonal Decomposition vs Sum-Score Flattening
====================================================================================
By: Anderson M. Rodriguez
https://github.com/amr28693/orthogonal_differential_adherence

Demonstrates that the Differential Adherence Index (Δ = intentional − unintentional)
captures statistically significant clinical heterogeneity that sum-scoring
(as used by MMAS and similar instruments) provably destroys.

Data: NeuroGerAd study (N=907 neurological patients)
      Prell et al. (2022) Sci Data 9, 734. doi:10.1038/s41597-022-01847-9
      Available: https://osf.io/kuaph/

Instrument: SAMS (Stendal Adherence to Medication Score) — 18 items, 0–4 Likert
      Sub-factors per validated 3-factor CFA (Prell et al. 2022):
        - Forgetting (unintentional): items 6, 14, 15, 16, 18
        - Intentional modification:  items 4, 7, 8, 9, 10, 11, 12, 13, 17
        - Missing knowledge:         items 1, 2, 3, 5

Usage:
    Place NeuroGerAd_Data_OSF.xlsx in the same folder as this script, then:
    $ python differential_adherence_analysis.py
    
    Following this initial script, run:
    $sensitivity_analysis.py

Outputs:
    fig1_3d_vs_sumscore.png         — 3D decomposition vs honest 1D sum-score strip
    fig2_group_heterogeneity.png    — Δ distributions + group comparisons
    fig3_information_loss.png       — Orthogonality + F-statistic comparison
    results_summary.txt             — All statistical results
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


# ============================================================
# CONFIGURATION
# ============================================================
DATA_FILE = 'NeuroGerAd_Data_OSF.xlsx'
FIGURE_DPI = 300

# SAMS sub-factor item mapping (validated 3-factor CFA, Prell et al. 2022)
# Scale: 0 = "never" to 4 = "most of the time" (higher = more nonadherent)
UNINTENTIONAL_ITEMS = ['sams_6', 'sams_14', 'sams_15', 'sams_16', 'sams_18']
INTENTIONAL_ITEMS   = ['sams_4', 'sams_7', 'sams_8', 'sams_9', 'sams_10',
                       'sams_11', 'sams_12', 'sams_13', 'sams_17']
KNOWLEDGE_ITEMS     = ['sams_1', 'sams_2', 'sams_3', 'sams_5']
DIAGNOSIS_COL       = 'diagnosis_collapsed'


# ============================================================
# 1. DATA LOADING
# ============================================================
def load_data():
    """Load NeuroGerAd Excel file and compute sub-factor means."""
    filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), DATA_FILE)
    if not os.path.exists(filepath):
        filepath = DATA_FILE
    if not os.path.exists(filepath):
        print(f"ERROR: Cannot find '{DATA_FILE}'")
        print(f"Place NeuroGerAd_Data_OSF.xlsx in the same folder as this script.")
        print(f"Download from: https://osf.io/kuaph/")
        sys.exit(1)

    print(f"Loading: {filepath}")
    raw = pd.read_excel(filepath)
    print(f"  Raw: {len(raw)} rows, {len(raw.columns)} columns")

    # Sub-factor means (means, not sums — comparable across unequal item counts)
    # Person-mean scoring: average available items per subscale (skipna=True).
    # Standard approach when items within a subscale are parallel indicators.
    # Participants are retained if they have >= 1 item answered per subscale.
    raw['unintentional'] = raw[UNINTENTIONAL_ITEMS].mean(axis=1, skipna=True)
    raw['intentional']   = raw[INTENTIONAL_ITEMS].mean(axis=1, skipna=True)
    raw['knowledge']     = raw[KNOWLEDGE_ITEMS].mean(axis=1, skipna=True)
    raw['group']         = raw[DIAGNOSIS_COL]

    df = raw.dropna(subset=['unintentional', 'intentional', 'group']).copy()
    print(f"  Valid: {len(df)} (dropped {len(raw) - len(df)} with all items missing in a subscale or missing diagnosis)")

    # Core computed variables
    df['delta']      = df['intentional'] - df['unintentional']
    df['mmas_proxy'] = df['intentional'] + df['unintentional']
    df['direction']  = np.where(df['delta'] >= 0,
                                'Intentional-dominant', 'Unintentional-dominant')
    return df


# ============================================================
# 2. DESCRIPTIVE STATISTICS
# ============================================================
def print_descriptives(df, log):
    """Print and log descriptive statistics."""
    lines = []
    lines.append("=" * 72)
    lines.append("DESCRIPTIVE STATISTICS")
    lines.append("=" * 72)
    lines.append(f"N = {len(df)}")
    lines.append(f"Diagnosis groups: {dict(df['group'].value_counts())}")
    lines.append("")
    lines.append("Sub-factor means (0-4 scale, higher = more nonadherent):")
    lines.append(f"  Unintentional (forgetting, 5 items):  "
                 f"M = {df['unintentional'].mean():.3f}, SD = {df['unintentional'].std():.3f}")
    lines.append(f"  Intentional (modification, 9 items):  "
                 f"M = {df['intentional'].mean():.3f}, SD = {df['intentional'].std():.3f}")
    lines.append(f"  Missing knowledge (4 items):          "
                 f"M = {df['knowledge'].mean():.3f}, SD = {df['knowledge'].std():.3f}")
    lines.append("")
    lines.append("Differential Adherence Index (delta = intentional - unintentional):")
    lines.append(f"  Mean = {df['delta'].mean():.3f}, SD = {df['delta'].std():.3f}")
    lines.append(f"  Range = [{df['delta'].min():.3f}, {df['delta'].max():.3f}]")
    lines.append(f"  Intentional-dominant (delta >= 0): "
                 f"{(df['delta'] >= 0).sum()} ({(df['delta'] >= 0).mean()*100:.1f}%)")
    lines.append(f"  Unintentional-dominant (delta < 0): "
                 f"{(df['delta'] < 0).sum()} ({(df['delta'] < 0).mean()*100:.1f}%)")
    lines.append("")
    lines.append("Sum-score proxy (MMAS-like = intentional + unintentional):")
    lines.append(f"  Mean = {df['mmas_proxy'].mean():.3f}, SD = {df['mmas_proxy'].std():.3f}")
    lines.append("")
    lines.append("-" * 72)
    lines.append("GROUP-LEVEL STATISTICS")
    lines.append("-" * 72)
    hdr = f"{'Group':<28} {'N':>5} {'Δ Mean':>8} {'Δ SD':>7} {'%Int':>6} {'Sum Mean':>9} {'Sum SD':>7}"
    lines.append(hdr)
    lines.append("-" * 72)
    for g in sorted(df['group'].unique()):
        sub = df[df['group'] == g]
        lines.append(f"{g:<28} {len(sub):>5} {sub['delta'].mean():>8.3f} "
                     f"{sub['delta'].std():>7.3f} "
                     f"{(sub['delta'] >= 0).mean()*100:>5.1f} "
                     f"{sub['mmas_proxy'].mean():>9.3f} "
                     f"{sub['mmas_proxy'].std():>7.3f}")
    lines.append("")

    text = "\n".join(lines)
    print(text)
    log.append(text)


# ============================================================
# 3. STATISTICAL TESTS
# ============================================================
def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


def run_statistics(df, log):
    """Run and log all statistical tests. Returns dict of key results."""
    lines = []
    lines.append("=" * 72)
    lines.append("STATISTICAL TESTS")
    lines.append("=" * 72)

    # Orthogonality
    r, p = stats.pearsonr(df['mmas_proxy'], df['delta'])
    lines.append("Orthogonality (Pearson correlation, sum vs delta):")
    lines.append(f"  r = {r:.3f}, p = {p:.4f}")
    lines.append(f"  Shared variance (r\u00b2): {r**2*100:.1f}%")
    lines.append(f"  -> {'Near-perfect orthogonality confirmed' if abs(r) < 0.1 else 'Some shared variance'}")
    lines.append("")

    # ANOVA
    groups = sorted(df['group'].unique())
    k, n = len(groups), len(df)
    grp_delta = [df[df['group'] == g]['delta'].values for g in groups]
    grp_mmas  = [df[df['group'] == g]['mmas_proxy'].values for g in groups]

    f_d, p_d = stats.f_oneway(*grp_delta)
    f_m, p_m = stats.f_oneway(*grp_mmas)
    eta_d = (f_d * (k-1)) / (f_d * (k-1) + (n-k))
    eta_m = (f_m * (k-1)) / (f_m * (k-1) + (n-k))

    lines.append("One-way ANOVA by diagnosis group:")
    lines.append(f"  delta (differential):  F({k-1},{n-k}) = {f_d:.2f}, "
                 f"p = {p_d:.4f} {sig_label(p_d)}, eta2 = {eta_d:.4f}")
    lines.append(f"  Sum (MMAS-like):       F({k-1},{n-k}) = {f_m:.2f},  "
                 f"p = {p_m:.4f} {sig_label(p_m)}, eta2 = {eta_m:.4f}")
    lines.append("")
    if p_d < 0.05 and p_m >= 0.05:
        lines.append("  * Groups differ significantly on delta but NOT on sum-score.")
        lines.append("    Sum-scoring destroys individual-level directionality (many-to-one projection).")
    elif eta_d > eta_m:
        lines.append(f"  -> delta captures {eta_d/max(eta_m,1e-6):.1f}x more "
                     f"between-group variance than sum-score")
    lines.append("")

    # Pairwise t-tests
    lines.append("Pairwise comparisons on delta (independent samples t-test):")
    lines.append(f"  {'Comparison':<50} {'t':>7} {'p':>8} {'d':>7} {'sig':>4}")
    lines.append("  " + "-" * 76)
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            g1, g2 = groups[i], groups[j]
            d1 = df[df['group'] == g1]['delta'].values
            d2 = df[df['group'] == g2]['delta'].values
            t_stat, p_val = stats.ttest_ind(d1, d2)
            pooled_sd = np.sqrt((d1.std()**2 + d2.std()**2) / 2)
            d_val = abs((d1.mean() - d2.mean()) / pooled_sd) if pooled_sd > 0 else 0
            label = f"{g1} vs {g2}"
            lines.append(f"  {label:<50} {t_stat:>7.2f} "
                         f"{p_val:>8.4f} {d_val:>7.2f} {sig_label(p_val):>4}")
    lines.append("")

    text = "\n".join(lines)
    print(text)
    log.append(text)

    return dict(r=r, r_p=p, f_delta=f_d, p_delta=p_d, eta_delta=eta_d,
                f_mmas=f_m, p_mmas=p_m, eta_mmas=eta_m, k=k, n=n)


# ============================================================
# 4. FIGURE 1: 3D Decomposition vs Honest 1D Sum-Score Strip
# ============================================================
def plot_figure1(df):
    fig = plt.figure(figsize=(16, 7))

    # --- Panel A: 3D orthogonal decomposition ---
    ax1 = fig.add_subplot(121, projection='3d')
    colors = np.where(df['delta'] >= 0, '#e74c3c', '#3498db')
    ax1.scatter(df['unintentional'], df['intentional'], df['delta'],
                c=colors, s=8, alpha=0.45, edgecolors='none')

    xx, yy = np.meshgrid(
        np.linspace(df['unintentional'].min(), df['unintentional'].max(), 5),
        np.linspace(df['intentional'].min(), df['intentional'].max(), 5))
    ax1.plot_surface(xx, yy, np.zeros_like(xx), alpha=0.06, color='gray')

    ax1.set_xlabel('Unintentional\n(Forgetting)', fontsize=10)
    ax1.set_ylabel('Intentional\n(Modification)', fontsize=10)
    ax1.set_zlabel(r'$\Delta$', fontsize=11)
    ax1.set_title(f'A) Orthogonal Decomposition (N={len(df)})',
                  fontsize=12, fontweight='bold')
    ax1.view_init(elev=22, azim=45)

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c',
               markersize=7, label=r'Intentional-dominant ($\Delta \geq 0$)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498db',
               markersize=7, label=r'Unintentional-dominant ($\Delta < 0$)')]
    ax1.legend(handles=legend_elements, loc='upper left', fontsize=7.5)

    # --- Panel B: What sum-scoring actually gives you ---
    ax2 = fig.add_subplot(122)
    jitter = np.random.default_rng(42).normal(0, 0.015, len(df))
    ax2.scatter(df['mmas_proxy'], jitter,
                c='#7f8c8d', s=6, alpha=0.25, edgecolors='none')

    ax2.set_xlabel('Sum-Score (Intentional + Unintentional)', fontsize=11)
    ax2.set_yticks([])
    ax2.set_ylim(-0.08, 0.08)
    ax2.set_title('B) Sum-Score Reduction (what you actually get)',
                  fontsize=12, fontweight='bold')
    ax2.text(0.97, 0.03, 'No color or direction.\nSimple single number.',
             transform=ax2.transAxes, ha='right', va='bottom',
             fontsize=9, color='#666', fontstyle='italic')

    plt.tight_layout()
    plt.savefig('fig1_3d_vs_sumscore.png', dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print("Saved: fig1_3d_vs_sumscore.png")


# ============================================================
# 5. FIGURE 2: Group-Level Heterogeneity
# ============================================================
def plot_figure2(df, sr):
    groups = sorted(df['group'].unique())
    gc = dict(zip(groups, plt.cm.Set2(np.linspace(0, 1, len(groups)))))
    gm = df.groupby('group')[['unintentional', 'intentional',
                               'delta', 'mmas_proxy']].mean()

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # --- A: delta distributions ---
    ax = axes[0]
    for g in groups:
        sub = df[df['group'] == g]
        ax.hist(sub['delta'], bins=30, alpha=0.4,
                label=f"{g} ({len(sub)})", color=gc[g], density=True)
    ax.axvline(0, color='black', linestyle='--', alpha=0.4, lw=0.8)
    ax.set_xlabel(r'$\Delta$ (intentional $-$ unintentional)', fontsize=10)
    ax.set_ylabel('Density', fontsize=10)
    ax.set_title(r'A) $\Delta$ Distribution by Diagnosis',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=6.5, loc='upper right', framealpha=0.7)

    # --- B: Group means — sum-score view ---
    ax = axes[1]
    for score in np.arange(0.3, 1.5, 0.15):
        xl = np.linspace(0, score, 100)
        yl = score - xl
        m = (xl >= 0.1) & (yl >= 0.05)
        if m.any():
            ax.plot(xl[m], yl[m], '-', color='#ccc', alpha=0.3, lw=0.7)
    for g in groups:
        row = gm.loc[g]
        ax.scatter(row['unintentional'], row['intentional'],
                   c=[gc[g]], s=160, zorder=5, edgecolors='black', linewidth=1)
        ax.annotate(g, (row['unintentional'] + 0.008,
                        row['intentional'] + 0.008), fontsize=7)
    ax.set_xlabel('Unintentional (mean)', fontsize=10)
    ax.set_ylabel('Intentional (mean)', fontsize=10)
    ax.set_title(f"B) Sum-Score View: F={sr['f_mmas']:.2f}, p={sr['p_mmas']:.3f}",
                 fontsize=11, fontweight='bold')

    # --- C: Group delta bar chart ---
    ax = axes[2]
    sorted_gm = gm.sort_values('delta')
    for i, (g, row) in enumerate(sorted_gm.iterrows()):
        ax.barh(i, row['delta'], color=gc[g], edgecolor='black',
                linewidth=0.5, height=0.6)
        n_g = len(df[df['group'] == g])
        offset = 0.008 if row['delta'] < 0 else -0.008
        ha = 'right' if row['delta'] < 0 else 'left'
        ax.text(row['delta'] + offset, i, f'n={n_g}',
                fontsize=7.5, va='center', ha=ha, color='#444')
    ax.axvline(0, color='black', linestyle='-', alpha=0.4, lw=0.8)
    ax.set_yticks(range(len(sorted_gm)))
    ax.set_yticklabels(sorted_gm.index, fontsize=9)
    ax.set_xlabel(r'Mean $\Delta$', fontsize=10)
    ax.set_title(r'C) $\Delta$ by Diagnosis: F={:.2f}, p={:.3f}*'.format(
                 sr['f_delta'], sr['p_delta']),
                 fontsize=11, fontweight='bold')

    plt.tight_layout()
    plt.savefig('fig2_group_heterogeneity.png', dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print("Saved: fig2_group_heterogeneity.png")


# ============================================================
# 6. FIGURE 3: Information Loss Quantification
# ============================================================
def plot_figure3(df, sr):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # --- A: Orthogonality scatter ---
    ax = axes[0]
    colors = np.where(df['delta'] >= 0, '#e74c3c', '#3498db')
    ax.scatter(df['mmas_proxy'], df['delta'], c=colors, s=6, alpha=0.3,
               edgecolors='none')
    ax.axhline(0, color='black', linestyle='--', alpha=0.25, lw=0.7)
    ax.text(0.03, 0.97,
            f"r = {sr['r']:.3f}\n{(1-sr['r']**2)*100:.1f}% independent",
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#ccc', alpha=0.85))
    ax.set_xlabel('Sum-Score', fontsize=11)
    ax.set_ylabel(r'$\Delta$', fontsize=11)
    ax.set_title(r'A) Orthogonality: Sum vs $\Delta$',
                 fontsize=12, fontweight='bold')

    # --- B: F-statistic comparison ---
    ax = axes[1]
    labels = ['Sum-Score', r'$\Delta$']
    f_vals = [sr['f_mmas'], sr['f_delta']]
    p_vals = [sr['p_mmas'], sr['p_delta']]
    bar_colors = ['#bdc3c7', '#2ecc71']

    bars = ax.bar(labels, f_vals, color=bar_colors, edgecolor='black',
                  linewidth=0.8, width=0.4)
    for bar, fv, pv in zip(bars, f_vals, p_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.06,
                f'F={fv:.2f}\np={pv:.3f} {sig_label(pv)}',
                ha='center', fontsize=10, fontweight='bold')

    f_crit = stats.f.ppf(0.95, sr['k'] - 1, sr['n'] - sr['k'])
    ax.axhline(y=f_crit, color='red', linestyle=':', alpha=0.4, lw=1,
               label=f'p=.05 (F={f_crit:.2f})')
    ax.set_ylabel('F-statistic (ANOVA by diagnosis)', fontsize=11)
    ax.set_title('B) Between-Group Signal', fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_ylim(0, max(f_vals) * 1.5)

    plt.tight_layout()
    plt.savefig('fig3_information_loss.png', dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print("Saved: fig3_information_loss.png")


# ============================================================
# MAIN
# ============================================================
def main():
    log = []
    log.append("Differential Adherence Analysis -- NeuroGerAd Real Data")
    log.append(f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}")
    log.append("")

    df = load_data()
    log.append(f"Data loaded: N={len(df)}, {df['group'].nunique()} diagnosis groups\n")

    print_descriptives(df, log)
    sr = run_statistics(df, log)

    print("\nGenerating figures...")
    plot_figure1(df)
    plot_figure2(df, sr)
    plot_figure3(df, sr)

    with open('results_summary.txt', 'w') as f:
        f.write("\n".join(log))
    print("Saved: results_summary.txt")
    print("\nDone.")


if __name__ == '__main__':
    main()
