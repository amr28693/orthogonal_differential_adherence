#!/usr/bin/env python3
"""
Sensitivity Analysis: Differential Adherence Under Varying Missing-Data Thresholds
===================================================================================
By: Anderson M. Rodriguez

Companion to differential_adherence_analysis.py. Reproduces Table 2 from:

    Rodriguez AM. Differential Adherence: Orthogonal Decomposition of Medication
    Nonadherence Reveals Clinically Significant Heterogeneity Destroyed by
    Sum-Score Instruments. Submitted, 2026.

Reruns omnibus ANOVA and the key pairwise contrast (movement disorder vs
neuromuscular) under four progressively stricter inclusion criteria.

Usage:
    python sensitivity_analysis.py

Requires NeuroGerAd_Data_OSF.xlsx in the same directory.
Download from: https://osf.io/kuaph/
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

# ============================================================
# CONFIGURATION (must match differential_adherence_analysis.py)
# ============================================================
DATA_FILE = 'NeuroGerAd_Data_OSF.xlsx'

UNINTENTIONAL_ITEMS = ['sams_6', 'sams_14', 'sams_15', 'sams_16', 'sams_18']
INTENTIONAL_ITEMS   = ['sams_4', 'sams_7', 'sams_8', 'sams_9', 'sams_10',
                       'sams_11', 'sams_12', 'sams_13', 'sams_17']
DIAGNOSIS_COL       = 'diagnosis_collapsed'

# Minimum items answered per subscale for each threshold
# Unintentional: 5 items total; Intentional: 9 items total
THRESHOLDS = [
    ("Person-mean (primary)",  0,   0  ),   # skipna=True, no minimum
    ("≥50% items",             3,   5  ),   # ceil(5*0.5)=3, ceil(9*0.5)=5
    ("≥80% items",             4,   8  ),   # ceil(5*0.8)=4, ceil(9*0.8)=8 (round up)
    ("Complete-case",          5,   9  ),   # all items answered
]

KEY_CONTRAST = ("movement disorder", "neuromuscular")


# ============================================================
# DATA LOADING
# ============================================================
def load_raw():
    """Load raw Excel file, return DataFrame with SAMS items and diagnosis."""
    filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), DATA_FILE)
    if not os.path.exists(filepath):
        filepath = DATA_FILE
    if not os.path.exists(filepath):
        print(f"ERROR: Cannot find '{DATA_FILE}'")
        print(f"Place NeuroGerAd_Data_OSF.xlsx in the same folder as this script.")
        print(f"Download from: https://osf.io/kuaph/")
        sys.exit(1)
    return pd.read_excel(filepath)


def apply_threshold(raw, min_unint, min_intent):
    """Score subscales and apply item-count threshold.

    For threshold=0 (person-mean), uses skipna=True with no minimum.
    For threshold>0, requires at least that many non-missing items per subscale.
    """
    df = raw.copy()

    # Count non-missing items per subscale
    df['n_unint']  = df[UNINTENTIONAL_ITEMS].notna().sum(axis=1)
    df['n_intent'] = df[INTENTIONAL_ITEMS].notna().sum(axis=1)

    # Person-mean scoring (always skipna=True; threshold enforced by masking)
    df['unintentional'] = df[UNINTENTIONAL_ITEMS].mean(axis=1, skipna=True)
    df['intentional']   = df[INTENTIONAL_ITEMS].mean(axis=1, skipna=True)

    # Apply threshold: set to NaN if below minimum item count
    if min_unint > 0:
        df.loc[df['n_unint'] < min_unint, 'unintentional'] = np.nan
    if min_intent > 0:
        df.loc[df['n_intent'] < min_intent, 'intentional'] = np.nan

    df['group'] = df[DIAGNOSIS_COL]
    df = df.dropna(subset=['unintentional', 'intentional', 'group']).copy()

    df['delta']      = df['intentional'] - df['unintentional']
    df['mmas_proxy'] = df['intentional'] + df['unintentional']

    return df


# ============================================================
# ANALYSIS
# ============================================================
def run_sensitivity(df):
    """Run omnibus ANOVA + key pairwise contrast. Return dict of results."""
    groups = sorted(df['group'].unique())
    k = len(groups)
    n = len(df)

    # Omnibus ANOVA
    grp_delta = [df[df['group'] == g]['delta'].values for g in groups]
    grp_mmas  = [df[df['group'] == g]['mmas_proxy'].values for g in groups]

    f_d, p_d = stats.f_oneway(*grp_delta)
    f_m, p_m = stats.f_oneway(*grp_mmas)

    # Key pairwise contrast: movement disorder vs neuromuscular
    g1, g2 = KEY_CONTRAST
    d1 = df[df['group'] == g1]['delta'].values
    d2 = df[df['group'] == g2]['delta'].values
    t_stat, p_pw = stats.ttest_ind(d1, d2)
    pooled_sd = np.sqrt((d1.std()**2 + d2.std()**2) / 2)
    d_cohen = abs((d1.mean() - d2.mean()) / pooled_sd) if pooled_sd > 0 else 0

    return dict(
        n=n, k=k,
        f_delta=f_d, p_delta=p_d,
        p_mmas=p_m,
        t_pw=t_stat, p_pw=p_pw, d_pw=d_cohen,
        df_pw=len(d1) + len(d2) - 2,
    )


# ============================================================
# MAIN
# ============================================================
def main():
    raw = load_raw()
    print(f"Loaded: {len(raw)} rows\n")

    results = []
    for label, min_u, min_i in THRESHOLDS:
        df = apply_threshold(raw, min_u, min_i)
        r = run_sensitivity(df)
        r['label'] = label
        results.append(r)
        print(f"  {label}: N={r['n']}, Δ ANOVA F={r['f_delta']:.2f} "
              f"p={r['p_delta']:.4f}, pairwise p={r['p_pw']:.4f} d={r['d_pw']:.2f}")

    # Format Table 2
    lines = []
    lines.append("")
    lines.append("=" * 90)
    lines.append("Table 2. Sensitivity Analysis Across Missing-Data Handling Approaches")
    lines.append("=" * 90)
    lines.append(f"{'Threshold':<28} {'N':>5}  {'Δ F':>7}  {'Δ P':>8}  "
                 f"{'Sum P':>8}  {'t(pw)':>7}  {'P(pw)':>8}  {'d':>5}")
    lines.append("-" * 90)
    for r in results:
        sig = "***" if r['p_delta'] < .001 else ("**" if r['p_delta'] < .01 else
              ("*" if r['p_delta'] < .05 else ""))
        lines.append(f"{r['label']:<28} {r['n']:>5}  {r['f_delta']:>7.2f}  "
                     f"{r['p_delta']:>8.4f}  {r['p_mmas']:>8.4f}  "
                     f"{r['t_pw']:>7.2f}  {r['p_pw']:>8.4f}  {r['d_pw']:>5.2f}  {sig}")
    lines.append("-" * 90)
    lines.append("Δ F, Δ P = omnibus ANOVA for Differential Adherence Index across 5 groups.")
    lines.append("Sum P = omnibus ANOVA P for sum-score. t(pw), P(pw), d = movement disorder")
    lines.append("vs neuromuscular pairwise contrast (Cohen |d|). All P values 2-tailed.")
    lines.append("* P<.05  ** P<.01  *** P<.001")
    lines.append("=" * 90)

    table_text = "\n".join(lines)
    print(table_text)

    # Save
    outfile = "sensitivity_results.txt"
    with open(outfile, 'w') as f:
        f.write(table_text)
    print(f"\nSaved: {outfile}")


if __name__ == '__main__':
    main()
