#!/usr/bin/env python
"""
ROC analysis for BiotinScore threshold optimization.

Reads pre-computed BiotinScore values from biotin_scoring.py output
(positive = nick-site-associated reads, negative = control reads),
generates an ROC curve at integer thresholds (step size = 1), and selects
the optimal threshold by maximizing Youden's J statistic.

Usage:
  python roc_curve.py --positive <pos1.xlsx> [pos2.xlsx ...] \
                      --negative <neg1.xlsx> [neg2.xlsx ...] [-o OUTPUT_DIR]
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Liberation Sans'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Liberation Sans'


def load_scores(xlsx_path):
    """
    Load BiotinScore values from biotin_scoring.py output.

    Reads the 'biotin_score' column from the 'per_read' sheet.

    Args:
        xlsx_path: Path to .xlsx output of biotin_scoring.py.

    Returns:
        1D numpy array of BiotinScore values.
    """
    df = pd.read_excel(xlsx_path, sheet_name="per_read")
    return df['biotin_score'].to_numpy(dtype=float)


def compute_roc_integer(pos_scores, neg_scores):
    """
    Compute ROC curve at integer thresholds (step size = 1).

    For each integer threshold t, sensitivity (TPR) is the fraction of
    positive reads with BiotinScore >= t, and false positive rate (FPR) is
    the fraction of negative reads with BiotinScore >= t.

    Args:
        pos_scores: 1D array of BiotinScore values (positive/nick-site reads).
        neg_scores: 1D array of BiotinScore values (negative/control reads).

    Returns:
        tuple: (fpr, tpr, thresholds, roc_auc) — all as numpy arrays / float.
    """
    n_pos = len(pos_scores)
    n_neg = len(neg_scores)

    t_min = int(np.floor(min(pos_scores.min(), neg_scores.min())))
    t_max = int(np.ceil(max(pos_scores.max(), neg_scores.max()))) + 1

    thresholds = np.arange(t_max, t_min - 1, -1)  # descending order
    tpr = np.array([(pos_scores >= t).sum() / n_pos for t in thresholds])
    fpr = np.array([(neg_scores >= t).sum() / n_neg for t in thresholds])

    # AUC via trapezoidal rule
    sorted_idx = np.argsort(fpr)
    fpr_sorted = fpr[sorted_idx]
    tpr_sorted = tpr[sorted_idx]
    roc_auc = float(np.trapz(tpr_sorted, fpr_sorted))

    return fpr, tpr, thresholds, roc_auc


def find_optimal_threshold(fpr, tpr, thresholds):
    """
    Find the optimal threshold by maximizing Youden's J statistic.

    Args:
        fpr: False positive rate array.
        tpr: True positive rate array.
        thresholds: Integer threshold array.

    Returns:
        tuple: (best_threshold, best_tpr, best_fpr)
    """
    j = tpr - fpr
    best_idx = np.argmax(j)
    return int(thresholds[best_idx]), tpr[best_idx], fpr[best_idx]


def plot_roc(fpr, tpr, roc_auc, best_int_fpr, best_int_tpr, output_dir):
    """
    Plot ROC curve with inset zoom and save as PNG and SVG.

    Args:
        fpr: False positive rate array.
        tpr: True positive rate array.
        roc_auc: Area under the ROC curve.
        best_int_fpr: FPR at the optimal threshold.
        best_int_tpr: TPR at the optimal threshold.
        output_dir: Directory for output files.
    """
    fig, ax = plt.subplots(figsize=(4, 4))

    ax.plot(fpr, tpr, color='#2c3e6b', lw=2.5)
    ax.plot([0, 1], [0, 1], color='gray', lw=1.2, linestyle='--')
    ax.plot(best_int_fpr, best_int_tpr, 'o', color='red', markersize=7, zorder=5)

    ax.set_xlim([-0.01, 1])
    ax.set_ylim([-0.01, 1.01])
    ax.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=14)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=14)
    ax.tick_params(labelsize=12)
    ax.set_aspect('equal')

    # Inset zoom
    axins = ax.inset_axes([0.50, 0.28, 0.42, 0.42])
    axins.plot(fpr, tpr, color='#2c3e6b', lw=2.5)
    axins.plot(best_int_fpr, best_int_tpr, 'o', color='red', markersize=6, zorder=5)

    x_max = min(max(best_int_fpr * 8, 0.06), 0.10)
    y_min = max(min(best_int_tpr - 0.05, 0.93), 0.90)
    axins.set_xlim([-0.001, x_max])
    axins.set_ylim([y_min, 1.005])
    axins.set_xlabel('FPR', fontsize=10)
    axins.set_ylabel('TPR', fontsize=10)
    axins.tick_params(labelsize=9)

    ax.text(0.72, 0.10, f"AUC = {roc_auc:.3f}",
            transform=ax.transAxes, fontsize=18, fontweight='bold',
            ha='center', va='center')

    for ext in ['png', 'svg']:
        path = os.path.join(output_dir, f"roc_curve.{ext}")
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"Saved: {path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="ROC analysis for BiotinScore threshold optimization.")
    parser.add_argument("--positive", required=True, nargs="+",
                        help="biotin_scoring.py output .xlsx "
                             "(nick-site-associated reads, one or more files)")
    parser.add_argument("--negative", required=True, nargs="+",
                        help="biotin_scoring.py output .xlsx "
                             "(control reads, one or more files)")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="Output directory (default: same as positive input)")
    args = parser.parse_args()

    output_dir = args.output_dir or os.path.dirname(os.path.abspath(args.positive[0]))

    # Load and pool pre-computed BiotinScore values
    pos_parts = []
    for path in args.positive:
        scores = load_scores(path)
        print(f"Loaded positive: {path} ({len(scores)} reads)")
        pos_parts.append(scores)
    pos_scores = np.concatenate(pos_parts)
    print(f"Total positive: {len(pos_scores)} reads")

    neg_parts = []
    for path in args.negative:
        scores = load_scores(path)
        print(f"Loaded negative: {path} ({len(scores)} reads)")
        neg_parts.append(scores)
    neg_scores = np.concatenate(neg_parts)
    print(f"Total negative: {len(neg_scores)} reads")

    # ROC analysis at integer thresholds
    fpr, tpr, thresholds, roc_auc = compute_roc_integer(pos_scores, neg_scores)
    print(f"AUC = {roc_auc:.6f}")

    # Optimal threshold (Youden's J)
    best_thresh, best_tpr, best_fpr = find_optimal_threshold(fpr, tpr, thresholds)
    print(f"Optimal threshold = {best_thresh}")
    print(f"  Sensitivity = {best_tpr * 100:.1f}%, "
          f"Specificity = {(1 - best_fpr) * 100:.1f}%")

    # Save ROC data
    xlsx_path = os.path.join(output_dir, "roc_data.xlsx")
    j_statistic = tpr - fpr
    df_roc = pd.DataFrame({
        "threshold": thresholds,
        "fpr (1-specificity)": fpr,
        "tpr (sensitivity)": tpr,
        "Youden's J": j_statistic,
    })
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as w:
        df_roc.to_excel(w, index=False, sheet_name="roc")
        pd.DataFrame({
            "auc": [roc_auc],
            "optimal_threshold": [best_thresh],
            "sensitivity": [best_tpr],
            "specificity": [1 - best_fpr],
        }).to_excel(w, index=False, sheet_name="summary")
    print(f"Saved: {xlsx_path}")

    # Plot
    plot_roc(fpr, tpr, roc_auc, best_fpr, best_tpr, output_dir)
    print("Done!")


if __name__ == "__main__":
    main()
