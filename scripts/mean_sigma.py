#!/usr/bin/env python
"""
Compute per-position mean and mean-2sigma threshold from per-base signal CSV.

Reads the output of extract_perbase_signal.py, applies a QC mask that excludes
positions near alignment boundaries, and computes summary statistics using
Welford's online algorithm for numerical stability.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple


# ============================================================================
# Welford online accumulator (column-wise, NaN-aware)
# ============================================================================

class WelfordVec:
    """Column-wise Welford online accumulator for streaming mean/variance."""

    def __init__(self, n_cols: int):
        self.count = np.zeros(n_cols, dtype=np.int64)
        self.mean = np.zeros(n_cols, dtype=np.float64)
        self.M2 = np.zeros(n_cols, dtype=np.float64)

    def update(self, x: np.ndarray):
        """Update with batch x of shape (r, C), ignoring NaN values."""
        valid = ~np.isnan(x)
        k = valid.sum(axis=0).astype(np.int64)
        if not np.any(k):
            return

        x_filled = np.where(valid, x, 0.0)
        sum_new = x_filled.sum(axis=0)
        count_prev = self.count.copy()
        mean_prev = self.mean.copy()
        M2_prev = self.M2.copy()

        self.count += k

        mean_batch = np.divide(sum_new, k, out=np.zeros_like(sum_new), where=(k > 0))
        sumsq_new = np.where(valid, x * x, 0.0).sum(axis=0)
        M2_batch = sumsq_new - (k * (mean_batch ** 2))

        n_total = self.count
        delta = mean_batch - mean_prev

        frac = np.divide(k, n_total, out=np.zeros_like(sum_new), where=(n_total > 0))
        self.mean = mean_prev + delta * frac

        cross = (delta ** 2) * np.where(n_total > 0, (count_prev * k) / n_total, 0.0)
        self.M2 = M2_prev + M2_batch + cross

    def finalize_std(self):
        """Compute sample standard deviation (ddof=1)."""
        sigma = np.full_like(self.mean, np.nan, dtype=np.float64)
        valid = self.count > 1
        sigma[valid] = np.sqrt(self.M2[valid] / (self.count[valid] - 1))
        return sigma


# ============================================================================
# CSV I/O helpers
# ============================================================================

def read_header_info(input_path: str) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Read the first two rows (position index and base type) from the CSV.

    Returns:
        (base_positions, base_types, total_cols, data_cols)
    """
    head = pd.read_csv(input_path, header=None, nrows=2,
                       na_values=["None"], skip_blank_lines=True)
    total_cols = head.shape[1]
    if total_cols < 5:
        raise ValueError("CSV must have at least 5 columns (4 metadata + data).")

    base_positions = head.iloc[0, 4:].astype(int).to_numpy()
    base_types = head.iloc[1, 4:].astype(str).to_numpy()
    data_cols = base_positions.shape[0]
    return base_positions, base_types, total_cols, data_cols


def gen_data_chunks(input_path: str, chunksize: int):
    """Yield data chunks starting from row 2 (after the two header rows)."""
    for chunk in pd.read_csv(
        input_path,
        header=None,
        skiprows=2,
        chunksize=chunksize,
        na_values=["None"],
        skip_blank_lines=True
    ):
        yield chunk


# ============================================================================
# QC mask
# ============================================================================

def apply_qc_mask(chunk: pd.DataFrame, data_cols: int, qc_bases: int,
                  base_positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Apply QC mask: exclude positions within qc_bases of alignment boundaries.

    Returns:
        (data_qc, meta_np, valid_mask):
            data_qc: (r, C) float32, out-of-range values set to NaN
            meta_np: (r, 4) string metadata
            valid_mask: (r, C) boolean mask
    """
    meta = chunk.iloc[:, 0:4]
    data_block = chunk.iloc[:, 4:]

    data_values = data_block.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)

    starts = meta.iloc[:, 1].astype(int).to_numpy().reshape(-1, 1)
    ends = meta.iloc[:, 2].astype(int).to_numpy().reshape(-1, 1)

    eff_starts = starts + qc_bases
    eff_ends = ends - qc_bases + 1

    pos_row = base_positions.reshape(1, -1)
    valid_mask = (pos_row >= eff_starts) & (pos_row < eff_ends)

    data_qc = np.where(valid_mask, data_values, np.nan).astype(np.float32)
    meta_np = meta.astype(str).to_numpy()

    return data_qc, meta_np, valid_mask


# ============================================================================
# Main computation
# ============================================================================

def compute_statistics(input_path: str, qc_bases: int, chunksize: int = 1000,
                       circular: str = "no") -> None:
    p = Path(input_path)
    if p.suffix.lower() != ".csv":
        raise ValueError("Only .csv files are supported.")

    base_positions, base_types, total_cols, data_cols = read_header_info(input_path)

    is_circular = (circular.lower() == "yes")
    if is_circular:
        genome_len = data_cols // 2
        if data_cols % 2 != 0:
            raise ValueError(f"circular=yes requires even data_cols ({data_cols}).")
        stats_cols = genome_len
    else:
        stats_cols = data_cols

    # Output file names
    stem = p.stem
    qc_tag = f"_QC{qc_bases}"
    out_qc_path = f"{stem}{qc_tag}.csv"
    out_mean_path = f"mean_{stem}{qc_tag}.csv"
    out_threshold_path = f"mean_2sigma_{stem}{qc_tag}.csv"

    # Single pass: Welford accumulation + QC CSV writing
    welford = WelfordVec(stats_cols)
    total_rows = 0
    empty_rows_after_qc = 0

    with open(out_qc_path, "w", encoding="utf-8") as f_qc:
        # Write header rows
        f_qc.write(",".join(["", "", "", ""] + [str(i) for i in base_positions.tolist()]) + "\n")
        f_qc.write(",".join(["", "", "", ""] + base_types.astype(str).tolist()) + "\n")

        for chunk in gen_data_chunks(input_path, chunksize):
            data_qc, meta_np, _ = apply_qc_mask(chunk, data_cols, qc_bases, base_positions)
            r = data_qc.shape[0]
            total_rows += r
            empty_rows_after_qc += int(np.sum(np.all(np.isnan(data_qc), axis=1)))

            # Write QC CSV rows
            for i in range(r):
                meta_str = [str(x) for x in meta_np[i, :].tolist()]
                data_str = [("nan" if np.isnan(v) else f"{v}") for v in data_qc[i, :].tolist()]
                f_qc.write(",".join(meta_str + data_str) + "\n")

            # Welford update
            if is_circular:
                data_circular = np.full((r * 2, genome_len), np.nan, dtype=np.float64)
                for j in range(genome_len):
                    data_circular[:r, j] = data_qc[:, j].astype(np.float64)
                    data_circular[r:, j] = data_qc[:, j + genome_len].astype(np.float64)
                welford.update(data_circular)
            else:
                welford.update(data_qc.astype(np.float64))

    # Finalize statistics
    count_by_pos = welford.count.astype(np.int64)
    mean_by_pos = welford.mean.astype(np.float64)
    sigma_by_pos = welford.finalize_std().astype(np.float64)

    mean_by_pos[count_by_pos == 0] = np.nan
    threshold_by_pos = mean_by_pos - 2 * sigma_by_pos
    threshold_by_pos[np.isnan(mean_by_pos) | np.isnan(sigma_by_pos)] = np.nan

    # Prepare output arrays (duplicate for circular genome)
    if is_circular:
        mean_out = np.concatenate([mean_by_pos, mean_by_pos])
        threshold_out = np.concatenate([threshold_by_pos, threshold_by_pos])
        count_out = np.concatenate([count_by_pos, count_by_pos])
    else:
        mean_out = mean_by_pos
        threshold_out = threshold_by_pos
        count_out = count_by_pos

    index_row = base_positions.tolist()
    base_row = base_types.astype(str).tolist()
    count_row = count_out.tolist()

    # Write mean CSV
    with open(out_mean_path, "w", encoding="utf-8") as f:
        f.write(",".join(map(str, mean_out.tolist())) + "\n")
        f.write(",".join(map(str, index_row)) + "\n")
        f.write(",".join(base_row) + "\n")
        f.write(",".join(map(str, count_row)) + "\n")

    # Write mean-2sigma threshold CSV
    with open(out_threshold_path, "w", encoding="utf-8") as f:
        f.write(",".join(map(str, threshold_out.tolist())) + "\n")
        f.write(",".join(map(str, index_row)) + "\n")
        f.write(",".join(base_row) + "\n")
        f.write(",".join(map(str, count_row)) + "\n")

    # Diagnostics
    print(f"[QC] rows total: {total_rows}, empty rows after QC: {empty_rows_after_qc}")
    print(f"[QC] empty columns: {int(np.sum(count_by_pos == 0))}/{data_cols}, "
          f"singleton columns: {int(np.sum(count_by_pos == 1))}")
    print(f"Output: {out_qc_path}, {out_mean_path}, {out_threshold_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute per-position mean and mean-2sigma threshold "
                    "from per-base signal CSV."
    )
    parser.add_argument("--input", required=True,
                        help="Input per-base signal CSV (output of extract_perbase_signal.py)")
    parser.add_argument("--qc-bases", type=int, default=0,
                        help="Number of bases to exclude from each alignment boundary (default: 0)")
    parser.add_argument("--chunk-size", type=int, default=1000,
                        help="Number of rows per processing chunk (default: 1000)")
    parser.add_argument("--circular", default="no", choices=["yes", "no"],
                        help="Circular genome mode for 2x-tiled references (default: no)")
    return parser.parse_args()


def main():
    args = parse_args()
    compute_statistics(args.input, args.qc_bases, args.chunk_size, args.circular)


if __name__ == "__main__":
    main()
