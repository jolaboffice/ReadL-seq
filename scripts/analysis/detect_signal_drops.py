#!/usr/bin/env python
"""
Per-position signal drop detection.

For each genomic position, counts the fraction of reads whose signal falls
below a pre-computed threshold (mu - 2*sigma), considering only positions
within the aligned region (excluding soft-clipped bases).
"""
import argparse

import numpy as np
import pandas as pd


def load_threshold(threshold_path):
    """
    Load pre-computed threshold values from CSV (output of mean_sigma.py).

    The threshold CSV has 4 rows: values, positions, bases, counts.

    Args:
        threshold_path: Path to threshold CSV file.

    Returns:
        Series with position index and threshold values.
    """
    df = pd.read_csv(threshold_path, header=None)
    threshold = df.iloc[0, :].astype(float)
    threshold.index = df.iloc[1, :].astype(int)
    return threshold


def detect_signal_drops(signal_path, threshold_path, chunksize=1000):
    """
    Detect per-position signal drops below threshold.

    For each read, only positions within the aligned region are considered.

    Args:
        signal_path: Per-base signal CSV (output of extract_perbase_signal.py).
        threshold_path: Pre-computed threshold CSV (output of mean_sigma.py).
        chunksize: Number of rows per processing chunk.
    """
    # Read header info
    df_header = pd.read_csv(signal_path, nrows=2, header=None)
    data_columns = df_header.iloc[0, 4:].astype(int)
    base = df_header.iloc[1, 4:]

    # Load threshold
    threshold = load_threshold(threshold_path)
    threshold_aligned = threshold.reindex(data_columns)

    # Accumulate counts in chunks
    biotin_signal_counts = pd.Series(0, index=data_columns, dtype=float)
    total_reads_per_base = pd.Series(0, index=data_columns, dtype=int)

    print("Processing data in chunks...")
    chunk_iter = pd.read_csv(
        signal_path,
        skiprows=2,
        header=None,
        chunksize=chunksize,
        na_values=["None"]
    )

    for chunk_num, chunk in enumerate(chunk_iter):
        if chunk_num % 10 == 0:
            print(f"Processing chunk {chunk_num + 1}...")

        chunk_read_info = chunk.iloc[:, 0:4]
        chunk_data = chunk.iloc[:, 4:]
        n_cols = min(len(chunk_data.columns), len(data_columns))
        chunk_data = chunk_data.iloc[:, :n_cols]
        chunk_data.columns = data_columns[:n_cols]
        chunk_data = chunk_data.apply(pd.to_numeric, errors="coerce")

        # Exclude soft-clipped regions: keep only aligned positions per read
        start_info = chunk_read_info.iloc[:, 1].astype(int).values
        end_info = chunk_read_info.iloc[:, 2].astype(int).values
        col_positions = chunk_data.columns.astype(int).values
        start_2d = start_info[:, np.newaxis]
        end_2d = end_info[:, np.newaxis]
        col_2d = col_positions[np.newaxis, :]
        mask_outside = (col_2d < start_2d) | (col_2d > end_2d)
        chunk_values = chunk_data.values
        chunk_values[mask_outside] = np.nan
        chunk_data = pd.DataFrame(chunk_values, index=chunk_data.index,
                                  columns=chunk_data.columns)

        # Compare against threshold
        common_cols = chunk_data.columns.intersection(threshold_aligned.index)
        if len(common_cols) > 0:
            chunk_aligned = chunk_data[common_cols]
            threshold_subset = threshold_aligned[common_cols]
            chunk_float = chunk_aligned.astype(float)
            drop_mask = (chunk_float < threshold_subset)
            biotin_signal_counts[common_cols] += drop_mask.astype(float).sum(skipna=False)
            total_reads_per_base[common_cols] += chunk_aligned.count()

    print("Processing complete.")
    ratio = biotin_signal_counts / total_reads_per_base

    # Align base info
    base_values = base.values
    n_final = len(data_columns)
    if len(base_values) >= n_final:
        base_aligned = list(base_values[:n_final])
    else:
        base_aligned = list(base_values) + [None] * (n_final - len(base_values))

    # Build result table
    total_aligned = total_reads_per_base.reindex(data_columns)
    drops_aligned = biotin_signal_counts.reindex(data_columns)
    ratio_aligned = ratio.reindex(data_columns)
    threshold_final = threshold_aligned.reindex(data_columns)

    result = pd.DataFrame({
        'pos': data_columns.values,
        'base': base_aligned,
        'total_reads': total_aligned.values,
        'biotin_signal_reads': drops_aligned.values,
        'ratio': ratio_aligned.values,
        'threshold': threshold_final.values
    })

    # Overall damage ratio
    total_bases = result.total_reads.sum()
    total_drops = result.biotin_signal_reads.sum()
    damage_ratio = total_drops / total_bases
    print(f"Damage ratio: {damage_ratio:.3e}")

    # Save output
    output_name = f"signal_drop_ratio_{damage_ratio:.3e}.csv".replace("+0", "").replace("+", "")
    result.to_csv(output_name, index=False)
    print(f"Output: {output_name}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect per-position signal drops below threshold (mu - 2*sigma)."
    )
    parser.add_argument("--signal", required=True,
                        help="Per-base signal CSV (output of extract_perbase_signal.py)")
    parser.add_argument("--threshold", required=True,
                        help="Threshold CSV (output of mean_sigma.py)")
    parser.add_argument("--chunk-size", type=int, default=1000,
                        help="Rows per processing chunk (default: 1000)")
    return parser.parse_args()


def main():
    args = parse_args()
    detect_signal_drops(args.signal, args.threshold, args.chunk_size)


if __name__ == "__main__":
    main()
