#!/usr/bin/env python
"""
Compute delta signal (mu - signal) for each read at each genomic position.

Reads a per-base signal CSV and a mean CSV, and outputs a CSV of the same
format where each data value is replaced by the deviation from the mean.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_mean(mean_path):
    """
    Load per-position mean values from CSV (output of mean_sigma.py).

    The mean CSV has 4 rows: values, positions, bases, counts.

    Args:
        mean_path: Path to mean CSV file.

    Returns:
        Series with position index and mean values.
    """
    df = pd.read_csv(mean_path, header=None)
    mean = df.iloc[0, :].astype(float)
    mean.index = df.iloc[1, :].astype(int)
    return mean


def compute_delta_signal(signal_path, mean_path, output_path, chunksize=2000):
    """
    Compute delta signal (mu - signal) and write to output CSV.

    The first 2 rows (position index and base type) are copied as-is.
    Columns 0-3 (metadata) are preserved; columns 4+ are replaced with
    mu - signal values.

    Args:
        signal_path: Per-base signal CSV (output of extract_perbase_signal.py).
        mean_path: Mean CSV (output of mean_sigma.py).
        output_path: Output CSV path.
        chunksize: Number of rows per processing chunk.
    """
    # Read header and position array
    head = pd.read_csv(signal_path, header=None, nrows=2)
    pos = head.iloc[0, 4:].astype(int).to_numpy()

    # Load mean and align to position order
    mean_series = load_mean(mean_path)
    mean_vec = mean_series.reindex(pos).to_numpy(dtype=np.float64)

    # Write header rows
    head.to_csv(output_path, index=False, header=False, mode='w')

    # Count total data rows for progress reporting
    total_rows = 0
    with open(signal_path, 'r', encoding='utf-8', newline='') as f:
        for _ in f:
            total_rows += 1
    data_rows = max(total_rows - 2, 0)

    # Process data in chunks
    processed = 0
    reader = pd.read_csv(
        signal_path,
        header=None,
        skiprows=2,
        chunksize=chunksize,
        na_values=["None"],
        keep_default_na=True
    )

    for chunk in reader:
        meta = chunk.iloc[:, :4]
        data = chunk.iloc[:, 4:]
        data = data.apply(pd.to_numeric, errors="coerce")
        data_vals = data.to_numpy(dtype=np.float64)

        # Broadcast: mu - signal (NaN propagates automatically)
        delta_vals = mean_vec - data_vals

        out_chunk = pd.concat(
            [meta.reset_index(drop=True),
             pd.DataFrame(delta_vals, columns=pos)], axis=1
        )
        out_chunk.to_csv(output_path, index=False, header=False, mode='a')

        processed += len(chunk)
        percent = (processed / data_rows * 100) if data_rows else 100.0
        print(f"Processed rows: {processed}/{data_rows} ({percent:.1f}%)", end="\r")

    print("\nDone.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute delta signal (mu - signal) for each read at each position."
    )
    parser.add_argument("--signal", required=True,
                        help="Per-base signal CSV (output of extract_perbase_signal.py)")
    parser.add_argument("--mean", required=True,
                        help="Mean CSV (output of mean_sigma.py)")
    parser.add_argument("--output", default=None,
                        help="Output CSV path (default: <signal_stem>_delta.csv)")
    parser.add_argument("--chunk-size", type=int, default=2000,
                        help="Rows per processing chunk (default: 2000)")
    return parser.parse_args()


def main():
    args = parse_args()
    output_path = args.output
    if output_path is None:
        p = Path(args.signal)
        output_path = str(p.with_name(p.stem + "_delta.csv"))

    compute_delta_signal(args.signal, args.mean, output_path, args.chunk_size)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
