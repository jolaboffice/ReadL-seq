#!/usr/bin/env python
"""
BiotinScore calculation and biotin-positive read detection.

For each read, computes BiotinScore = sum of delta-signal values in a window
around the 3' aligned end, and classifies reads as biotin-positive based on
a score threshold. Outputs per-read and per-position summary tables.

Usage:
  python biotin_scoring.py <input.csv> -b <BAM> [-o OUT_XLSX]
      [--biotin-score-threshold THRESH] [--circular {no,yes}]
"""
import argparse
import os

import numpy as np
import pandas as pd
import pysam

# ======================== Parameters ========================
CHUNKSIZE = 2000

# biotin_score parameters
# Scoring window: [anchor - START_OFFSET, anchor + min(END_OFFSET, sc_len))
BIOTIN_SCORE_START_OFFSET = 3
BIOTIN_SCORE_END_OFFSET = 100
BIOTIN_SCORE_THRESH = 6.0

# total_reads and biotin hit assignment range
TOTAL_READS_SOFTCLIP_LIMIT = 10
BIOTIN_POSITIONS_RANGE = 10
# ============================================================

META_COLS = ['read_id', 'start', 'end', 'length']


def biotin_score(signal, anchor, window=100, sc_len=0):
    """
    Compute BiotinScore for a single read.

    Sums delta-signal values in [anchor - START_OFFSET, anchor + effective_window),
    where effective_window = min(window, sc_len).

    Args:
        signal: 1D numpy array (delta-signal with NaN replaced by 0).
        anchor: Index in signal array corresponding to the 3' aligned end.
        window: Maximum window size beyond anchor (default 100).
        sc_len: 3' soft-clip length; effective window is min(window, sc_len).

    Returns:
        float: BiotinScore (sum of delta-signal values in the range).
    """
    effective_window = min(window, sc_len + 1)
    start = max(anchor - BIOTIN_SCORE_START_OFFSET, 0)
    end = min(anchor + effective_window, len(signal))
    return float(np.sum(signal[start:end]))


def detect_biotin_runs(input_csv_path, bam_path, circular="no",
                       biotin_score_thresh=None):
    """
    Detect biotin-positive reads using BiotinScore.

    Args:
        input_csv_path: Path to delta-signal CSV file.
        bam_path: Path to BAM file (for 3' soft-clip length).
        circular: "no" or "yes" for circular genome mode.
        biotin_score_thresh: BiotinScore threshold for biotin classification.

    Returns:
        tuple: (df_per_read, df_per_position, base_row_full, position_axis, cols)
    """
    if biotin_score_thresh is None:
        biotin_score_thresh = BIOTIN_SCORE_THRESH

    header_df = pd.read_csv(input_csv_path, nrows=0)
    cols = header_df.columns.tolist()
    if not all(c in cols for c in META_COLS):
        raise ValueError(f"CSV must contain columns: {META_COLS}")
    position_columns = [c for c in cols if c not in META_COLS]

    # Row 1: reference base row
    base_row_full = pd.read_csv(
        input_csv_path, header=None, names=cols, skiprows=1, nrows=1
    ).iloc[0]

    # Position axis (column names are integer positions)
    try:
        position_axis = np.array([int(p) for p in position_columns], dtype=int)
    except Exception:
        raise ValueError("Position column names must be integer strings.")

    n_pos = len(position_axis)

    # Circular genome handling
    is_circular = (circular.lower() == "yes")
    if is_circular:
        genome_len = len(position_axis) // 2
        if len(position_axis) % 2 != 0:
            raise ValueError(
                f"circular=yes requires even number of positions ({len(position_axis)}).")
        stats_position_axis = position_axis[:genome_len]
    else:
        stats_position_axis = position_axis
        genome_len = None

    # Load 3' soft-clip lengths from BAM
    print("[BiotinScore] Reading BAM for 3' soft-clip lengths...")
    sc_3prime_info = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped:
            continue
        cig = r.cigartuples
        if not cig:
            continue
        # Forward strand: 3' soft-clip is the trailing S operation
        sc_3prime_len = int(cig[-1][1]) if cig[-1][0] == 4 else 0
        sc_3prime_info[r.query_name] = sc_3prime_len
    bam.close()
    print(f"[BiotinScore] Loaded 3' soft-clip info for {len(sc_3prime_info)} reads")

    # Accumulators
    read_ids_all, starts_all, ends_all, lengths_all = [], [], [], []
    biotin_scores_all = []
    biotin_positions_all = []
    is_biotin_all = []
    softclip_3prime_lengths = []
    end_plus_softclip = []

    # Per-position counting arrays
    total_counts = np.zeros(n_pos, dtype=np.int64)
    biotin_counts = np.zeros(n_pos, dtype=np.int64)

    processed = 0
    reader = pd.read_csv(
        input_csv_path,
        header=None,
        names=cols,
        skiprows=2,
        chunksize=CHUNKSIZE,
        na_values=["None", ""],
        keep_default_na=True,
        engine="python",
        on_bad_lines="error"
    )

    for chunk in reader:
        starts = pd.to_numeric(chunk['start'], errors='raise')
        ends = pd.to_numeric(chunk['end'], errors='raise')
        lengths = pd.to_numeric(chunk['length'], errors='coerce')

        data_vals = chunk[position_columns].apply(
            pd.to_numeric, errors="coerce"
        ).to_numpy(dtype=float)

        for i in range(len(chunk)):
            row = data_vals[i, :]
            read_id = chunk['read_id'].iat[i]
            read_start = int(starts.iat[i])
            read_end = int(ends.iat[i])
            read_len = (int(lengths.iat[i]) if not np.isnan(lengths.iat[i])
                        else (read_end - read_start))

            # 3' soft-clip length
            sc_3p_len = sc_3prime_info.get(str(read_id), 0)
            end_plus_sc = read_end + sc_3p_len

            # Effective end: end + min(softclip_3prime, TOTAL_READS_SOFTCLIP_LIMIT)
            real_end = read_end + min(sc_3p_len, TOTAL_READS_SOFTCLIP_LIMIT)

            # Count total reads per position: [start, real_end)
            si = np.searchsorted(position_axis, read_start, side='left')
            ei = np.searchsorted(position_axis, real_end, side='right')
            if si < ei:
                total_counts[si:ei] += 1

            # Find anchor index for read_end
            anchor_idx = np.searchsorted(position_axis, read_end, side='left')
            if anchor_idx >= n_pos or position_axis[anchor_idx] != read_end:
                anchor_idx = max(0, min(anchor_idx, n_pos - 1))

            signal_array = np.nan_to_num(row, nan=0.0)

            # Compute BiotinScore
            score = biotin_score(
                signal_array,
                anchor_idx,
                window=BIOTIN_SCORE_END_OFFSET,
                sc_len=sc_3p_len,
            )

            # Classify as biotin-positive
            is_biotin = (score >= biotin_score_thresh)

            # Biotin hit assignment and position reporting
            if is_biotin:
                # Per-position biotin hits: [end - 10, real_end)
                biotin_range_start = read_end - BIOTIN_POSITIONS_RANGE
                bi = np.searchsorted(position_axis, biotin_range_start, side='left')
                be = np.searchsorted(position_axis, real_end, side='right')
                if bi < be:
                    biotin_counts[bi:be] += 1

                # Per-read biotin positions: positive signal within end +/- 10
                analysis_start = read_end - BIOTIN_POSITIONS_RANGE
                analysis_end = read_end + BIOTIN_POSITIONS_RANGE

                a_si = np.searchsorted(position_axis, analysis_start, side='left')
                a_ei = np.searchsorted(position_axis, analysis_end, side='right')
                if a_si < a_ei:
                    sub_signal = signal_array[a_si:a_ei]
                    sub_finite = np.isfinite(row[a_si:a_ei])
                    bp_mask = (sub_signal > 0) & sub_finite
                    biotin_positions = position_axis[a_si:a_ei][bp_mask].tolist()
                else:
                    biotin_positions = []
            else:
                biotin_positions = []

            # Store results
            read_ids_all.append(read_id)
            starts_all.append(read_start)
            ends_all.append(read_end)
            lengths_all.append(read_len)
            biotin_scores_all.append(score)
            biotin_positions_all.append(sorted(biotin_positions))
            is_biotin_all.append(is_biotin)
            softclip_3prime_lengths.append(sc_3p_len)
            end_plus_softclip.append(end_plus_sc)

        processed += len(chunk)
        print(f"[BiotinScore] processed rows: {processed}", end="\r")
    print()

    # Per-read DataFrame
    df_per_read = pd.DataFrame({
        'read_id': read_ids_all,
        'start': starts_all,
        'end': ends_all,
        'length': lengths_all,
        'softclip_3prime_length': softclip_3prime_lengths,
        'end_plus_softclip': end_plus_softclip,
        'biotin_score': biotin_scores_all,
        'is_biotin': is_biotin_all,
        'biotin_count': [len(peaks) for peaks in biotin_positions_all],
        'biotin_positions': [
            ', '.join(map(str, peaks)) if peaks else ''
            for peaks in biotin_positions_all
        ]
    })

    # Per-position DataFrame
    if is_circular:
        total_arr = total_counts[:genome_len] + total_counts[genome_len:]
        biotin_arr = biotin_counts[:genome_len] + biotin_counts[genome_len:]
        positions = stats_position_axis
        bases = np.array([base_row_full.get(str(int(p)), '-') for p in positions])
    else:
        total_arr = total_counts
        biotin_arr = biotin_counts
        positions = position_axis
        bases = np.array([base_row_full.get(str(int(p)), '-') for p in positions])

    ratio = np.divide(biotin_arr.astype(float), total_arr.astype(float),
                      out=np.zeros(len(positions), dtype=float),
                      where=(total_arr > 0))
    df_per_position = pd.DataFrame({
        'position': positions,
        'base': bases,
        'total_reads': total_arr,
        'biotin_hits': biotin_arr,
        'biotin_ratio': ratio
    })

    return df_per_read, df_per_position, base_row_full, position_axis, cols


def main():
    ap = argparse.ArgumentParser(
        description="BiotinScore calculation and biotin-positive read detection.")
    ap.add_argument("input_csv",
                    help="Delta-signal CSV (row 0: header, row 1: base row)")
    ap.add_argument("-b", "--bam", required=True,
                    help="BAM file path (for 3' soft-clip length)")
    ap.add_argument("-o", "--out", default=None,
                    help="Output .xlsx path")
    ap.add_argument("--circular", default="no", choices=["no", "yes"],
                    help="Circular genome mode (default: no)")
    ap.add_argument("--biotin-score-threshold", type=float,
                    default=BIOTIN_SCORE_THRESH,
                    help=f"BiotinScore threshold (default: {BIOTIN_SCORE_THRESH})")
    args = ap.parse_args()

    df_per_read, df_per_position, _, _, _ = detect_biotin_runs(
        args.input_csv,
        bam_path=args.bam,
        circular=args.circular,
        biotin_score_thresh=args.biotin_score_threshold,
    )

    # Output path
    if args.out:
        out_xlsx = args.out
    else:
        basename = os.path.splitext(os.path.basename(args.input_csv))[0]
        input_dir = os.path.dirname(args.input_csv)
        thresh_val = (int(args.biotin_score_threshold)
                      if args.biotin_score_threshold == int(args.biotin_score_threshold)
                      else args.biotin_score_threshold)
        out_xlsx = os.path.join(input_dir, f"{basename}_BiotinScore_{thresh_val}.xlsx")

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as w:
        df_per_read.to_excel(w, index=False, sheet_name="per_read")
        df_per_position.to_excel(w, index=False, sheet_name="per_position")
    print(f"[Done] Saved: {out_xlsx}")


if __name__ == "__main__":
    main()
