#!/usr/bin/env python
"""
Per-base nanopore signal extraction with soft-clip support.

Extracts per-base representative signal values from nanopore sequencing data
using Remora's signal-to-reference mapping, extended to soft-clipped regions.
"""
import argparse
import os
import random
import sys

import numpy as np
import pandas as pd
import pod5
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from remora import io, refine_signal_map
from remora import data_chunks as DC
from remora.io import compute_base_space_sig_coords

CIGAR_SOFT_CLIP = 4


# ============================================================================
# Utility functions
# ============================================================================

def get_seq(reference, strand, start, end):
    """
    Get sequence from reference file.

    Args:
        reference: Reference file path.
        strand: Strand of the sequence ("+" or "-").
        start: Start position (0-based).
        end: End position (0-based, exclusive).

    Returns:
        Sequence string (upper case, 5' to 3').
    """
    for record in SeqIO.parse(reference, "fasta"):
        sequence = str(record.seq).upper()
        break  # Use only the first record

    seq_slice = sequence[start:end]

    if strand == "-" or (isinstance(strand, str) and strand.lower() == "reverse"):
        seq_slice = str(Seq(seq_slice).reverse_complement())

    return seq_slice


def get_softclip_lengths(cigar_tuples):
    """
    Extract soft-clip lengths from CIGAR tuples.

    Args:
        cigar_tuples: List of CIGAR tuples.

    Returns:
        (sc_5prime, sc_3prime): 5' and 3' soft-clip lengths.
    """
    sc_5prime = 0
    sc_3prime = 0

    if cigar_tuples and len(cigar_tuples) > 0:
        if cigar_tuples[0][0] == CIGAR_SOFT_CLIP:
            sc_5prime = cigar_tuples[0][1]
        if cigar_tuples[-1][0] == CIGAR_SOFT_CLIP:
            sc_3prime = cigar_tuples[-1][1]

    return sc_5prime, sc_3prime


def refine_ref_to_signal(ref_to_signal, io_read, sig_map_refiner, read_id, seq, sig):
    """
    Refine ref_to_signal mapping using Remora's SigMapRefiner.

    Args:
        ref_to_signal: Reference-to-signal mapping array (before refinement).
        io_read: IORead object.
        sig_map_refiner: SigMapRefiner object.
        read_id: Read ID string.
        seq: Reference sequence for refinement.
        sig: Fallback signal array.

    Returns:
        (refined_ref_to_signal, normalized_sig). Returns (ref_to_signal, sig) on failure.
    """
    if sig_map_refiner is None or not sig_map_refiner.is_loaded:
        return ref_to_signal, sig

    if ref_to_signal is None or len(ref_to_signal) == 0:
        return ref_to_signal, sig if sig is not None else np.array([])

    if len(ref_to_signal) == 1:
        return ref_to_signal, np.array([])

    try:
        sig_start = ref_to_signal[0]
        sig_end = ref_to_signal[-1]

        if sig_start >= sig_end:
            return ref_to_signal, np.array([])

        trim_dacs = io_read.dacs[sig_start:sig_end]
        shift_seq_to_sig = ref_to_signal - sig_start

        remora_read = DC.RemoraRead(
            dacs=trim_dacs,
            shift=io_read.shift_dacs_to_norm,
            scale=io_read.scale_dacs_to_norm,
            seq_to_sig_map=shift_seq_to_sig,
            str_seq=seq,
            read_id=read_id,
        )

        remora_read.refine_signal_mapping(sig_map_refiner)

        refined_ref_to_signal = remora_read.seq_to_sig_map + sig_start
        normalized_sig = (trim_dacs - remora_read.shift) / remora_read.scale

        return refined_ref_to_signal, normalized_sig

    except Exception as e:
        print(f"Refinement failed ({read_id}): {e}")
        return ref_to_signal, sig


# ============================================================================
# Signal extraction functions
# ============================================================================

def extract_matched_region_signal(ref_to_signal, read_reg, ref_reg_io):
    """
    Extract signal mapping for the matched (aligned) region.

    Args:
        ref_to_signal: Full reference-to-signal mapping.
        read_reg: Read region object.
        ref_reg_io: Reference region from IORead.

    Returns:
        (ref_to_signal_region, seq_to_sig_map):
            Region mapping and zero-based adjusted mapping.
    """
    if ref_reg_io.strand == "+":
        reg_st = max(0, read_reg.ref_reg.start - ref_reg_io.start)
        reg_en = read_reg.ref_reg.end - ref_reg_io.start
    else:
        reg_st = max(0, ref_reg_io.end - read_reg.ref_reg.end)
        reg_en = ref_reg_io.end - read_reg.ref_reg.start

    if reg_st >= 0 and reg_en < len(ref_to_signal):
        ref_to_signal = ref_to_signal[reg_st:reg_en + 1].copy()

        if len(ref_to_signal) == 0:
            seq_to_sig_map = np.array([])
        elif len(ref_to_signal) == 1:
            seq_to_sig_map = np.array([0])
        else:
            seq_to_sig_map = ref_to_signal - ref_to_signal[0]
        return ref_to_signal, seq_to_sig_map
    else:
        return np.array([]), np.array([])


def extract_softclip_signals(sig_full, sig_start_idx, sig_end_idx):
    """
    Extract signal segments for 5' and 3' soft-clip regions.

    Args:
        sig_full: Full normalized signal array.
        sig_start_idx: Signal index at the start of the matched region.
        sig_end_idx: Signal index at the end of the matched region.

    Returns:
        (sig_sc_5prime, sig_sc_3prime): 5' and 3' soft-clip signal arrays.
    """
    sig_sc_5prime = sig_full[0:sig_start_idx]
    sig_sc_3prime = sig_full[sig_end_idx:]
    return sig_sc_5prime, sig_sc_3prime


def process_softclip_xcoord(query_to_signal, sc_len, io_read, is_5prime, sig_sc=None):
    """
    Compute signal mapping for a soft-clip region using the query-to-signal array.

    The query-to-signal mapping spans the full basecalled query, including
    soft-clipped bases. The prefix and suffix entries, as defined by the
    5' and 3' soft-clip lengths, directly provide signal-index boundaries.

    Args:
        query_to_signal: Query-to-signal mapping array.
        sc_len: Soft-clip length (from CIGAR string).
        io_read: IORead object.
        is_5prime: True for 5' soft-clip, False for 3' soft-clip.
        sig_sc: Soft-clip signal array (for length validation).

    Returns:
        (ref_to_signal_sc, seq_to_sig_map_sc, sig_sc_out):
            Signal mapping, zero-based mapping, and normalized signal
            for the soft-clip region.
    """
    if sc_len > 0 and (sig_sc is None or len(sig_sc) > 0):
        if is_5prime:
            ref_to_signal_sc = query_to_signal[0:sc_len + 1]
        else:
            ref_to_signal_sc = query_to_signal[-sc_len - 1:]

        if len(ref_to_signal_sc) > 0:
            sc_sig_start = ref_to_signal_sc[0]
            sc_sig_end = ref_to_signal_sc[-1]
            sig_sc_out = io_read.norm_signal[sc_sig_start:sc_sig_end]
        else:
            sig_sc_out = sig_sc if sig_sc is not None else np.array([])

        if len(ref_to_signal_sc) == 0:
            return ref_to_signal_sc, np.array([]), np.array([])
        elif len(ref_to_signal_sc) == 1:
            seq_to_sig_map_sc = np.array([0])
            return ref_to_signal_sc, seq_to_sig_map_sc, sig_sc if sig_sc is not None else np.array([])

        seq_to_sig_map_sc = ref_to_signal_sc - ref_to_signal_sc[0]
        return ref_to_signal_sc, seq_to_sig_map_sc, sig_sc_out
    else:
        return np.array([]), np.array([]), np.array([])


def process_reverse_strand(sig, sig_sc_5prime, sig_sc_3prime,
                           seq_to_sig_map, seq_to_sig_map_sc_5prime,
                           seq_to_sig_map_sc_3prime):
    """
    Reverse signal and coordinate arrays for reverse-strand reads.

    Signals are flipped and 5'/3' soft-clip assignments are swapped.

    Args:
        sig: Matched region signal.
        sig_sc_5prime: 5' soft-clip signal.
        sig_sc_3prime: 3' soft-clip signal.
        seq_to_sig_map: Matched region mapping.
        seq_to_sig_map_sc_5prime: 5' soft-clip mapping.
        seq_to_sig_map_sc_3prime: 3' soft-clip mapping.

    Returns:
        (sig, sig_left, sig_right, seq_to_sig_map, map_left, map_right):
            Reversed and swapped arrays.
    """
    sig = sig[::-1]
    sig_sc_5prime = sig_sc_5prime[::-1]
    sig_sc_3prime = sig_sc_3prime[::-1]

    seq_to_sig_map = seq_to_sig_map[-1] - seq_to_sig_map[::-1]
    if len(seq_to_sig_map_sc_5prime) > 0:
        seq_to_sig_map_sc_5prime = seq_to_sig_map_sc_5prime[-1] - seq_to_sig_map_sc_5prime[::-1]
    if len(seq_to_sig_map_sc_3prime) > 0:
        seq_to_sig_map_sc_3prime = seq_to_sig_map_sc_3prime[-1] - seq_to_sig_map_sc_3prime[::-1]

    return (sig, sig_sc_3prime, sig_sc_5prime,
            seq_to_sig_map, seq_to_sig_map_sc_3prime, seq_to_sig_map_sc_5prime)


def calculate_xcoords(seq_to_sig_map, seq_to_sig_map_left, seq_to_sig_map_right,
                      read_reg, ref_reg_io):
    """
    Calculate reference-space coordinates for matched and soft-clip regions.

    Args:
        seq_to_sig_map: Matched region mapping.
        seq_to_sig_map_left: Left soft-clip mapping.
        seq_to_sig_map_right: Right soft-clip mapping.
        read_reg: Read region object.
        ref_reg_io: Reference region object.

    Returns:
        (xCoord, xCoord_left, xCoord_right): Reference-space coordinates.
    """
    xCoord_base = compute_base_space_sig_coords(seq_to_sig_map)
    xCoord = xCoord_base + read_reg.ref_reg.start

    if len(seq_to_sig_map_left) > 0:
        xCoord_left = compute_base_space_sig_coords(seq_to_sig_map_left)
        xCoord_left += read_reg.ref_reg.start - len(seq_to_sig_map_left) + 1
    else:
        xCoord_left = np.array([])

    if len(seq_to_sig_map_right) > 0:
        xCoord_right = compute_base_space_sig_coords(seq_to_sig_map_right)
        xCoord_right += read_reg.ref_reg.end
    else:
        xCoord_right = np.array([])

    return xCoord, xCoord_left, xCoord_right


# ============================================================================
# Main
# ============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract per-base nanopore signal values with soft-clip support."
    )
    parser.add_argument("--bam", required=True, help="Sorted BAM file path")
    parser.add_argument("--pod5", required=True, help="POD5 file path")
    parser.add_argument("--reference", required=True, help="Reference FASTA file path")
    parser.add_argument("--level-table", required=True,
                        help="k-mer expected level table path for SigMapRefiner")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    parser.add_argument("--strand", default="forward", choices=["forward", "reverse"],
                        help="Strand to analyze (default: forward)")
    parser.add_argument("--target-start", type=int, default=None,
                        help="Target start position, 1-based (default: 1)")
    parser.add_argument("--target-end", type=int, default=None,
                        help="Target end position, 1-based (default: reference length)")
    parser.add_argument("--chunk-size", type=int, default=500,
                        help="Number of rows per CSV write chunk (default: 500)")
    parser.add_argument("--max-reads", type=int, default=None,
                        help="Maximum number of reads to process (default: all)")
    return parser.parse_args()


def main():
    args = parse_args()

    # Read reference sequence
    for record in SeqIO.parse(args.reference, "fasta"):
        sequence = str(record.seq).upper()
        ctg_name = record.id
        ref_length = len(sequence)

    # Set analysis range (1-based positions)
    target_start = args.target_start if args.target_start is not None else 1
    target_end = args.target_end if args.target_end is not None else ref_length

    position_columns = list(range(target_start, target_end + 1))
    columns = ['read_id', 'start', 'end', 'length'] + position_columns
    base_row = ['base', 'start', 'end', 'length'] + [sequence[pos - 1] for pos in position_columns]

    # Create output directory and write header
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    pd.DataFrame([base_row], columns=columns).to_csv(args.output, index=False)

    # Initialize signal map refiner
    print("Initializing SigMapRefiner...")
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=args.level_table,
        scale_iters=0,
        do_fix_guage=True,
    )
    print("Done.\n")

    # Set reference region
    direction = "+" if args.strand == "forward" else "-"
    ref_reg = io.RefRegion(ctg=ctg_name, strand=direction, start=0, end=ref_length)
    pod5_fh = pod5.Reader(args.pod5)
    bam_fh = pysam.AlignmentFile(args.bam)

    # Fetch BAM reads covering the reference region
    print("Fetching read signal mappings...")
    sample_bam_reads = io.get_reg_bam_reads(ref_reg, bam_fh)
    if len(sample_bam_reads) == 0:
        raise ValueError("No reads covering region")
    if args.max_reads is not None and len(sample_bam_reads) > args.max_reads:
        sample_bam_reads = random.sample(sample_bam_reads, args.max_reads)

    # Create IORead objects
    io_reads, skipped_count = io.get_io_reads(
        sample_bam_reads, pod5_fh, reverse_signal=False, missing_ok=True
    )
    print(f"Skipped {skipped_count} reads due to move table/signal mismatch or missing POD5 records")

    # Extract reference regions
    read_ref_regs = [[io_read.extract_ref_reg(ref_reg) for io_read in io_reads]]

    # Build lookup dictionaries
    print("Building IORead lookup table...")
    io_read_dict = {io_read.read_id: io_read for io_read in io_reads}
    print(f"Loaded IORead info for {len(io_read_dict)} reads")

    print("Building read CIGAR lookup table...")
    read_cigar_dict = {}
    for bam_read in sample_bam_reads:
        if bam_read.is_unmapped:
            continue
        read_cigar_dict[bam_read.query_name] = bam_read.cigartuples
    print(f"Loaded CIGAR info for {len(read_cigar_dict)} reads")

    # Signal extraction
    print("Extracting signals...")
    count = 0
    all_rows = []
    for sample_reads in read_ref_regs:
        total_reads = len(sample_reads)
        for read_reg in sample_reads:
            count += 1
            read_id = read_reg.read_id
            matched_start_base = read_reg.ref_reg.start  # 0-based
            matched_end_base = read_reg.ref_reg.end      # 0-based exclusive

            # Get full signal and query_to_signal from IORead
            sig_full = None
            query_to_signal = None
            ref_reg_io = None
            if read_id in io_read_dict:
                io_read = io_read_dict[read_id]
                sig_full = io_read.norm_signal
                query_to_signal = io_read.query_to_signal
                ref_reg_io = io_read.ref_reg
            else:
                sig_full = read_reg.norm_signal

            # Get soft-clip lengths from CIGAR
            soft_clip_5prime_len = 0
            soft_clip_3prime_len = 0
            cigar = None
            if read_id in read_cigar_dict:
                cigar = read_cigar_dict[read_id]
                if cigar:
                    soft_clip_5prime_len, soft_clip_3prime_len = get_softclip_lengths(cigar)

            sig = None
            xCoord_full = None
            signal_final = None

            ref_to_signal = None
            if query_to_signal is not None and cigar is not None and ref_reg_io is not None:
                # Build reference-to-signal mapping
                knots = DC.make_sequence_coordinate_mapping2(cigar)
                ref_to_signal = DC.map_ref_to_signal(
                    query_to_signal=query_to_signal,
                    ref_to_query_knots=knots
                )

                # Refine mapping
                match_start = ref_reg_io.start
                match_end = ref_reg_io.end
                matched_seq = get_seq(
                    args.reference, ref_reg_io.strand, match_start, match_end
                )

                if len(ref_to_signal) > 0:
                    sig_fallback = sig_full[ref_to_signal[0]:ref_to_signal[-1]]
                else:
                    sig_fallback = np.array([])

                ref_to_signal, sig = refine_ref_to_signal(
                    ref_to_signal=ref_to_signal,
                    io_read=io_read,
                    sig_map_refiner=sig_map_refiner,
                    read_id=read_id,
                    seq=matched_seq,
                    sig=sig_fallback
                )

                # Extract matched region signal
                ref_to_signal_region, seq_to_sig_map = extract_matched_region_signal(
                    ref_to_signal, read_reg, ref_reg_io
                )
                ref_to_signal = ref_to_signal_region

                if len(sig) > 0 and len(ref_to_signal) > 0:
                    sig_start_idx = ref_to_signal[0]
                    sig_end_idx = ref_to_signal[-1]

                    # Extract soft-clip signals
                    sig_sc_5prime_raw, sig_sc_3prime_raw = extract_softclip_signals(
                        sig_full, sig_start_idx, sig_end_idx
                    )

                    # Compute soft-clip mappings
                    _, seq_to_sig_map_sc_5prime, sig_sc_5prime = process_softclip_xcoord(
                        query_to_signal, soft_clip_5prime_len, io_read,
                        is_5prime=True, sig_sc=sig_sc_5prime_raw
                    )
                    _, seq_to_sig_map_sc_3prime, sig_sc_3prime = process_softclip_xcoord(
                        query_to_signal, soft_clip_3prime_len, io_read,
                        is_5prime=False, sig_sc=sig_sc_3prime_raw
                    )

                    # Default: forward strand assignment
                    sig_left = sig_sc_5prime
                    sig_right = sig_sc_3prime
                    seq_to_sig_map_left = seq_to_sig_map_sc_5prime
                    seq_to_sig_map_right = seq_to_sig_map_sc_3prime

                    # Reverse strand: flip and swap
                    if ref_reg_io.strand == "-":
                        (sig, sig_left, sig_right,
                         seq_to_sig_map, seq_to_sig_map_left, seq_to_sig_map_right) = \
                            process_reverse_strand(
                                sig, sig_sc_5prime, sig_sc_3prime,
                                seq_to_sig_map, seq_to_sig_map_sc_5prime,
                                seq_to_sig_map_sc_3prime
                            )

                    # Calculate reference-space coordinates
                    xCoord, xCoord_left, xCoord_right = calculate_xcoords(
                        seq_to_sig_map, seq_to_sig_map_left, seq_to_sig_map_right,
                        read_reg, ref_reg_io
                    )

                    # Concatenate full read coordinates and signal
                    xCoord_full = np.concatenate([
                        xCoord_left if len(xCoord_left) > 0 else np.array([]),
                        xCoord if len(xCoord) > 0 else np.array([]),
                        xCoord_right if len(xCoord_right) > 0 else np.array([])
                    ])
                    signal_final = np.concatenate([
                        sig_left if len(sig_left) > 0 else np.array([]),
                        sig if len(sig) > 0 else np.array([]),
                        sig_right if len(sig_right) > 0 else np.array([])
                    ])
                else:
                    sig = np.array([])
                    xCoord_full = np.array([])
                    signal_final = np.array([])
            else:
                # IORead not found; fall back to read_reg signal
                sig_full = read_reg.norm_signal
                sig = read_reg.norm_signal
                xCoord_full = read_reg.ref_sig_coords
                signal_final = sig

            start_base = int(xCoord_full[0])       # 0-based
            end_base = int(xCoord_full[-1]) + 1     # 1-based

            base_medians = {}
            if (xCoord_full is not None and signal_final is not None
                    and len(xCoord_full) > 0 and len(signal_final) > 0):
                for base in range(max(start_base + 1, position_columns[0]),
                                  min(end_base + 1, position_columns[-1]) + 1):
                    mask = (xCoord_full >= base - 1) & (xCoord_full < base)
                    if not np.any(mask):
                        continue
                    y_sub = signal_final[mask]
                    n = len(y_sub)

                    # Temporally central datapoint as representative value
                    mid = n // 2
                    if n % 2 == 0:
                        base_median = (y_sub[mid - 1] + y_sub[mid]) / 2
                    else:
                        base_median = y_sub[mid]

                    base_medians[base] = base_median

            row = [read_id, matched_start_base + 1, matched_end_base,
                   matched_end_base - matched_start_base]
            for pos in position_columns:
                row.append(base_medians.get(pos, np.nan))
            all_rows.append(row)

            percent = count / total_reads * 100
            sys.stdout.write(f"\r{count} reads processed ({percent:.1f}%)" + " " * 20)
            sys.stdout.flush()

            if len(all_rows) >= args.chunk_size:
                df_chunk = pd.DataFrame(all_rows, columns=columns)
                df_chunk.to_csv(args.output, mode='a', header=False, index=False)
                all_rows.clear()

    # Flush remaining rows
    if all_rows:
        df_chunk = pd.DataFrame(all_rows, columns=columns)
        df_chunk.to_csv(args.output, mode='a', header=False, index=False)
        all_rows.clear()

    print(f"\nTotal {count} reads processed.")
    print(f"\nOutput saved to: {args.output}")


if __name__ == "__main__":
    main()
