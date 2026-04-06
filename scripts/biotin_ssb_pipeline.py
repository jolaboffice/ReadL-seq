#!/usr/bin/env python
"""
Integrated SSB detection pipeline for nanopore sequencing data.

Detects single-strand break (SSB) sites by computing BiotinScore for each
aligned read in a single pass through the data:

  POD5/BAM → signal extraction → k-mer baseline level lookup
           → delta-signal computation → BiotinScore → SSB classification

Instead of requiring control sequencing data as a baseline, this pipeline
uses the canonical k-mer signal level table from Remora (levels.txt) to
estimate expected signal levels at each position. The delta-signal
(predicted - observed) is then summed over a scoring window near the
read's 3' end to produce the BiotinScore.

The pipeline supports FASTA references with multiple contigs. Each contig
is processed independently, producing per-position CSV files (coverage and
biotin-positive hit counts) and a combined per-read XLSX file.

Outputs:
  <output_dir>/<contig>.csv   — per-position coverage and biotin hit counts
  <output_dir>/<output>.xlsx  — per-read BiotinScore and classification

Dependencies:
  pod5, pysam, biopython, remora, numpy, pandas, openpyxl

Usage:
  python biotin_ssb_pipeline.py \\
      --bam <BAM> --pod5 <POD5> --reference <FASTA> \\
      --level-table <LEVELS.txt> \\
      [--region chr:start-end] \\
      [--threshold 8.0] [--output OUTPUT_DIR]
"""
import argparse
import os
import random
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pod5
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from remora import io, refine_signal_map
from remora import data_chunks as DC
from remora.io import compute_base_space_sig_coords

# ======================== Constants ========================

# CIGAR operation code for soft-clipping
CIGAR_SOFT_CLIP = 4

# BiotinScore scoring window parameters
BIOTIN_SCORE_START_OFFSET = 3    # positions upstream of 3' end to start scoring
BIOTIN_SCORE_END_OFFSET = 100    # maximum positions downstream into soft-clip
BIOTIN_SCORE_THRESH = 8.0        # default threshold for biotin-positive classification

# Per-position counting: how far soft-clip extends coverage
TOTAL_READS_SOFTCLIP_LIMIT = 10

# Range around 3' end for identifying individual biotin positions
BIOTIN_POSITIONS_RANGE = 10

# K-mer model parameters (9-mer with center at index 6)
KMER_CENTER_IDX = 6
KMER_SIZE = 9


# ============================================================================
# Signal extraction utilities
# ============================================================================

def get_softclip_lengths(cigar_tuples: List[Tuple[int, int]]) -> Tuple[int, int]:
    """Extract 5' and 3' soft-clip lengths from CIGAR tuples.

    Returns:
        Tuple of (5-prime soft-clip length, 3-prime soft-clip length).
    """
    sc_5prime = 0
    sc_3prime = 0
    if cigar_tuples and len(cigar_tuples) > 0:
        if cigar_tuples[0][0] == CIGAR_SOFT_CLIP:
            sc_5prime = cigar_tuples[0][1]
        if cigar_tuples[-1][0] == CIGAR_SOFT_CLIP:
            sc_3prime = cigar_tuples[-1][1]
    return sc_5prime, sc_3prime


def refine_ref_to_signal(ref_to_signal: np.ndarray, io_read, sig_map_refiner,
                         read_id: str, seq: str,
                         sig: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Refine reference-to-signal mapping using Remora's SigMapRefiner.

    Applies iterative signal-to-reference alignment refinement to improve
    the accuracy of per-base signal segmentation.

    Args:
        ref_to_signal: Initial reference-to-signal index mapping.
        io_read: Remora IORead object containing raw signal data.
        sig_map_refiner: Remora SigMapRefiner instance.
        read_id: Read identifier (for error reporting).
        seq: Reference sequence for the aligned region.
        sig: Fallback normalized signal array.

    Returns:
        Tuple of (refined ref_to_signal mapping, normalized signal array).
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


def extract_matched_region_signal(ref_to_signal: np.ndarray, read_reg,
                                  ref_reg_io) -> Tuple[np.ndarray, np.ndarray]:
    """Extract signal mapping for the matched (aligned) region.

    Trims the full reference-to-signal mapping to cover only the aligned
    portion of the read, accounting for strand orientation.

    Returns:
        Tuple of (trimmed ref_to_signal, zero-based seq_to_sig_map).
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


def process_softclip_xcoord(query_to_signal: np.ndarray, sc_len: int, io_read,
                            is_5prime: bool,
                            sig_sc: Optional[np.ndarray] = None
                            ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute signal mapping for a soft-clip region.

    Extracts the query-to-signal indices corresponding to a soft-clipped
    segment and produces the normalized signal for that region.

    Args:
        query_to_signal: Full query-to-signal index mapping.
        sc_len: Length of the soft-clip in bases.
        io_read: Remora IORead object.
        is_5prime: True for 5' soft-clip, False for 3'.
        sig_sc: Pre-extracted soft-clip signal (optional fallback).

    Returns:
        Tuple of (ref_to_signal, seq_to_sig_map, signal) for the soft-clip.
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


def reverse_strand_arrays(sig: np.ndarray, sig_sc_5prime: np.ndarray,
                          sig_sc_3prime: np.ndarray,
                          seq_to_sig_map: np.ndarray,
                          seq_to_sig_map_sc_5prime: np.ndarray,
                          seq_to_sig_map_sc_3prime: np.ndarray
                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                     np.ndarray, np.ndarray, np.ndarray]:
    """Reverse signal and coordinate arrays for reverse-strand reads.

    Nanopore reads are sequenced 5' to 3', so reverse-strand reads must have
    their signal and coordinate arrays flipped to align with the forward
    reference orientation. The 5' and 3' soft-clip assignments are also
    swapped.

    Returns:
        Tuple of (sig, sig_left, sig_right, seq_to_sig_map,
        seq_to_sig_map_left, seq_to_sig_map_right) with reversed orientation.
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


def calculate_xcoords(seq_to_sig_map: np.ndarray, seq_to_sig_map_left: np.ndarray,
                      seq_to_sig_map_right: np.ndarray,
                      read_reg, ref_reg_io
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate reference-space x-coordinates for signal visualization.

    Maps signal indices back to reference genome positions for the matched
    region and both soft-clip flanks.

    Returns:
        Tuple of (matched coords, left flank coords, right flank coords).
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
# K-mer baseline level functions
# ============================================================================

def load_kmer_levels(level_table_path: str) -> Dict[str, float]:
    """Load the canonical k-mer signal level table (Remora levels.txt).

    This table serves as a control-sequencing-free baseline for
    delta-signal estimation. Each entry maps a 9-mer to its expected
    (canonical) signal level, enabling estimation of signal deviations
    without requiring a separate control sequencing dataset.

    The file format is a two-column TSV: k-mer <TAB> expected_level.

    Returns:
        Dict mapping each 9-mer string to its canonical signal level.
    """
    kmer_levels = {}
    with open(level_table_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                kmer_levels[parts[0].upper()] = float(parts[1])
    return kmer_levels


def get_predicted_level(pos_0based: int, sequence: str,
                        kmer_levels: Dict[str, float]) -> float:
    """Get the canonical baseline signal level for a genomic position.

    Extracts the 9-mer centered on the given position from the reference
    sequence and looks up its canonical signal level from the Remora
    k-mer level table.

    Args:
        pos_0based: 0-based position in reference coordinates.
        sequence: Reference nucleotide sequence.
        kmer_levels: Canonical k-mer signal level table (from levels.txt).

    Returns:
        Canonical signal level, or 0.0 if the k-mer is not found or the
        position is too close to the sequence boundary.
    """
    start = pos_0based - KMER_CENTER_IDX
    end = start + KMER_SIZE
    if start < 0 or end > len(sequence):
        return 0.0
    kmer = sequence[start:end].upper()
    return kmer_levels.get(kmer, 0.0)


# ============================================================================
# Region parsing
# ============================================================================

def parse_regions(region_strs: Optional[List[str]]
                  ) -> Optional[List[Tuple[str, Optional[int], Optional[int]]]]:
    """Parse one or more genomic region strings.

    Supported formats per entry:
        'chr1:1000-2000'  - specific range (1-based input, converted to 0-based)
        'chr1'            - entire contig

    Args:
        region_strs: List of region strings, or None.

    Returns:
        List of (contig, start, end) tuples, or None if no regions specified.
        Coordinates are 0-based half-open [start, end).
    """
    if region_strs is None:
        return None
    regions = []
    for r in region_strs:
        if ":" in r:
            ctg, coords = r.split(":", 1)
            start, end = coords.split("-")
            regions.append((ctg, int(start) - 1, int(end)))
        else:
            regions.append((r, None, None))
    return regions


def validate_regions(regions: List[Tuple[str, Optional[int], Optional[int]]],
                     references: Dict[str, str]) -> None:
    """Validate parsed regions against the loaded reference sequences.

    Checks that each contig name exists in the reference and that
    coordinate ranges are within bounds. Exits with an error message
    if validation fails.
    """
    ref_info = "\n".join(
        f"    {name} ({len(seq)} bp)" for name, seq in references.items()
    )
    for ctg, start, end in regions:
        if ctg not in references:
            print(f"\n[Error] Contig '{ctg}' not found in reference FASTA.")
            print(f"  Available contigs:\n{ref_info}")
            sys.exit(1)
        if start is not None:
            ctg_len = len(references[ctg])
            if start < 0 or end > ctg_len or start >= end:
                print(f"\n[Error] Invalid coordinates for '{ctg}:{start+1}-{end}'.")
                print(f"  '{ctg}' length is {ctg_len} bp (valid range: 1-{ctg_len}).")
                sys.exit(1)


# ============================================================================
# Per-read signal extraction
# ============================================================================

def extract_read_signal(read_reg, io_read_dict: Dict, read_cigar_dict: Dict,
                        sequence: str, sig_map_refiner
                        ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray],
                                   int, int]:
    """Extract per-base signal coordinates and normalized signal for one read.

    Performs reference-to-signal mapping, Remora refinement, soft-clip
    handling, and reverse-strand flipping to produce a complete per-base
    signal trace covering the full read extent (including soft-clips).

    Args:
        read_reg: Remora read-reference region object.
        io_read_dict: Dict mapping read_id to Remora IORead objects.
        read_cigar_dict: Dict mapping read_id to CIGAR tuples.
        sequence: Reference nucleotide sequence for the current contig.
        sig_map_refiner: Remora SigMapRefiner instance.

    Returns:
        Tuple of (xCoord_full, signal_final, sc_5prime_len, sc_3prime_len).
        Returns (None, None, 0, 0) if signal extraction fails.
    """
    read_id = read_reg.read_id

    # Get signal and mappings from IORead
    if read_id in io_read_dict:
        io_read = io_read_dict[read_id]
        sig_full = io_read.norm_signal
        query_to_signal = io_read.query_to_signal
        ref_reg_io = io_read.ref_reg
    else:
        sig_full = read_reg.norm_signal
        query_to_signal = None
        ref_reg_io = None

    # Get soft-clip lengths from CIGAR
    sc_5prime_len = 0
    sc_3prime_len = 0
    cigar = read_cigar_dict.get(read_id, None)
    if cigar:
        sc_5prime_len, sc_3prime_len = get_softclip_lengths(cigar)

    # Fast path: no detailed mapping available
    if query_to_signal is None or cigar is None or ref_reg_io is None:
        return read_reg.ref_sig_coords, read_reg.norm_signal, sc_5prime_len, sc_3prime_len

    # Build reference-to-signal mapping via CIGAR
    knots = DC.make_sequence_coordinate_mapping2(cigar)
    ref_to_signal = DC.map_ref_to_signal(
        query_to_signal=query_to_signal,
        ref_to_query_knots=knots
    )

    # Prepare matched sequence for refinement
    match_start = ref_reg_io.start
    match_end = ref_reg_io.end
    seq_slice = sequence[match_start:match_end]
    if ref_reg_io.strand == "-":
        matched_seq = str(Seq(seq_slice).reverse_complement())
    else:
        matched_seq = seq_slice

    sig_fallback = (sig_full[ref_to_signal[0]:ref_to_signal[-1]]
                    if len(ref_to_signal) > 0 else np.array([]))

    # Refine signal-to-reference alignment
    ref_to_signal, sig = refine_ref_to_signal(
        ref_to_signal=ref_to_signal,
        io_read=io_read,
        sig_map_refiner=sig_map_refiner,
        read_id=read_id,
        seq=matched_seq,
        sig=sig_fallback
    )

    # Trim to matched region
    ref_to_signal, seq_to_sig_map = extract_matched_region_signal(
        ref_to_signal, read_reg, ref_reg_io
    )

    if len(sig) == 0 or len(ref_to_signal) == 0:
        return None, None, sc_5prime_len, sc_3prime_len

    sig_start_idx = ref_to_signal[0]
    sig_end_idx = ref_to_signal[-1]

    # Extract soft-clip signals
    sig_sc_5prime_raw = sig_full[0:sig_start_idx]
    sig_sc_3prime_raw = sig_full[sig_end_idx:]

    _, seq_to_sig_map_sc_5prime, sig_sc_5prime = process_softclip_xcoord(
        query_to_signal, sc_5prime_len, io_read,
        is_5prime=True, sig_sc=sig_sc_5prime_raw
    )
    _, seq_to_sig_map_sc_3prime, sig_sc_3prime = process_softclip_xcoord(
        query_to_signal, sc_3prime_len, io_read,
        is_5prime=False, sig_sc=sig_sc_3prime_raw
    )

    # Assign flanks (forward strand default)
    sig_left = sig_sc_5prime
    sig_right = sig_sc_3prime
    seq_to_sig_map_left = seq_to_sig_map_sc_5prime
    seq_to_sig_map_right = seq_to_sig_map_sc_3prime

    # Reverse strand: flip signals and swap flanks
    if ref_reg_io.strand == "-":
        (sig, sig_left, sig_right,
         seq_to_sig_map, seq_to_sig_map_left, seq_to_sig_map_right) = \
            reverse_strand_arrays(
                sig, sig_sc_5prime, sig_sc_3prime,
                seq_to_sig_map, seq_to_sig_map_sc_5prime,
                seq_to_sig_map_sc_3prime
            )

    # Calculate reference-space coordinates
    xCoord, xCoord_left, xCoord_right = calculate_xcoords(
        seq_to_sig_map, seq_to_sig_map_left, seq_to_sig_map_right,
        read_reg, ref_reg_io
    )

    # Concatenate full read extent
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

    return xCoord_full, signal_final, sc_5prime_len, sc_3prime_len


# ============================================================================
# BiotinScore computation
# ============================================================================

def compute_base_medians(xCoord_full: np.ndarray, signal_final: np.ndarray,
                         start: int, end: int) -> Dict[int, float]:
    """Compute per-base median signal values over a position range.

    For each 1-based position in [start, end], selects signal data points
    whose reference coordinate falls within [pos-1, pos) and computes
    the median.

    Args:
        xCoord_full: Reference-space coordinates for each signal data point.
        signal_final: Normalized signal values corresponding to xCoord_full.
        start: First 1-based position (inclusive).
        end: Last 1-based position (inclusive).

    Returns:
        Dict mapping 1-based position to median signal value.
    """
    medians = {}
    for base_pos in range(start, end + 1):
        mask = (xCoord_full >= base_pos - 1) & (xCoord_full < base_pos)
        if not np.any(mask):
            continue
        y_sub = signal_final[mask]
        n = len(y_sub)
        mid = n // 2
        if n % 2 == 0:
            median_val = (y_sub[mid - 1] + y_sub[mid]) / 2
        else:
            median_val = y_sub[mid]
        medians[base_pos] = float(median_val)
    return medians


def score_read(read_end: int, sc_3p_len: int, ref_length: int,
               xCoord_full: np.ndarray, signal_final: np.ndarray,
               sequence: str, kmer_levels: Dict[str, float],
               threshold: float) -> Tuple[float, bool, List[int]]:
    """Compute BiotinScore for a single read and classify it.

    The scoring window spans from (read_end - START_OFFSET) to
    (read_end + effective_window), where effective_window is bounded
    by both END_OFFSET and the actual 3' soft-clip length.

    Args:
        read_end: 1-based inclusive end position of the aligned region.
        sc_3p_len: 3-prime soft-clip length in bases.
        ref_length: Total reference contig length.
        xCoord_full: Reference-space coordinates for the read's signal.
        signal_final: Normalized signal values for the read.
        sequence: Reference nucleotide sequence for the contig.
        kmer_levels: Canonical k-mer signal level table (from levels.txt).
        threshold: BiotinScore threshold for positive classification.

    Returns:
        Tuple of (score, is_biotin, biotin_positions).
    """
    # Define scoring window (1-based positions)
    anchor_pos = read_end
    effective_window = min(BIOTIN_SCORE_END_OFFSET, sc_3p_len + 1)
    win_start = max(anchor_pos - BIOTIN_SCORE_START_OFFSET, 1)
    win_end = min(anchor_pos + effective_window - 1, ref_length)

    # Compute base medians for scoring + biotin position windows

    median_start = max(min(win_start, read_end - BIOTIN_POSITIONS_RANGE), 1)
    median_end = min(max(win_end, read_end + BIOTIN_POSITIONS_RANGE), ref_length)
    base_medians = compute_base_medians(xCoord_full, signal_final,
                                        median_start, median_end)

    # Sum delta-signal (canonical baseline - observed) across the scoring window
    score = 0.0
    for pos in range(win_start, win_end + 1):
        observed = base_medians.get(pos)
        if observed is None:
            continue
        predicted = get_predicted_level(pos - 1, sequence, kmer_levels)
        score += predicted - observed

    is_biotin = (score >= threshold)

    # Identify individual biotin positions (positive delta near 3' end)
    biotin_positions = []
    if is_biotin:
        bp_start = max(read_end - BIOTIN_POSITIONS_RANGE, 1)
        bp_end = min(read_end + BIOTIN_POSITIONS_RANGE, ref_length)
        for pos in range(bp_start, bp_end + 1):
            observed = base_medians.get(pos)
            if observed is None:
                continue
            predicted = get_predicted_level(pos - 1, sequence, kmer_levels)
            if predicted - observed > 0:
                biotin_positions.append(pos)

    return score, is_biotin, biotin_positions


# ============================================================================
# Per-position output
# ============================================================================

def write_per_position_csv(total_counts: Dict[int, int],
                           biotin_counts: Dict[int, int],
                           sequence: str, ref_length: int,
                           csv_path: str) -> None:
    """Write per-position coverage and biotin hit counts to CSV.

    Each row contains a 1-based position, the reference base, total read
    coverage, biotin-positive hit count, and the biotin ratio.

    Args:
        total_counts: Dict mapping 1-based position to total read count.
        biotin_counts: Dict mapping 1-based position to biotin hit count.
        sequence: Reference nucleotide sequence.
        ref_length: Total reference contig length.
        csv_path: Output CSV file path.
    """
    all_positions = sorted(set(total_counts.keys()) | set(biotin_counts.keys()))
    if all_positions:
        pos_arr = np.array(all_positions)
        base_arr = np.array([
            sequence[p - 1] if 0 < p <= ref_length else '-'
            for p in all_positions
        ])
        total_arr = np.array([total_counts[p] for p in all_positions])
        biotin_arr = np.array([biotin_counts[p] for p in all_positions])
        ratio_arr = np.divide(
            biotin_arr.astype(float), total_arr.astype(float),
            out=np.zeros(len(all_positions), dtype=float),
            where=(total_arr > 0)
        )
        df = pd.DataFrame({
            'position': pos_arr,
            'base': base_arr,
            'total_reads': total_arr,
            'biotin_hits': biotin_arr,
            'biotin_ratio': ratio_arr,
        })
    else:
        df = pd.DataFrame(
            columns=['position', 'base', 'total_reads', 'biotin_hits', 'biotin_ratio']
        )
    df.to_csv(csv_path, index=False)


# ============================================================================
# Main pipeline
# ============================================================================

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Integrated SSB detection pipeline: "
                    "signal extraction, predicted level lookup, delta computation, "
                    "and BiotinScore classification."
    )
    parser.add_argument("--bam", required=True, help="Sorted BAM file path")
    parser.add_argument("--pod5", required=True, help="POD5 file path")
    parser.add_argument("--reference", required=True, help="Reference FASTA file path")
    parser.add_argument("--level-table", required=True,
                        help="Canonical k-mer signal level table from Remora "
                             "(levels.txt, TSV: kmer<TAB>level)")
    parser.add_argument("--region", nargs="+", default=None,
                        help="Analysis region(s). Formats: 'chr1:1000-2000' or 'chr1'. "
                             "Multiple regions can be specified.")
    parser.add_argument("--threshold", type=float, default=BIOTIN_SCORE_THRESH,
                        help=f"BiotinScore threshold (default: {BIOTIN_SCORE_THRESH})")
    parser.add_argument("--max-reads", type=int, default=None,
                        help="Maximum number of reads to process")
    parser.add_argument("--output", default=None,
                        help="Output directory name (also used as XLSX filename prefix)")
    return parser.parse_args()


def process_contig(ctg_name: str, sequence: str, region_map: Dict,
                   bam_fh, pod5_fh, sig_map_refiner,
                   kmer_levels: Dict[str, float], threshold: float,
                   max_reads: Optional[int], out_dir: str
                   ) -> Tuple[List[dict], int, int]:
    """Process all reads for a single contig.

    Extracts signals, computes BiotinScore for each read, writes
    per-position CSV, and returns per-read results.

    Args:
        ctg_name: Contig name.
        sequence: Reference nucleotide sequence for the contig.
        region_map: Dict of contig to (start, end) bounds or None.
        bam_fh: Open pysam AlignmentFile handle.
        pod5_fh: Open pod5 Reader handle.
        sig_map_refiner: Remora SigMapRefiner instance.
        kmer_levels: Canonical k-mer signal level table (from levels.txt).
        threshold: BiotinScore threshold.
        max_reads: Maximum reads to process (None for all).
        out_dir: Output directory path.

    Returns:
        Tuple of (per-read results list, total read count, biotin count).
    """
    ref_length = len(sequence)
    ref_reg = io.RefRegion(ctg=ctg_name, strand="+", start=0, end=ref_length)

    # Fetch reads
    print(f"[Extract] Fetching reads for contig '{ctg_name}' ({ref_length} bp)...")
    sample_bam_reads = io.get_reg_bam_reads(ref_reg, bam_fh)
    if len(sample_bam_reads) == 0:
        print(f"  No reads for contig '{ctg_name}', skipping.")
        return [], 0, 0
    print(f"  Found {len(sample_bam_reads)} reads")

    if max_reads is not None and len(sample_bam_reads) > max_reads:
        sample_bam_reads = random.sample(sample_bam_reads, max_reads)
        print(f"  Subsampled to {max_reads} reads")

    # Create IORead objects
    io_reads, skipped_count = io.get_io_reads(
        sample_bam_reads, pod5_fh, reverse_signal=False, missing_ok=True
    )
    print(f"  Skipped {skipped_count} reads (move table/signal mismatch or missing POD5)")

    # Extract reference regions and build lookup tables
    read_ref_regs = [io_read.extract_ref_reg(ref_reg) for io_read in io_reads]
    print(f"[Extract] Building lookup tables for '{ctg_name}'...")
    io_read_dict = {io_read.read_id: io_read for io_read in io_reads}
    read_cigar_dict = {}
    for bam_read in sample_bam_reads:
        if bam_read.is_unmapped:
            continue
        read_cigar_dict[bam_read.query_name] = bam_read.cigartuples
    print(f"  IORead: {len(io_read_dict)}, CIGAR: {len(read_cigar_dict)}")

    # Process each read
    print(f"[Score] Computing BiotinScore for contig '{ctg_name}'...")
    total_reads = len(read_ref_regs)
    results = []
    total_counts = defaultdict(int)
    biotin_counts = defaultdict(int)

    for count, read_reg in enumerate(read_ref_regs, 1):
        # Signal extraction
        xCoord_full, signal_final, sc_5p_len, sc_3p_len = extract_read_signal(
            read_reg, io_read_dict, read_cigar_dict, sequence, sig_map_refiner
        )
        if xCoord_full is None or len(xCoord_full) == 0:
            continue

        # Read alignment info (1-based positions)
        read_start = read_reg.ref_reg.start + 1
        read_end = read_reg.ref_reg.end          # 1-based inclusive
        read_len = read_end - read_start + 1

        # BiotinScore computation
        score, is_biotin, biotin_positions = score_read(
            read_end, sc_3p_len, ref_length,
            xCoord_full, signal_final,
            sequence, kmer_levels, threshold
        )

        # Store per-read result
        results.append({
            'read_id': read_reg.read_id,
            'contig': ctg_name,
            'start': read_start,
            'end': read_end,
            'length': read_len,
            'softclip_3prime_length': sc_3p_len,
            'end_plus_softclip': read_end + sc_3p_len,
            'biotin_score': round(score, 4),
            'is_biotin': is_biotin,
            'biotin_count': len(biotin_positions),
            'biotin_positions': ', '.join(map(str, biotin_positions)) if biotin_positions else '',
        })

        # Per-position counting
        real_end = min(read_end + min(sc_3p_len, TOTAL_READS_SOFTCLIP_LIMIT), ref_length)
        for pos in range(read_start, real_end + 1):
            total_counts[pos] += 1
        if is_biotin:
            biotin_range_start = read_end - BIOTIN_POSITIONS_RANGE
            for pos in range(biotin_range_start, real_end + 1):
                biotin_counts[pos] += 1

        # Progress
        if count % 100 == 0 or count == total_reads:
            pct = count / total_reads * 100
            sys.stdout.write(f"\r  {count}/{total_reads} reads ({pct:.1f}%)")
            sys.stdout.flush()

    print()  # newline after progress

    # Write per-position CSV
    n_biotin = sum(1 for r in results if r['is_biotin'])
    csv_path = os.path.join(out_dir, f"{ctg_name}.csv")
    write_per_position_csv(total_counts, biotin_counts, sequence, ref_length, csv_path)

    print(f"  [{ctg_name}] {len(results)} reads, {n_biotin} biotin-positive")
    print(f"    per-position -> {csv_path}")

    return results, len(results), n_biotin


def main() -> None:
    """Run the integrated SSB detection pipeline.

    Pipeline steps:
      1. Load reference FASTA and canonical k-mer signal level table.
      2. Initialize Remora SigMapRefiner for signal alignment refinement.
      3. For each contig, extract per-read signals from POD5/BAM.
      4. Compute per-base signal medians and delta-signal (baseline - observed).
      5. Sum delta-signals in the scoring window to produce BiotinScore.
      6. Classify reads as biotin-positive if BiotinScore >= threshold.
      7. Output per-position CSV and per-read XLSX files.
    """
    args = parse_args()

    # ---- Load reference ----
    print("[Init] Loading reference...")
    references = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(args.reference, "fasta")
    }
    for ctg_name, seq in references.items():
        print(f"  Reference: {ctg_name} ({len(seq)} bp)")
    print(f"  Total contigs: {len(references)}")

    # ---- Load canonical k-mer signal level table (baseline) ----
    print("[Init] Loading canonical k-mer signal level table...")
    kmer_levels = load_kmer_levels(args.level_table)
    print(f"  Loaded {len(kmer_levels)} k-mer entries")

    # ---- Initialize SigMapRefiner ----
    print("[Init] Initializing SigMapRefiner...")
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=args.level_table,
        scale_iters=0,
        do_fix_guage=True,
    )
    print("  Done.")

    # ---- Parse and validate regions ----
    regions = parse_regions(args.region)
    if regions is not None:
        validate_regions(regions, references)
        for ctg, start, end in regions:
            if start is not None:
                print(f"[Init] Region: {ctg}:{start + 1}-{end}")
            else:
                print(f"[Init] Region: {ctg} (entire contig)")

    # ---- Open BAM/POD5 ----
    pod5_fh = pod5.Reader(args.pod5)
    bam_fh = pysam.AlignmentFile(args.bam)

    # ---- Determine contigs to process ----
    if regions is not None:
        region_map = {}
        for ctg, start, end in regions:
            region_map[ctg] = (start, end) if start is not None else None
        contigs_to_process = {ctg: references[ctg] for ctg in region_map}
    else:
        region_map = {}
        contigs_to_process = references

    # ---- Output directory ----
    if args.output:
        out_dir_name = args.output
    else:
        thresh_val = (int(args.threshold)
                      if args.threshold == int(args.threshold)
                      else args.threshold)
        out_dir_name = f"BiotinScore_{thresh_val}"
    out_dir = os.path.join(os.path.dirname(args.bam) or ".", out_dir_name)
    os.makedirs(out_dir, exist_ok=True)
    print(f"[Init] Output directory: {out_dir}")

    # ---- Process each contig ----
    all_results = []
    total_reads_all = 0
    total_biotin_all = 0

    for ctg_name, sequence in contigs_to_process.items():
        results, n_reads, n_biotin = process_contig(
            ctg_name, sequence, region_map,
            bam_fh, pod5_fh, sig_map_refiner,
            kmer_levels, args.threshold, args.max_reads, out_dir
        )
        all_results.extend(results)
        total_reads_all += n_reads
        total_biotin_all += n_biotin

    # ---- Write per-read XLSX (all contigs combined) ----
    xlsx_path = os.path.join(out_dir, f"{out_dir_name}.xlsx")
    print("[Output] Writing per-read results...")
    df_per_read = pd.DataFrame(all_results)
    if len(df_per_read) == 0:
        df_per_read = pd.DataFrame(columns=[
            'read_id', 'contig', 'start', 'end', 'length',
            'softclip_3prime_length', 'end_plus_softclip',
            'biotin_score', 'is_biotin', 'biotin_count', 'biotin_positions'
        ])
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as w:
        df_per_read.to_excel(w, index=False, sheet_name="per_read")

    # ---- Summary ----
    if total_reads_all > 0:
        print(f"[Done] {total_reads_all} reads processed, {total_biotin_all} biotin-positive "
              f"({total_biotin_all / total_reads_all * 100:.1f}%)")
    else:
        print("[Done] No reads.")
    print(f"  per-read -> {xlsx_path}")


if __name__ == "__main__":
    main()
