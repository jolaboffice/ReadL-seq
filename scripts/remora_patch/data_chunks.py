def make_sequence_coordinate_mapping2(cigar):
    """Maps an element in reference to every element in basecalls using
    alignment in `cigar`.

    Args:
        cigar (list): "cigartuples" representing alignment

    Returns:
        array shape (ref_len,). [x_0, x_1, ..., x_(ref_len)]
            such that read_seq[x_i] <> ref_seq[i]. Note that ref_len is derived
            from the cigar input.
    """
    ops, lens = map(np.array, zip(*cigar))
    assert ops.min() >= 0 and ops.max() <= 8, "invalid cigar op(s)"
    assert lens.min() >= 0, "cigar lengths may not be negative"

    is_match = MATCH_OPS[ops]
    match_counts = lens[is_match]
    offsets = np.array([match_counts, np.ones_like(match_counts)])

    # TODO remove knots around ambiguous indels (e.g. left justified HPs)
    # note this requires the ref and query sequences
    ref_knots = np.cumsum(np.where(REF_OPS[ops], lens, 0))
    ref_knots = np.concatenate(
        [[0], (ref_knots[is_match] - offsets).T.flatten(), [ref_knots[-1]]]
    )
    query_knots = np.cumsum(np.where(QUERY_OPS[ops], lens, 0))
    query_knots = np.concatenate(
        [[0], (query_knots[is_match] - offsets).T.flatten(), [query_knots[is_match][-1]]]
    )
    knots = np.interp(np.arange(ref_knots[-1] + 1), ref_knots, query_knots)

    return knots

