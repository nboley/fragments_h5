"""Regression test for float32 cumsum overflow in get_g_or_c_cumsum.

one_hot_encode_sequences returns float32.  Prior to commit 95c76f5 (2026-03-09),
get_g_or_c_cumsum accumulated in float32.  float32 has a 24-bit mantissa, so once
the running G/C count reaches 2**24 = 16,777,216, adding 1.0 rounds back to the
same value (ties-to-even).  Every subsequent cumsum[stop] - cumsum[start] is then
0, causing GC to silently become 0.0 for the remainder of the contig.

In production this zeroed GC past ~34-43 Mb on every chromosome (shipped in
releases v2.2.1 through v2.6.0 for ~4 months).

This test drives the accumulator past 2**24 with a synthetic all-G contig and
verifies that GC is computed correctly for fragments beyond the overflow threshold.
"""

import os
import tempfile

import numpy as np
import pysam
import pytest

from fragments_h5.fragment import get_g_or_c_cumsum

# 2**24 = 16_777_216.  In an all-GC sequence every base increments the
# accumulator, so we need at least 2**24 + 1 bases to cross the threshold.
# We add a small margin (100 bases) to place test fragments past the boundary.
OVERFLOW_THRESHOLD = 2**24  # 16_777_216
MARGIN = 100
CONTIG_LENGTH = OVERFLOW_THRESHOLD + MARGIN  # 16_777_316


@pytest.fixture(scope="module")
def overflow_fasta():
    """Create a temporary FASTA with a single all-G contig."""
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "overflow_test.fa")
        with open(fasta_path, "w") as f:
            f.write(">chrTest\n")
            seq = "G" * CONTIG_LENGTH
            # Write in 80-char lines per FASTA convention
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")

        # Index the FASTA for pysam
        pysam.faidx(fasta_path)
        yield fasta_path


@pytest.mark.slow
def test_gc_cumsum_no_float32_overflow(overflow_fasta):
    """GC must be correct for fragments past the float32 overflow threshold.

    With an all-G contig, every fragment should have GC = 1.0 regardless of
    position.  The buggy float32 accumulator would yield GC = 0.0 for any
    fragment starting at or past base 16,777,216.
    """
    g_or_c_cumsum, offset = get_g_or_c_cumsum(overflow_fasta, "chrTest")
    assert offset == 0

    # The cumsum array includes a leading pad base ('a'), so its length is
    # CONTIG_LENGTH + 1.
    assert len(g_or_c_cumsum) == CONTIG_LENGTH + 1

    # --- Verify dtype is float64 (cheap guard against reintroduction) ---
    assert g_or_c_cumsum.dtype == np.float64, (
        f"Expected float64 cumsum but got {g_or_c_cumsum.dtype}; "
        "float32 overflows at 2**24 G/C bases"
    )

    # --- Behavioral test: GC fraction for fragments past the threshold ---
    # Test several fragments straddling and past the overflow boundary.
    # For an all-G contig, cumsum[stop] - cumsum[start] should equal
    # (stop - start) exactly, giving GC = 1.0.
    test_fragments = [
        # (start, stop) — genomic coordinates
        (0, 200),                                          # beginning
        (OVERFLOW_THRESHOLD - 200, OVERFLOW_THRESHOLD),    # just before threshold
        (OVERFLOW_THRESHOLD, OVERFLOW_THRESHOLD + 50),     # right at threshold
        (OVERFLOW_THRESHOLD + 10, OVERFLOW_THRESHOLD + 60),  # past threshold
        (CONTIG_LENGTH - 50, CONTIG_LENGTH),               # end of contig
    ]

    for frag_start, frag_stop in test_fragments:
        frag_len = frag_stop - frag_start
        gc_count = g_or_c_cumsum[frag_stop] - g_or_c_cumsum[frag_start]
        gc_fraction = gc_count / frag_len

        assert gc_fraction == pytest.approx(1.0, abs=1e-9), (
            f"Fragment [{frag_start}, {frag_stop}) (len={frag_len}): "
            f"expected GC=1.0 but got {gc_fraction:.10f} "
            f"(gc_count={gc_count}, cumsum[stop]={g_or_c_cumsum[frag_stop]}, "
            f"cumsum[start]={g_or_c_cumsum[frag_start]})"
        )

    # --- Verify the accumulator value at the end is exact ---
    # For an all-G sequence of length L, cumsum[-1] should be exactly L
    # (the pad 'a' contributes 0 to G/C).  A float32 accumulator would
    # stick at 16777216.
    assert g_or_c_cumsum[-1] == pytest.approx(CONTIG_LENGTH, abs=0), (
        f"Final cumsum value should be {CONTIG_LENGTH} but got "
        f"{g_or_c_cumsum[-1]} (stuck at 2**24 = {2**24}?)"
    )
