"""Regression test for float32 cumsum overflow in get_g_or_c_cumsum.

one_hot_encode_sequences returns float32.  Prior to commit 95c76f5 (2026-03-09),
get_g_or_c_cumsum accumulated in float32.  float32 has a 24-bit mantissa, so once
the running G/C count reaches 2**24 = 16,777,216, adding 1.0 rounds back to the
same value (ties-to-even).  Every subsequent cumsum[stop] - cumsum[start] is then
0, causing GC to silently become 0.0 for the remainder of the contig.

In production this zeroed GC past ~34-43 Mb on every chromosome (shipped in
releases v2.2.1 through v2.6.0 for ~4 months).

Tests are split into fast (unmarked) and slow (@pytest.mark.slow):

- Fast tests use a tiny 80-bp non-uniform contig to verify dtype, offset
  arithmetic, and sub-region correctness.  These run in milliseconds and are
  never skipped by -m "not slow".
- The slow test drives the accumulator past 2**24 with a synthetic ~16.8M-base
  contig and verifies GC for fragments beyond the overflow threshold.
  It requires ~676 MB RSS.
"""

import os
import tempfile

import numpy as np
import pysam
import pytest

from fragments_h5.fragment import get_g_or_c_cumsum

# ---------------------------------------------------------------------------
# Constants for the overflow (slow) test
# ---------------------------------------------------------------------------

# 2**24 = 16_777_216.  In an all-GC sequence every base increments the
# accumulator, so we need at least 2**24 + 1 bases to cross the threshold.
OVERFLOW_THRESHOLD = 2**24  # 16_777_216
MARGIN = 100
CONTIG_LENGTH = OVERFLOW_THRESHOLD + MARGIN  # 16_777_316

# A short run of non-GC bases past the overflow boundary.  With an all-G
# contig, cumsum[i] == i, so a constant-offset indexing error is invisible:
# cumsum[stop-k] - cumsum[start-k] == stop-start for any k.  The A-patch
# breaks that invariance.
#
# Note the patch only defeats an offset error for fragments whose *endpoint*
# cuts through it.  A fragment that fully contains the patch still counts the
# same 20 non-G/C bases after a small shift, so it is NOT offset-sensitive.
# Both cases are asserted below, for different reasons.
A_PATCH_START = OVERFLOW_THRESHOLD + 50
A_PATCH_LENGTH = 20  # A-patch covers [A_PATCH_START, A_PATCH_START + 20)


# ---------------------------------------------------------------------------
# Helpers and fixtures
# ---------------------------------------------------------------------------

def _write_fasta(path, name, seq):
    """Write a single-contig FASTA and index it."""
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")
    pysam.faidx(path)


@pytest.fixture(scope="module")
def overflow_fasta():
    """A mostly-G contig crossing the float32 overflow boundary.

    All G except for A_PATCH_LENGTH A's starting at A_PATCH_START.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "overflow_test.fa")
        seq = (
            "G" * A_PATCH_START
            + "A" * A_PATCH_LENGTH
            + "G" * (CONTIG_LENGTH - A_PATCH_START - A_PATCH_LENGTH)
        )
        _write_fasta(fasta_path, "chrTest", seq)
        yield fasta_path


@pytest.fixture(scope="module")
def tiny_fasta():
    """A small non-uniform contig for fast tests: G*30 + A*20 + G*30 (80 bp)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "tiny_test.fa")
        _write_fasta(fasta_path, "chrTiny", "G" * 30 + "A" * 20 + "G" * 30)
        yield fasta_path


# ---------------------------------------------------------------------------
# Fast tests (unmarked — always run, even with -m "not slow")
# ---------------------------------------------------------------------------

def test_gc_cumsum_dtype_is_float64(tiny_fasta):
    """Cumsum dtype must be float64, not float32.

    This is a cheap early-warning against reintroduction of the float32
    accumulator, but it is NOT sufficient alone: a refactor to
    .cumsum().astype(numpy.float64) (accumulate in float32, then cast)
    would pass this check while producing fully saturated garbage past
    2**24 bases.  See test_gc_cumsum_no_float32_overflow for the true
    regression guard.
    """
    g_or_c_cumsum, _ = get_g_or_c_cumsum(tiny_fasta, "chrTiny")
    assert g_or_c_cumsum.dtype == np.float64, (
        f"Expected float64 cumsum but got {g_or_c_cumsum.dtype}; "
        "float32 overflows at 2**24 G/C bases"
    )


def test_gc_cumsum_correctness_nonuniform(tiny_fasta):
    """GC fractions must be correct on a non-uniform contig.

    An all-G contig has cumsum[i] == i, making cumsum[stop] - cumsum[start]
    equal stop - start regardless of any constant offset error.  This
    non-uniform contig (G*30 + A*20 + G*30) catches such indexing bugs.

    All comparisons use abs=0 (exact equality) because the cumsum contains
    integer-valued float64s and the divisions produce exact results for
    these specific fractions.
    """
    g_or_c_cumsum, offset = get_g_or_c_cumsum(tiny_fasta, "chrTiny")
    assert offset == 0
    assert len(g_or_c_cumsum) == 81  # 80 bases + 1 leading pad

    # Fragment in first G-run: [0, 10) -> GC = 10/10 = 1.0
    gc = (g_or_c_cumsum[10] - g_or_c_cumsum[0]) / 10
    assert gc == pytest.approx(1.0, abs=0), f"All-G fragment: expected 1.0, got {gc}"

    # Fragment in A-run: [30, 50) -> GC = 0/20 = 0.0
    gc = (g_or_c_cumsum[50] - g_or_c_cumsum[30]) / 20
    assert gc == pytest.approx(0.0, abs=0), f"All-A fragment: expected 0.0, got {gc}"

    # Spanning G->A boundary: [20, 40) -> 10G + 10A -> GC = 10/20 = 0.5
    gc = (g_or_c_cumsum[40] - g_or_c_cumsum[20]) / 20
    assert gc == pytest.approx(0.5, abs=0), f"Boundary fragment: expected 0.5, got {gc}"

    # Whole contig: [0, 80) -> 60G + 20A -> GC = 60/80 = 0.75
    gc = (g_or_c_cumsum[80] - g_or_c_cumsum[0]) / 80
    assert gc == pytest.approx(0.75, abs=0), f"Whole contig: expected 0.75, got {gc}"


def test_gc_cumsum_subregion_offset(tiny_fasta):
    """Sub-region queries must agree with whole-contig queries for the same fragment.

    get_g_or_c_cumsum(region_start, region_stop) returns a shorter cumsum and
    an offset; callers subtract the offset from genomic coordinates:
        gc = (cumsum[stop - offset] - cumsum[start - offset]) / length

    The original bug report wrongly blamed this code path for the GC corruption,
    so pinning its correctness has diagnostic value beyond regression safety.
    """
    full_cumsum, full_offset = get_g_or_c_cumsum(tiny_fasta, "chrTiny")
    assert full_offset == 0

    # Sub-region covering bases [20, 60) — spans both G/A boundaries
    region_start, region_stop = 20, 60
    sub_cumsum, sub_offset = get_g_or_c_cumsum(
        tiny_fasta, "chrTiny",
        region_start=region_start, region_stop=region_stop,
    )
    assert sub_offset == region_start

    test_fragments = [
        (20, 40),  # 10G + 10A -> GC = 0.5
        (30, 50),  # 20A -> GC = 0.0
        (40, 60),  # 10A + 10G -> GC = 0.5
        (25, 55),  # 5G + 20A + 5G -> GC = 10/30
    ]
    for frag_start, frag_stop in test_fragments:
        frag_len = frag_stop - frag_start
        gc_full = (full_cumsum[frag_stop] - full_cumsum[frag_start]) / frag_len
        gc_sub = (
            sub_cumsum[frag_stop - sub_offset]
            - sub_cumsum[frag_start - sub_offset]
        ) / frag_len
        # Integer-valued float64 cumsum — exact equality is appropriate
        assert gc_sub == pytest.approx(gc_full, abs=0), (
            f"Fragment [{frag_start}, {frag_stop}): "
            f"whole-contig GC={gc_full}, sub-region GC={gc_sub}"
        )


# ---------------------------------------------------------------------------
# Slow test (requires ~676 MB RSS)
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_gc_cumsum_no_float32_overflow(overflow_fasta):
    """GC must be correct for fragments past the float32 overflow threshold.

    Drives the accumulator past 2**24 with a synthetic ~16.8M-base contig.
    Requires ~676 MB RSS.

    The contig is all-G except for a short A-patch past the overflow boundary,
    making the test offset-aware: a constant-offset indexing error would
    produce the wrong GC for the A-patch fragment.

    Why this behavioral test is load-bearing and must not be deleted:
    A refactor to .cumsum().astype(numpy.float64) (accumulate in float32,
    cast afterward) leaves dtype == float64, passing the dtype check, while
    producing fully saturated garbage past 2**24 bases.  The dtype check is
    a cheap proxy that a one-line reordering defeats.  This test is the true
    regression guard against float32 overflow.
    """
    g_or_c_cumsum, offset = get_g_or_c_cumsum(overflow_fasta, "chrTest")
    assert offset == 0
    assert len(g_or_c_cumsum) == CONTIG_LENGTH + 1

    # --- Fragments in all-G regions ---
    all_g_fragments = [
        (0, 200),                                        # beginning
        (OVERFLOW_THRESHOLD - 200, OVERFLOW_THRESHOLD),  # just before threshold
        (OVERFLOW_THRESHOLD, OVERFLOW_THRESHOLD + 50),   # at threshold, before A-patch
        (A_PATCH_START + A_PATCH_LENGTH, CONTIG_LENGTH),  # after A-patch to end
    ]
    for frag_start, frag_stop in all_g_fragments:
        frag_len = frag_stop - frag_start
        gc_count = g_or_c_cumsum[frag_stop] - g_or_c_cumsum[frag_start]
        gc_fraction = gc_count / frag_len
        assert gc_fraction == pytest.approx(1.0, abs=0), (
            f"All-G fragment [{frag_start}, {frag_stop}): "
            f"expected GC=1.0 but got {gc_fraction:.10f}"
        )

    # --- Fragment fully containing the A-patch ---
    # [A_PATCH_START - 10, A_PATCH_START + A_PATCH_LENGTH + 10):
    #   10 G + 20 A + 10 G = 40 bases, 20 G/C -> expected GC = 0.5
    # This confirms non-G/C bases are still correctly excluded past 2**24 (a
    # saturated accumulator yields 0.0).  It is deliberately NOT an offset
    # check: shifting the window by <=10 keeps all 20 A's inside it, so the
    # count is unchanged.
    patch_start = A_PATCH_START - 10
    patch_stop = A_PATCH_START + A_PATCH_LENGTH + 10
    patch_len = patch_stop - patch_start  # 40
    gc_fraction = (g_or_c_cumsum[patch_stop] - g_or_c_cumsum[patch_start]) / patch_len
    assert gc_fraction == pytest.approx(0.5, abs=0), (
        f"Patch-containing fragment [{patch_start}, {patch_stop}) (len={patch_len}): "
        f"expected GC=0.5 but got {gc_fraction:.10f}"
    )

    # --- Offset-sensitive fragment whose start cuts through the A-patch ---
    # [A_PATCH_START + 10, A_PATCH_START + 50) = [.., CONTIG_LENGTH):
    #   10 A + 30 G = 40 bases, 30 G/C -> expected GC = 0.75
    # A one-base index shift gives 11 A + 29 G = 0.725, so this assertion
    # fails on any constant-offset error past the overflow boundary.
    cut_start = A_PATCH_START + 10
    cut_stop = A_PATCH_START + 50  # == CONTIG_LENGTH
    cut_len = cut_stop - cut_start  # 40
    gc_fraction = (g_or_c_cumsum[cut_stop] - g_or_c_cumsum[cut_start]) / cut_len
    assert gc_fraction == pytest.approx(0.75, abs=0), (
        f"Patch-cutting fragment [{cut_start}, {cut_stop}) (len={cut_len}): "
        f"expected GC=0.75 but got {gc_fraction:.10f} — offset check failed"
    )

    # --- Final accumulator value ---
    # Total G/C bases = CONTIG_LENGTH - A_PATCH_LENGTH (A's are not G/C)
    expected_total = CONTIG_LENGTH - A_PATCH_LENGTH
    assert g_or_c_cumsum[-1] == pytest.approx(expected_total, abs=0), (
        f"Final cumsum should be {expected_total} but got "
        f"{g_or_c_cumsum[-1]} (stuck at 2**24 = {2**24}?)"
    )
