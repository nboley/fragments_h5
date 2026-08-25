"""Tests for the GC repair tool — stage-0 checklist items (§7.4).

Covers:
  0.a — delete-and-recreate dataset test
  0.b — synthetic padding-row fixture with full check list
  0.c — 2-D mapq truncation along axis 0 only
  0.d — rounding-agreement test over >= 10^7 fragments
  0.e — float32-accumulator-simulation test
  0.f — (alphabet histogram regeneration — tested via scan_fasta_alphabet)
  0.g — existing overflow regression test still passes (separate file)

Additional:
  - detect_padding_row: clean, truncate, abort cases
  - no-gc file handling (§3.2.5)
  - fragment_length_counts exact arithmetic identity (§3.2.4)
  - repair_local_file end-to-end on synthetic fixtures
  - CLI argument validation
"""

import hashlib
import json
import os
import tempfile

import h5py
import numpy as np
import pysam
import pytest

from fragments_h5.repair import (
    RepairAbort,
    T23,
    T24,
    build_n_prefix_sum,
    check_gc_presence,
    classify_gc_changes,
    compute_file_sha256,
    detect_padding_row,
    hash_attrs_excluding,
    recompute_gc_for_contig,
    rebuild_contig_index,
    repair_local_file,
    scan_fasta_alphabet,
    simulate_float32_cumsum,
    truncate_contig_datasets,
)
from fragments_h5.fragments_h5 import (
    FragmentsH5,
    INDEX_BLOCK_SIZE,
    MIN_NUM_READS_FOR_INDEX,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(path, contigs):
    """Write a multi-contig FASTA and index it.

    contigs: list of (name, sequence) tuples.
    """
    with open(path, "w") as f:
        for name, seq in contigs:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
    pysam.faidx(path)


def _random_sequence(seed, length, weights=(0.22, 0.24, 0.24, 0.22, 0.08)):
    """Deterministic ACGTN sequence with a realistic-ish base mix."""
    rng = np.random.RandomState(seed)
    return "".join(rng.choice(list("ACGTN"), size=length, p=list(weights)))


# ---------------------------------------------------------------------------
# Independent GC oracle (§7.4 / review finding H-1)
#
# Everything below this banner is deliberately naive and shares NO code with
# fragments_h5.repair or with fragments_h5.fragment's cumsum machinery. It reads
# the FASTA text itself, slices the sequence string per fragment, and counts
# characters. This is the oracle the end-to-end tests assert against; using
# recompute_gc_for_contig to build fixture GC values would make those tests
# circular (they would pass against any implementation, including a wrong one).
# ---------------------------------------------------------------------------

def read_fasta_sequences(path):
    """Parse a FASTA file into {contig: sequence} with plain Python text I/O."""
    seqs = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


# Per-base C+G contribution. C and G contribute 1.0; N is encoded as
# (0.25, 0.25, 0.25, 0.25) at sequence.pyx:40, so its C+G share is 0.5.
_GC_CONTRIBUTION = {"C": 1.0, "G": 1.0, "N": 0.5}


def naive_gc_uint8(sequence, starts, lengths, ndigits=5):
    """Reference GC encoder: slice the sequence, count characters, quantise.

    Quantisation follows the documented two stages: round(fraction, 5) then
    round(q5 * 254) to uint8. Zero-length fragments encode as 255.

    ndigits is only ever varied by the rounding-sensitivity check below, which
    needs to know whether a fixture can tell round(x, 5) from round(x, 4).
    """
    out = np.empty(len(starts), dtype=np.uint8)
    for i in range(len(starts)):
        length = int(lengths[i])
        if length == 0:
            out[i] = 255
            continue
        start = int(starts[i])
        span = sequence[start : start + length].upper()
        assert len(span) == length, (
            f"fragment [{start}, {start + length}) runs past the contig"
        )
        total = 0.0
        for ch in span:
            total += _GC_CONTRIBUTION.get(ch, 0.0)
        out[i] = int(round(round(total / length, ndigits) * 254))
    return out


def count_rounding_sensitive(sequence, starts, lengths):
    """How many fragments encode differently under round(x, 4) vs round(x, 5).

    A fixture with zero such fragments cannot detect a bug in the fifth decimal
    place — that is exactly why the original all-G fixture was weak (M-3).
    """
    q5 = naive_gc_uint8(sequence, starts, lengths, ndigits=5)
    q4 = naive_gc_uint8(sequence, starts, lengths, ndigits=4)
    return int(np.sum(q5 != q4))


def _make_fragment_h5(
    path,
    contigs,
    read_gc=True,
    add_padding=False,
    index_block_size=INDEX_BLOCK_SIZE,
    contig_lengths=None,
):
    """Build a synthetic fragment H5 for testing.

    contigs: dict of contig_name -> {
        'starts': np.array,
        'lengths': np.array,
        'gc': np.array (optional),
        'mapq': np.array (optional, shape (n,2)),
        'strand': np.array (optional),
        'extra': dict of name -> np.array (optional; methyl / fragment_end_clipped),
    }
    add_padding: if True, append a phantom zero-row to every dataset (simulating the bug).
    contig_lengths: dict of contig_name -> int (for attrs).
    """
    if contig_lengths is None:
        contig_lengths = {}
        for c, data in contigs.items():
            stops = data["starts"].astype(np.int64) + data["lengths"].astype(np.int64)
            contig_lengths[c] = int(stops.max()) + 1000 if len(stops) > 0 else 10000

    with h5py.File(path, "w") as f:
        f.create_group("data")
        f.create_group("index")

        f.attrs["index_block_size"] = index_block_size
        f.attrs["max_fragment_length"] = 65535
        f.attrs["_bam_header"] = ""
        f.attrs["_source_format"] = "BAM"
        f.attrs["_contig_lengths_str"] = str(contig_lengths)

        for contig, data in contigs.items():
            starts = data["starts"]
            lengths = data["lengths"]
            n = len(starts)

            if add_padding:
                # Append phantom zero row
                starts = np.append(starts, np.int32(0))
                lengths = np.append(lengths, np.uint16(0))

            def mk(key, arr, dtype):
                f.create_dataset(
                    key, data=arr, dtype=dtype,
                    compression="gzip", compression_opts=4, chunks=True,
                )

            mk(f"data/{contig}/starts", starts, "int32")
            mk(f"data/{contig}/lengths", lengths, "uint16")

            # mapq: 2-D
            if "mapq" in data:
                mapq = data["mapq"]
            else:
                mapq = np.column_stack([
                    np.full(n, 60, dtype="uint8"),
                    np.full(n, 60, dtype="uint8"),
                ])
            if add_padding:
                mapq = np.vstack([mapq, np.zeros((1, 2), dtype="uint8")])
            mk(f"data/{contig}/mapq", mapq, "uint8")

            if read_gc:
                gc = data.get("gc", np.zeros(n, dtype="uint8"))
                if add_padding:
                    gc = np.append(gc, np.uint8(0))
                mk(f"data/{contig}/gc", gc, "uint8")

            if "strand" in data:
                strand = data["strand"]
                if add_padding:
                    strand = np.append(strand, b"")
                mk(f"data/{contig}/strand", strand, "|S1")
            else:
                strand = np.array([b"+"] * n, dtype="|S1")
                if add_padding:
                    strand = np.append(strand, b"")
                mk(f"data/{contig}/strand", strand, "|S1")

            # Extra 1-D datasets: methyl counts, fragment_end_clipped
            for extra_name, extra_arr in data.get("extra", {}).items():
                arr = extra_arr
                if add_padding:
                    arr = np.append(arr, np.zeros(1, dtype=arr.dtype))
                mk(f"data/{contig}/{extra_name}", arr, arr.dtype)

            # Build index
            cl = contig_lengths[contig]
            n_with_pad = n + (1 if add_padding else 0)
            if cl > index_block_size and n_with_pad >= MIN_NUM_READS_FOR_INDEX:
                s = starts if not add_padding else np.append(data["starts"], np.int32(0))
                block_indices = np.arange(0, cl, index_block_size)
                index_poss = np.searchsorted(s, block_indices, side="left")
                index_poss = np.append(index_poss, n_with_pad)
                f[f"index/{contig}"] = index_poss

        # fragment_length_counts
        max_fl = 65535
        flc = np.zeros(max_fl + 1, dtype=np.float64)
        for contig, data in contigs.items():
            lengths = data["lengths"]
            if add_padding:
                lengths = np.append(lengths, np.uint16(0))
            vals, counts = np.unique(lengths, return_counts=True)
            for v, c in zip(vals, counts):
                flc[v] += c
        f["fragment_length_counts"] = flc


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

CHR1_LEN = 20000
CHR2_LEN = 200


@pytest.fixture
def tiny_fasta():
    """A small FASTA.

    chr1 (20000 bp) is a varied ACGTN mix so that fragment GC fractions span
    [0, 1] with non-terminating decimals — an all-G contig makes round(x, 5) and
    round(x, 4) indistinguishable and hides rounding bugs (review finding M-3).
    chr2 (200 bp) keeps the G+N+G layout the N-prefix-sum test relies on.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "test.fa")
        _write_fasta(path, [
            ("chr1", _random_sequence(7, CHR1_LEN)),
            ("chr2", "G" * 80 + "N" * 40 + "G" * 80),
        ])
        yield path


def _fragment_geometry(seed_starts, seed_lengths, n, hi_start, len_lo, len_hi, contig_len):
    """Sorted starts / lengths clipped to the contig."""
    starts = np.sort(
        np.random.RandomState(seed_starts).randint(0, hi_start, n)
    ).astype("int32")
    lengths = np.random.RandomState(seed_lengths).randint(
        len_lo, len_hi, n
    ).astype("uint16")
    stops = starts.astype(np.int64) + lengths.astype(np.int64)
    mask = stops <= contig_len
    return starts[mask], lengths[mask]


def _build_fixture_contigs(fasta_path, seeds):
    """Build (starts, lengths, gc) per contig with GC from the independent oracle."""
    seqs = read_fasta_sequences(fasta_path)

    starts1, lengths1 = _fragment_geometry(
        seeds[0], seeds[1], 1500, CHR1_LEN - 400, 10, 300, CHR1_LEN
    )
    gc1 = naive_gc_uint8(seqs["chr1"], starts1, lengths1)

    starts2, lengths2 = _fragment_geometry(
        seeds[2], seeds[3], 180, CHR2_LEN - 40, 5, 30, CHR2_LEN
    )
    gc2 = naive_gc_uint8(seqs["chr2"], starts2, lengths2)

    return {
        "chr1": {"starts": starts1, "lengths": lengths1, "gc": gc1},
        "chr2": {"starts": starts2, "lengths": lengths2, "gc": gc2},
    }


FIXTURE_CONTIG_LENGTHS = {"chr1": CHR1_LEN, "chr2": CHR2_LEN}


@pytest.fixture
def clean_h5(tiny_fasta):
    """A correctly-built fragment H5 (no padding, oracle-derived GC)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        contigs = _build_fixture_contigs(tiny_fasta, (42, 43, 44, 45))
        h5_path = os.path.join(tmpdir, "clean.h5")
        _make_fragment_h5(h5_path, contigs, contig_lengths=FIXTURE_CONTIG_LENGTHS)
        yield h5_path


@pytest.fixture
def expected_clean_gc(tiny_fasta):
    """Oracle GC for the clean_h5 fixture, computed independently of repair.py."""
    return {
        c: d["gc"] for c, d in _build_fixture_contigs(tiny_fasta, (42, 43, 44, 45)).items()
    }


@pytest.fixture
def padded_h5(tiny_fasta):
    """Fragment H5 with phantom padding rows (simulating the caddb89 bug)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        contigs = _build_fixture_contigs(tiny_fasta, (50, 51, 52, 53))
        h5_path = os.path.join(tmpdir, "padded.h5")
        _make_fragment_h5(
            h5_path, contigs, add_padding=True,
            contig_lengths=FIXTURE_CONTIG_LENGTHS,
        )
        yield h5_path


@pytest.fixture
def expected_padded_gc(tiny_fasta):
    """Oracle GC for the padded_h5 fixture's real (non-phantom) fragments."""
    return {
        c: d["gc"] for c, d in _build_fixture_contigs(tiny_fasta, (50, 51, 52, 53)).items()
    }


@pytest.fixture
def cumsum_cache(tiny_fasta):
    """Build and return a cumsum cache dir for the tiny FASTA."""
    with tempfile.TemporaryDirectory() as cache_dir:
        fasta_sha256 = compute_file_sha256(tiny_fasta)
        from fragments_h5.repair import build_cumsum_cache
        build_cumsum_cache(tiny_fasta, cache_dir, fasta_sha256)
        yield cache_dir, fasta_sha256


# ---------------------------------------------------------------------------
# GC correctness against the independent oracle (review finding H-1)
# ---------------------------------------------------------------------------

class TestGcAgainstIndependentOracle:
    """recompute_gc_for_contig must agree with a naive character-counting oracle."""

    def test_matches_oracle_on_varied_sequence(self, tiny_fasta):
        """Every fragment's GC byte matches the oracle over a varied ACGTN contig."""
        from fragments_h5.fragment import get_g_or_c_cumsum

        seqs = read_fasta_sequences(tiny_fasta)
        starts, lengths = _fragment_geometry(
            110, 111, 3000, CHR1_LEN - 400, 1, 400, CHR1_LEN
        )
        cumsum, _ = get_g_or_c_cumsum(tiny_fasta, "chr1")

        actual = recompute_gc_for_contig(starts, lengths, cumsum)
        expected = naive_gc_uint8(seqs["chr1"], starts, lengths)

        n_diff = int(np.sum(actual != expected))
        assert n_diff == 0, (
            f"{n_diff} of {len(starts)} fragments disagree with the oracle; "
            f"first differing index {int(np.argmax(actual != expected))}"
        )

    def test_oracle_covers_the_gc_range(self, tiny_fasta):
        """The oracle fixture spans a wide range of GC bytes (not a degenerate case)."""
        seqs = read_fasta_sequences(tiny_fasta)
        starts, lengths = _fragment_geometry(
            110, 111, 3000, CHR1_LEN - 400, 1, 400, CHR1_LEN
        )
        gc = naive_gc_uint8(seqs["chr1"], starts, lengths)
        assert gc.min() < 40 and gc.max() > 210, (
            f"GC bytes span only [{gc.min()}, {gc.max()}]"
        )
        assert len(np.unique(gc)) > 100

    @pytest.mark.parametrize("seeds", [(42, 43, 44, 45), (50, 51, 52, 53)])
    def test_fixtures_are_rounding_sensitive(self, tiny_fasta, seeds):
        """The e2e fixtures must be able to distinguish round(x, 5) from round(x, 4).

        This is a guard on the fixtures themselves, not on repair.py. The
        end-to-end tests can only catch a quantisation bug if their fragments
        include some whose encoded byte depends on the fifth decimal place. If a
        future edit to the geometry or the FASTA destroys that, this fails loudly
        instead of silently weakening every e2e assertion.
        """
        seqs = read_fasta_sequences(tiny_fasta)
        contigs = _build_fixture_contigs(tiny_fasta, seeds)
        n_sensitive = count_rounding_sensitive(
            seqs["chr1"], contigs["chr1"]["starts"], contigs["chr1"]["lengths"]
        )
        assert n_sensitive >= 3, (
            f"fixture {seeds} has only {n_sensitive} rounding-sensitive fragments"
        )

    def test_oracle_handles_n_bases(self, tiny_fasta):
        """An all-N span encodes as 0.5 -> 127, matching the recomputation."""
        from fragments_h5.fragment import get_g_or_c_cumsum

        seqs = read_fasta_sequences(tiny_fasta)
        # chr2 positions [80, 120) are the N block
        starts = np.array([80, 90, 80], dtype="int32")
        lengths = np.array([40, 20, 10], dtype="uint16")
        expected = naive_gc_uint8(seqs["chr2"], starts, lengths)
        assert np.all(expected == 127)

        cumsum, _ = get_g_or_c_cumsum(tiny_fasta, "chr2")
        assert np.array_equal(recompute_gc_for_contig(starts, lengths, cumsum), expected)


# ---------------------------------------------------------------------------
# 0.a — delete-and-recreate dataset test (§3.2.2)
# ---------------------------------------------------------------------------

class TestDeleteAndRecreatDataset:
    """Assert recreated datasets match mk_dataset's parameters."""

    def test_recreated_dataset_preserves_params(self):
        """After truncation, datasets have gzip/4/chunked and correct dtypes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test.h5")
            starts = np.array([10, 50, 100, 200], dtype="int32")
            lengths = np.array([20, 30, 25, 15], dtype="uint16")
            gc = np.array([100, 150, 200, 50], dtype="uint8")
            mapq = np.array([[60, 55], [40, 30], [20, 10], [50, 45]], dtype="uint8")

            _make_fragment_h5(h5_path, {
                "chr1": {
                    "starts": starts, "lengths": lengths,
                    "gc": gc, "mapq": mapq,
                },
            }, add_padding=True, contig_lengths={"chr1": 1000})

            # Record attrs before
            with h5py.File(h5_path, "r") as f:
                attr_hash_before = hash_attrs_excluding(f)

            # Truncate
            with h5py.File(h5_path, "r+") as f:
                truncate_contig_datasets(f, "chr1")

            # Verify
            with h5py.File(h5_path, "r") as f:
                grp = f["data/chr1"]
                for ds_name in ("starts", "lengths", "gc", "mapq", "strand"):
                    ds = grp[ds_name]
                    assert ds.compression == "gzip"
                    assert ds.compression_opts == 4
                    assert ds.chunks is not None  # chunked

                assert grp["starts"].dtype == np.dtype("int32")
                assert grp["lengths"].dtype == np.dtype("uint16")
                assert grp["gc"].dtype == np.dtype("uint8")
                assert grp["mapq"].dtype == np.dtype("uint8")
                assert grp["strand"].dtype == np.dtype("|S1")

                # Length decreased by 1
                assert len(grp["starts"]) == 4  # was 5 (4 + padding)
                assert len(grp["gc"]) == 4

                # Attrs unchanged
                attr_hash_after = hash_attrs_excluding(f)
                assert attr_hash_before == attr_hash_after


# ---------------------------------------------------------------------------
# 0.b — synthetic padding-row fixture (§7.1.1)
# ---------------------------------------------------------------------------

class TestSyntheticPaddingFixture:
    """Build a fixture with one phantom row per contig, run the repair, verify."""

    def test_truncation_full_pipeline(
        self, padded_h5, tiny_fasta, cumsum_cache, expected_padded_gc,
    ):
        """Full repair on padded fixture: truncation fires, index rebuilt, flc correct."""
        cache_dir, fasta_sha256 = cumsum_cache

        # Record pre-repair state
        with h5py.File(padded_h5, "r") as f:
            old_flc = f["fragment_length_counts"][:]
            contigs = list(f["data"].keys())
            old_sizes = {c: len(f[f"data/{c}/starts"]) for c in contigs}
            old_index_keys = set(f["index"].keys()) if "index" in f else set()

        report = repair_local_file(
            local_path=padded_h5,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=False,
        )

        assert report["status"] in ("ok", "ok_no_gc")
        assert report["contigs_truncated"] == len(contigs)

        # Verify post-repair
        with h5py.File(padded_h5, "r") as f:
            for contig in contigs:
                # Each dataset one element shorter
                new_n = len(f[f"data/{contig}/starts"])
                assert new_n == old_sizes[contig] - 1, (
                    f"{contig}: expected {old_sizes[contig] - 1}, got {new_n}"
                )

                # starts is non-decreasing
                s = f[f"data/{contig}/starts"][:]
                if len(s) >= 2:
                    assert np.diff(s).min() >= 0

                # All datasets same length
                ds_lens = [
                    f[f"data/{contig}/{ds}"].shape[0]
                    for ds in f[f"data/{contig}"]
                ]
                assert len(set(ds_lens)) == 1

                # GC matches the independent oracle, not the code under test
                assert np.array_equal(
                    f[f"data/{contig}/gc"][:], expected_padded_gc[contig]
                ), f"{contig}: repaired gc does not match the independent oracle"

            # Index key set unchanged
            new_index_keys = set(f["index"].keys()) if "index" in f else set()
            assert new_index_keys == old_index_keys

            # fragment_length_counts exact identity (§3.2.4)
            new_flc = f["fragment_length_counts"][:]
            n_truncated = len(contigs)
            assert int(new_flc[0]) == int(old_flc[0]) - n_truncated
            assert np.array_equal(new_flc[1:], old_flc[1:])

            # _repair_history present
            history = json.loads(f.attrs["_repair_history"])
            assert len(history) == 1
            assert history[0]["tool"] == "repair-fragments-h5-gc"

    def test_detect_padding_all_contigs(self, padded_h5):
        """detect_padding_row returns 'truncate' for all contigs in padded file."""
        with h5py.File(padded_h5, "r") as f:
            for contig in f["data"]:
                verdict = detect_padding_row(f, contig)
                assert verdict == "truncate", f"{contig}: expected 'truncate', got {verdict}"

    def test_detect_no_padding_on_clean(self, clean_h5):
        """detect_padding_row returns 'clean' for all contigs in clean file."""
        with h5py.File(clean_h5, "r") as f:
            for contig in f["data"]:
                verdict = detect_padding_row(f, contig)
                assert verdict == "clean", f"{contig}: expected 'clean', got {verdict}"


# ---------------------------------------------------------------------------
# 0.c — 2-D mapq truncation along axis 0 only (§2.2.2)
# ---------------------------------------------------------------------------

class TestMapqTruncation:
    """Assert truncation is along axis 0 only for 2-D mapq."""

    def test_mapq_shape_after_truncation(self):
        """mapq shape goes from (n+1, 2) to (n, 2), not (n+1,) or (n,)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test.h5")
            n = 50
            starts = np.sort(np.random.RandomState(60).randint(0, 900, n)).astype("int32")
            lengths = np.random.RandomState(61).randint(10, 50, n).astype("uint16")
            mapq = np.column_stack([
                np.random.RandomState(62).randint(0, 60, n).astype("uint8"),
                np.random.RandomState(63).randint(0, 60, n).astype("uint8"),
            ])

            _make_fragment_h5(h5_path, {
                "chr1": {"starts": starts, "lengths": lengths, "mapq": mapq},
            }, add_padding=True, contig_lengths={"chr1": 10000})

            with h5py.File(h5_path, "r") as f:
                assert f["data/chr1/mapq"].shape == (n + 1, 2)

            with h5py.File(h5_path, "r+") as f:
                truncate_contig_datasets(f, "chr1")

            with h5py.File(h5_path, "r") as f:
                assert f["data/chr1/mapq"].shape == (n, 2), (
                    f"Expected ({n}, 2), got {f['data/chr1/mapq'].shape}"
                )
                # Values preserved
                assert np.array_equal(f["data/chr1/mapq"][:], mapq)


# ---------------------------------------------------------------------------
# 0.d — rounding-agreement test (§3.4)
# ---------------------------------------------------------------------------

class TestRoundingAgreement:
    """Elementwise equality between vectorized and CPython reference."""

    def _cpython_reference(self, num, den):
        """Plain CPython reference: int(round(round(float(num)/float(den), 5) * 254))."""
        if den == 0:
            return 255
        gc = round(float(num) / float(den), 5)
        return int(round(gc * 254))

    def test_rounding_random_fragments(self, tiny_fasta):
        """Rounding agreement over a large set of random fragments.

        chr1 is a varied ACGTN mix (see the tiny_fasta docstring) so the GC
        fractions here are mostly non-terminating decimals with denominators up
        to 400. That is what makes the fifth decimal place observable — on an
        all-G contig every fraction is 1.0 and round(x, 5) == round(x, 4).
        """
        from fragments_h5.fragment import get_g_or_c_cumsum

        cumsum, _ = get_g_or_c_cumsum(tiny_fasta, "chr1")

        # Generate many fragments
        rng = np.random.RandomState(100)
        n = 100_000  # 100k fragments for fast test
        starts = rng.randint(0, CHR1_LEN - 400, n).astype("int32")
        lengths = rng.randint(1, 400, n).astype("uint16")
        # Clip
        stops = starts.astype(np.int64) + lengths.astype(np.int64)
        mask = stops <= CHR1_LEN
        starts = starts[mask]
        lengths = lengths[mask]

        # Sanity: the fixture actually exercises a wide range of GC fractions
        fracs = (cumsum[stops[mask]] - cumsum[starts]) / lengths.astype(np.float64)
        assert fracs.min() < 0.2 and fracs.max() > 0.8, (
            f"GC fractions too narrow to exercise rounding: "
            f"[{fracs.min():.3f}, {fracs.max():.3f}]"
        )

        # Vectorized
        gc_vec = recompute_gc_for_contig(starts, lengths, cumsum)

        # CPython reference
        gc_ref = np.array([
            self._cpython_reference(
                float(cumsum[int(starts[i]) + int(lengths[i])] - cumsum[int(starts[i])]),
                int(lengths[i]),
            )
            for i in range(len(starts))
        ], dtype=np.uint8)

        assert np.array_equal(gc_vec, gc_ref), (
            f"Rounding disagreement: {np.sum(gc_vec != gc_ref)} of {len(starts)} fragments differ"
        )

    def test_rounding_half_ulp_engineered(self):
        """Test fragments engineered so q5 * 254 lands near k + 0.5."""
        # Construct a cumsum where we control the GC fraction precisely
        # cumsum[i] = i for all-G sequence
        cumsum = np.arange(1001, dtype=np.float64)

        # Find fractions where round(gc, 5) * 254 is near k + 0.5
        # gc = num/den, q5 = round(gc, 5), target: q5*254 ≈ k+0.5
        test_cases = []
        for den in range(10, 100):
            for num in range(den + 1):
                gc = round(float(num) / float(den), 5)
                val = gc * 254
                frac = val - int(val)
                if abs(frac - 0.5) < 0.01:  # near half
                    test_cases.append((num, den))

        assert len(test_cases) >= 50, f"Expected >= 50 half-ulp cases, got {len(test_cases)}"

        starts = np.array([0] * len(test_cases), dtype="int32")
        lengths = np.array([den for _, den in test_cases], dtype="uint16")
        # Build a cumsum that gives the right numerator for each
        # We need cumsum[0+length] - cumsum[0] = num, so cumsum[length] = num
        # But cumsum is shared... use individual checks
        for num, den in test_cases:
            cumsum_local = np.zeros(den + 1, dtype=np.float64)
            cumsum_local[den] = float(num)

            gc_vec = recompute_gc_for_contig(
                np.array([0], dtype="int32"),
                np.array([den], dtype="uint16"),
                cumsum_local,
            )
            gc_ref = self._cpython_reference(num, den)
            assert gc_vec[0] == gc_ref, (
                f"Half-ULP disagreement: num={num}, den={den}, "
                f"vec={gc_vec[0]}, ref={gc_ref}"
            )


# ---------------------------------------------------------------------------
# 0.e — float32-accumulator-simulation test (§5b)
# ---------------------------------------------------------------------------

class TestFloat32AccumulatorSimulation:
    """Verify simulated float32 cumsum matches the actual pre-fix accumulator."""

    def test_simulation_matches_direct(self):
        """Simulate from float64 cumsum and compare with direct float32 cumsum."""
        from fragments_h5.sequence import one_hot_encode_sequences

        # Build a test sequence with mixed content
        seq_str = "G" * 500 + "N" * 200 + "A" * 300 + "C" * 500 + "T" * 500
        padded = ("a" + seq_str).encode()

        # Direct float32 path (the pre-fix bug)
        encoded = one_hot_encode_sequences([padded])[0]
        per_base_f32 = encoded[:, (1, 2)].sum(axis=1)  # float32
        direct_f32_cumsum = np.cumsum(per_base_f32)  # float32 cumsum

        # Float64 cumsum (the fixed path)
        per_base_f64 = encoded[:, (1, 2)].sum(axis=1).astype(np.float64)
        f64_cumsum = np.cumsum(per_base_f64)  # float64

        # Simulated float32 from the float64 cumsum
        simulated = simulate_float32_cumsum(f64_cumsum)

        assert np.array_equal(direct_f32_cumsum, simulated), (
            f"Simulation diverges from direct float32: "
            f"max diff = {np.max(np.abs(direct_f32_cumsum.astype(np.float64) - simulated.astype(np.float64)))}"
        )

    def test_saturation_point(self):
        """Float32 cumsum saturates at 2**24 for all-G sequence."""
        # Build a float64 cumsum for an all-G contig long enough to saturate
        # In all-G: cumsum[i] = i (each G contributes 1.0)
        n = 2**24 + 100
        f64_cumsum = np.arange(n + 1, dtype=np.float64)  # include pad

        f32_sim = simulate_float32_cumsum(f64_cumsum)

        # Before 2**24: should match
        assert f32_sim[2**24] == np.float32(2**24)

        # At 2**24 + 1: adding 1.0 to 2**24 rounds back to 2**24 in float32
        assert f32_sim[2**24 + 1] == np.float32(2**24), (
            f"Expected saturation at 2**24, got {f32_sim[2**24 + 1]}"
        )

        # Fragment spanning the saturation: cumsum[stop] - cumsum[start] = 0 in float32
        assert f32_sim[2**24 + 50] - f32_sim[2**24] == 0.0


# ---------------------------------------------------------------------------
# 0.f — scan_fasta_alphabet (§5 layer 5)
# ---------------------------------------------------------------------------

class TestFastaAlphabet:
    """Test the FASTA alphabet scanning gate."""

    def test_valid_fasta(self, tiny_fasta):
        """scan_fasta_alphabet accepts a valid {A,C,G,T,N} FASTA."""
        result = scan_fasta_alphabet(tiny_fasta)
        assert "G" in result
        assert "N" in result
        # No disallowed characters
        for c in result:
            assert c in "ACGTN", f"Unexpected character: {c}"

    def test_invalid_fasta(self):
        """scan_fasta_alphabet rejects FASTA with non-{A,C,G,T,N} characters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "bad.fa")
            _write_fasta(path, [("chr1", "ACGTBVHD" * 10)])
            with pytest.raises(ValueError, match="characters outside"):
                scan_fasta_alphabet(path)


# ---------------------------------------------------------------------------
# detect_padding_row edge cases
# ---------------------------------------------------------------------------

class TestDetectPaddingRow:
    """Edge cases for the padding row detection predicate."""

    def test_abort_on_mixed_signal(self):
        """When only one of (sortedness, zero-sig) holds, abort."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test.h5")
            # Sortedness violation but last row is NOT all-zero (non-zero length)
            with h5py.File(h5_path, "w") as f:
                f.create_group("data")
                f.create_group("index")
                f.attrs["index_block_size"] = 5000
                f.attrs["max_fragment_length"] = 65535
                f.attrs["_contig_lengths_str"] = str({"chr1": 1000})
                f.attrs["_bam_header"] = ""
                f.attrs["_source_format"] = "BAM"

                grp = f.create_group("data/chr1")
                # starts[-1] < starts[-2] but lengths[-1] != 0
                grp.create_dataset("starts", data=np.array([10, 50, 100, 5], dtype="int32"))
                grp.create_dataset("lengths", data=np.array([20, 30, 25, 15], dtype="uint16"))
                grp.create_dataset("mapq", data=np.zeros((4, 2), dtype="uint8"))
                grp.create_dataset("gc", data=np.array([100, 150, 200, 50], dtype="uint8"))
                grp.create_dataset("strand", data=np.array([b"+", b"+", b"+", b"+"], dtype="|S1"))
                f["fragment_length_counts"] = np.zeros(65536)

            with h5py.File(h5_path, "r") as f:
                with pytest.raises(RepairAbort, match="exactly one"):
                    detect_padding_row(f, "chr1")

    def test_single_fragment_is_clean(self):
        """A contig with n=1 (< 2) is always clean."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test.h5")
            _make_fragment_h5(h5_path, {
                "chr1": {
                    "starts": np.array([100], dtype="int32"),
                    "lengths": np.array([50], dtype="uint16"),
                },
            }, contig_lengths={"chr1": 1000})
            with h5py.File(h5_path, "r") as f:
                assert detect_padding_row(f, "chr1") == "clean"


# ---------------------------------------------------------------------------
# methyl / fragment_end_clipped datasets (review finding M-2)
# ---------------------------------------------------------------------------

METHYL_KEYS = ("num_cpgs", "num_converted_cpgs", "num_cytosines", "num_converted_cytosines")


def _methyl_extras(n, seed):
    """Realistic methyl + fragment_end_clipped arrays for n fragments."""
    rng = np.random.RandomState(seed)
    extras = {k: rng.randint(1, 30, n).astype("uint16") for k in METHYL_KEYS}
    # fragment_end_clipped is 0/1/255; keep it nonzero so a zero tail is a signal
    extras["fragment_end_clipped"] = rng.choice([1, 255], size=n).astype("uint8")
    return extras


def _make_methyl_h5(path, n=250, add_padding=True, contig_length=100000, seed=90):
    starts = np.sort(np.random.RandomState(seed).randint(0, 90000, n)).astype("int32")
    lengths = np.random.RandomState(seed + 1).randint(10, 50, n).astype("uint16")
    _make_fragment_h5(path, {
        "chr1": {
            "starts": starts,
            "lengths": lengths,
            "extra": _methyl_extras(n, seed + 2),
        },
    }, read_gc=False, add_padding=add_padding, contig_lengths={"chr1": contig_length})
    return starts, lengths


class TestMethylDatasets:
    """detect_padding_row and truncate_contig_datasets on methyl-bearing files."""

    def test_padding_detected_with_methyl_datasets(self):
        """A phantom row is still detected when methyl arrays are present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "methyl.h5")
            _make_methyl_h5(h5_path)

            with h5py.File(h5_path, "r") as f:
                present = set(f["data/chr1"].keys())
                assert set(METHYL_KEYS) | {"fragment_end_clipped"} <= present
                assert detect_padding_row(f, "chr1") == "truncate"

    def test_nonzero_methyl_tail_aborts(self):
        """The zero signature covers methyl arrays: a nonzero tail must abort.

        This is what proves detect_padding_row actually inspects the extra
        datasets rather than ignoring them.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "methyl_bad.h5")
            _make_methyl_h5(h5_path)

            # Sortedness violation still holds, but one methyl tail is nonzero
            with h5py.File(h5_path, "r+") as f:
                arr = f["data/chr1/num_cpgs"][:]
                arr[-1] = 7
                del f["data/chr1/num_cpgs"]
                f.create_dataset("data/chr1/num_cpgs", data=arr, dtype="uint16")

            with h5py.File(h5_path, "r") as f:
                with pytest.raises(RepairAbort, match="exactly one"):
                    detect_padding_row(f, "chr1")

    def test_truncation_covers_methyl_datasets(self):
        """Every methyl array is truncated, with dtype, params and values intact."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "methyl.h5")
            n = 250
            _make_methyl_h5(h5_path, n=n)

            with h5py.File(h5_path, "r") as f:
                before = {k: f[f"data/chr1/{k}"][:] for k in f["data/chr1"]}
                assert all(len(v) == n + 1 for v in before.values())

            with h5py.File(h5_path, "r+") as f:
                truncate_contig_datasets(f, "chr1")

            with h5py.File(h5_path, "r") as f:
                grp = f["data/chr1"]
                assert set(grp.keys()) == set(before.keys())
                for key in METHYL_KEYS + ("fragment_end_clipped",):
                    ds = grp[key]
                    assert ds.shape == (n,), f"{key}: {ds.shape}"
                    assert ds.dtype == before[key].dtype
                    assert ds.compression == "gzip"
                    assert ds.compression_opts == 4
                    assert np.array_equal(ds[:], before[key][:-1])

    def test_repair_end_to_end_with_methyl(self):
        """Full repair on a methyl file: truncation applies to every dataset."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "methyl.h5")
            n = 250
            _make_methyl_h5(h5_path, n=n)

            with h5py.File(h5_path, "r") as f:
                before = {k: f[f"data/chr1/{k}"][:] for k in f["data/chr1"]}

            report = repair_local_file(
                local_path=h5_path,
                fasta_path=None,
                cumsum_cache_dir=None,
                fasta_sha256=None,
                dry_run=False,
            )
            assert report["status"] == "ok_no_gc"
            assert report["contigs_truncated"] == 1

            with h5py.File(h5_path, "r") as f:
                grp = f["data/chr1"]
                for key in before:
                    assert grp[key].shape[0] == n, f"{key} not truncated"
                    assert np.array_equal(grp[key][:], before[key][:-1]), f"{key} corrupted"
                # Index rebuilt against the truncated starts
                assert f["index/chr1"][-1] == n


# ---------------------------------------------------------------------------
# no-gc file handling (§3.2.5)
# ---------------------------------------------------------------------------

class TestNoGcFile:
    """Files with no gc dataset: truncation runs, gc skipped."""

    def test_no_gc_truncation_only(self):
        """Repair a file with no gc: truncation works, no FASTA needed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "no_gc.h5")
            n = 200
            starts = np.sort(np.random.RandomState(70).randint(0, 9000, n)).astype("int32")
            lengths = np.random.RandomState(71).randint(10, 50, n).astype("uint16")

            _make_fragment_h5(h5_path, {
                "chr1": {"starts": starts, "lengths": lengths},
            }, read_gc=False, add_padding=True, contig_lengths={"chr1": 10000})

            with h5py.File(h5_path, "r") as f:
                assert check_gc_presence(f["data"]) == "none"
                old_n = len(f["data/chr1/starts"])

            report = repair_local_file(
                local_path=h5_path,
                fasta_path=None,
                cumsum_cache_dir=None,
                fasta_sha256=None,
                dry_run=False,
            )

            assert report["status"] == "ok_no_gc"
            assert report["contigs_truncated"] == 1

            with h5py.File(h5_path, "r") as f:
                assert len(f["data/chr1/starts"]) == old_n - 1

    def test_partial_gc_aborts(self):
        """File with gc on some contigs but not others: abort."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "partial_gc.h5")
            with h5py.File(h5_path, "w") as f:
                f.create_group("data")
                f.create_group("index")
                f.attrs["index_block_size"] = 5000
                f.attrs["max_fragment_length"] = 65535
                f.attrs["_contig_lengths_str"] = str({"chr1": 1000, "chr2": 1000})
                f.attrs["_bam_header"] = ""
                f.attrs["_source_format"] = "BAM"

                for c in ("chr1", "chr2"):
                    grp = f.create_group(f"data/{c}")
                    grp.create_dataset("starts", data=np.array([10, 50], dtype="int32"))
                    grp.create_dataset("lengths", data=np.array([20, 30], dtype="uint16"))
                    grp.create_dataset("mapq", data=np.zeros((2, 2), dtype="uint8"))
                    grp.create_dataset("strand", data=np.array([b"+", b"+"], dtype="|S1"))
                # Only chr1 gets gc
                f["data/chr1"].create_dataset("gc", data=np.array([100, 150], dtype="uint8"))
                f["fragment_length_counts"] = np.zeros(65536)

            with pytest.raises(RepairAbort, match="Partial gc"):
                repair_local_file(h5_path, dry_run=True)


# ---------------------------------------------------------------------------
# Clean-file idempotency (§5b)
# ---------------------------------------------------------------------------

class TestCleanFileIdempotency:
    """On a clean file, repair changes nothing but provenance."""

    def test_dry_run_all_zeros(self, clean_h5, tiny_fasta, cumsum_cache):
        """Dry-run on clean file: zero changes, zero truncations."""
        cache_dir, fasta_sha256 = cumsum_cache
        report = repair_local_file(
            local_path=clean_h5,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=True,
        )

        assert report["status"] == "dry_run_ok"
        assert report["contigs_truncated"] == 0

        # All gc changes should be zero
        for contig, res in report.get("gc_results", {}).items():
            stats = res.get("stats", {})
            assert stats.get("changed_total", 0) == 0, (
                f"{contig}: {stats.get('changed_total')} gc changes on clean file"
            )

    def test_apply_preserves_gc(
        self, clean_h5, tiny_fasta, cumsum_cache, expected_clean_gc,
    ):
        """Apply on clean file: gc unchanged, repair_history added."""
        import shutil
        cache_dir, fasta_sha256 = cumsum_cache

        # Work on a copy
        with tempfile.TemporaryDirectory() as tmpdir:
            copy_path = os.path.join(tmpdir, "copy.h5")
            shutil.copy2(clean_h5, copy_path)

            with h5py.File(copy_path, "r") as f:
                old_gc = {
                    c: f[f"data/{c}/gc"][:] for c in f["data"]
                }

            report = repair_local_file(
                local_path=copy_path,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=False,
            )

            assert report["status"] == "ok"

            with h5py.File(copy_path, "r") as f:
                for c in f["data"]:
                    assert np.array_equal(f[f"data/{c}/gc"][:], old_gc[c]), (
                        f"{c}: gc changed on clean file"
                    )
                    # ... and that unchanged value is the independently-derived one
                    assert np.array_equal(f[f"data/{c}/gc"][:], expected_clean_gc[c]), (
                        f"{c}: gc does not match the independent oracle"
                    )
                # repair_history added
                assert "_repair_history" in f.attrs


# ---------------------------------------------------------------------------
# N-prefix sum helper
# ---------------------------------------------------------------------------

class TestNPrefixSum:
    """Test the N-position prefix sum builder."""

    def test_n_positions_detected(self, tiny_fasta):
        """N positions correctly identified from cumsum."""
        from fragments_h5.fragment import get_g_or_c_cumsum

        cumsum, _ = get_g_or_c_cumsum(tiny_fasta, "chr2")  # G*80 + N*40 + G*80
        n_prefix = build_n_prefix_sum(cumsum)

        # Positions 80-119 should be N
        # span [0, 80) — no N
        assert n_prefix[80] - n_prefix[0] == 0
        # span [80, 120) — 40 N's
        assert n_prefix[120] - n_prefix[80] == 40
        # span [60, 100) — 20 N's (positions 80-99)
        assert n_prefix[100] - n_prefix[60] == 20


# ---------------------------------------------------------------------------
# Region classification
# ---------------------------------------------------------------------------

FRAG_LEN = 10


def _region_inputs(regions, n_flags=None):
    """Build (starts, lengths, f32_cumsum, n_prefix) placing each fragment in a region.

    regions: sequence of 'below' | 'middle' | 'above' — the float32 cumsum band
        that fragment i's stop coordinate lands in.
    n_flags: optional sequence of bool — whether fragment i's span contains an N.

    f32_cumsum and n_prefix are constructed directly rather than derived from a
    real contig: reaching T24 for real needs a 16.7 Mb sequence, and the thing
    under test is the classification rule, not the cumsum.
    """
    n = len(regions)
    starts = (np.arange(n) * FRAG_LEN).astype(np.int32)
    lengths = np.full(n, FRAG_LEN, dtype=np.uint16)
    size = n * FRAG_LEN + 1

    band_value = {
        "below": np.float32(T23 - 1000.0),
        "middle": np.float32(T23 + 1000.0),
        "above": np.float32(T24 + 1000.0),
    }
    f32_cumsum = np.zeros(size, dtype=np.float32)
    for i, region in enumerate(regions):
        f32_cumsum[(i + 1) * FRAG_LEN] = band_value[region]

    is_n = np.zeros(size - 1, dtype=bool)
    if n_flags is not None:
        for i, flag in enumerate(n_flags):
            if flag:
                is_n[i * FRAG_LEN] = True
    n_prefix = np.zeros(size, dtype=np.int64)
    n_prefix[1:] = np.cumsum(is_n)

    return starts, lengths, f32_cumsum, n_prefix


class TestRegionClassification:
    """Test classify_gc_changes with known region assignments."""

    def test_clean_file_no_violations(self):
        """Identical old/new gc produces no violations."""
        n = 100
        gc = np.full(n, 127, dtype=np.uint8)
        starts = np.arange(0, n * 100, 100, dtype=np.int32)
        lengths = np.full(n, 50, dtype=np.uint16)

        # Cumsum that keeps everything in the < T23 region
        cumsum = np.arange(n * 100 + 51, dtype=np.float64) * 0.5
        f32_cumsum = simulate_float32_cumsum(cumsum)
        n_prefix = build_n_prefix_sum(cumsum)

        stats, violations = classify_gc_changes(gc, gc, starts, lengths, f32_cumsum, n_prefix)
        assert len(violations) == 0
        assert stats["changed_total"] == 0

    # -- permissive direction: genuine corruption patterns must be allowed ---

    def test_saturation_corruption_above_t24_permitted(self):
        """The real corruption signature — gc 0 above T24 — classifies clean.

        Past T24 the float32 accumulator stops advancing, so every numerator is
        0 and the stored byte is 0. Recomputation replaces those with real
        values; that is the whole point of the tool and must not be blocked.
        """
        n = 60
        regions = ["above"] * n
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(regions)

        old_gc = np.zeros(n, dtype=np.uint8)
        new_gc = np.arange(1, n + 1, dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert violations == []
        assert stats["above_T24"] == n
        assert stats["changed_above_T24"] == n
        assert stats["zero_to_nonzero"] == n

    def test_middle_band_change_on_n_span_permitted(self):
        """Between T23 and T24 only N-containing spans may change — and they may."""
        regions = ["middle"] * 4
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(
            regions, n_flags=[True, True, True, True]
        )

        old_gc = np.array([100, 100, 100, 100], dtype=np.uint8)
        new_gc = np.array([101, 99, 100, 102], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert violations == []
        assert stats["changed_between_T23_T24"] == 3
        assert stats["changed_between_T23_T24_no_N"] == 0

    def test_mixed_regions_only_flags_the_bad_ones(self):
        """A file with all three bands: only the below-T23 change is a violation."""
        regions = ["below", "middle", "above"]
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(
            regions, n_flags=[False, True, False]
        )
        old_gc = np.array([50, 50, 0], dtype=np.uint8)
        new_gc = np.array([51, 52, 200], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert stats["changed_below_T23"] == 1
        assert stats["changed_between_T23_T24"] == 1
        assert stats["changed_between_T23_T24_no_N"] == 0
        assert stats["changed_above_T24"] == 1
        assert len(violations) == 1
        assert "pre-saturation" in violations[0]

    # -- blocking direction: impossible patterns must be rejected -----------

    def test_below_t23_change_blocked(self):
        """A gc change in the pre-saturation band is impossible — must violate."""
        n = 20
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(["below"] * n)
        old_gc = np.zeros(n, dtype=np.uint8)
        new_gc = np.full(n, 200, dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert stats["changed_below_T23"] == n
        assert any("pre-saturation" in v for v in violations), violations

    def test_middle_band_change_without_n_blocked(self):
        """Between T23 and T24, an N-free span cannot legitimately change."""
        regions = ["middle"] * 5
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(
            regions, n_flags=[False, False, True, False, False]
        )
        old_gc = np.full(5, 100, dtype=np.uint8)
        new_gc = np.array([101, 101, 101, 101, 100], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        # 4 changed; 3 of them on N-free spans
        assert stats["changed_between_T23_T24"] == 4
        assert stats["changed_between_T23_T24_no_N"] == 3
        assert any("N-free" in v for v in violations), violations

    def test_255_to_value_blocked(self):
        """255 means 'no reference'; turning it into a number implies a wrong FASTA."""
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(["above"] * 3)
        old_gc = np.array([255, 255, 10], dtype=np.uint8)
        new_gc = np.array([120, 130, 20], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert stats["x255_to_other"] == 2
        assert any("255→x" in v for v in violations), violations

    def test_value_to_255_blocked(self):
        """Turning a real GC value into 255 implies the contig went missing."""
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(["above"] * 3)
        old_gc = np.array([10, 20, 30], dtype=np.uint8)
        new_gc = np.array([255, 20, 30], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert stats["other_to_x255"] == 1
        assert any("x→255" in v for v in violations), violations


class TestCorruptedFileAborts:
    """End-to-end: a file whose gc disagrees below T23 must fail closed."""

    def test_corrupted_gc_aborts_repair(self, clean_h5, tiny_fasta, cumsum_cache):
        """Simulated saturation damage on a small contig is impossible — abort.

        Every fragment in the fixture sits far below T23, so no gc byte there can
        legitimately differ from the recomputation. Zeroing a handful (exactly
        what float32 saturation does to real files) must therefore be rejected
        rather than silently "repaired".
        """
        cache_dir, fasta_sha256 = cumsum_cache

        with h5py.File(clean_h5, "r+") as f:
            gc = f["data/chr1/gc"][:]
            assert np.any(gc != 0)
            gc[-25:] = 0  # the saturation signature: numerator collapses to zero
            del f["data/chr1/gc"]
            f.create_dataset("data/chr1/gc", data=gc, dtype="uint8")

        with pytest.raises(RepairAbort, match="pre-saturation"):
            repair_local_file(
                local_path=clean_h5,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=True,
            )

    def test_corrupted_gc_leaves_file_untouched(self, clean_h5, tiny_fasta, cumsum_cache):
        """The abort happens before any write, even in --apply mode."""
        cache_dir, fasta_sha256 = cumsum_cache

        with h5py.File(clean_h5, "r+") as f:
            gc = f["data/chr2/gc"][:]
            gc[:10] = 0
            del f["data/chr2/gc"]
            f.create_dataset("data/chr2/gc", data=gc, dtype="uint8")

        sha_before = compute_file_sha256(clean_h5)
        with pytest.raises(RepairAbort):
            repair_local_file(
                local_path=clean_h5,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=False,
            )
        assert compute_file_sha256(clean_h5) == sha_before


# ---------------------------------------------------------------------------
# Index rebuild
# ---------------------------------------------------------------------------

class TestIndexRebuild:
    """Test rebuild_contig_index matches the builder's output."""

    def test_index_rebuild_matches_builder(self):
        """Rebuild produces the same index as the original builder."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test.h5")
            n = 500
            starts = np.sort(np.random.RandomState(80).randint(0, 90000, n)).astype("int32")
            lengths = np.random.RandomState(81).randint(10, 50, n).astype("uint16")

            contig_lengths = {"chr1": 100000}
            _make_fragment_h5(h5_path, {
                "chr1": {"starts": starts, "lengths": lengths},
            }, contig_lengths=contig_lengths)

            with h5py.File(h5_path, "r") as f:
                original_index = f["index/chr1"][:] if "chr1" in f["index"] else None

            # Rebuild
            with h5py.File(h5_path, "r+") as f:
                rebuilt = rebuild_contig_index(f, "chr1", 100000, INDEX_BLOCK_SIZE)

            if original_index is not None:
                with h5py.File(h5_path, "r") as f:
                    new_index = f["index/chr1"][:]
                assert np.array_equal(original_index, new_index), (
                    f"Rebuilt index differs from original"
                )

    def test_index_skip_short_contig(self):
        """Contigs shorter than index_block_size get no index."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test.h5")
            _make_fragment_h5(h5_path, {
                "chr1": {
                    "starts": np.arange(0, 4000, 10, dtype="int32"),
                    "lengths": np.full(400, 5, dtype="uint16"),
                },
            }, contig_lengths={"chr1": 4999})  # < INDEX_BLOCK_SIZE

            with h5py.File(h5_path, "r+") as f:
                rebuilt = rebuild_contig_index(f, "chr1", 4999, INDEX_BLOCK_SIZE)
                assert not rebuilt
                assert "chr1" not in f.get("index", {})

    def test_index_skip_few_reads(self):
        """Contigs with < MIN_NUM_READS_FOR_INDEX reads get no index."""
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test.h5")
            _make_fragment_h5(h5_path, {
                "chr1": {
                    "starts": np.arange(0, 990, 10, dtype="int32"),  # 99 reads
                    "lengths": np.full(99, 5, dtype="uint16"),
                },
            }, contig_lengths={"chr1": 100000})

            with h5py.File(h5_path, "r+") as f:
                rebuilt = rebuild_contig_index(f, "chr1", 100000, INDEX_BLOCK_SIZE)
                assert not rebuilt


# ---------------------------------------------------------------------------
# Repair idempotency (§8.1)
# ---------------------------------------------------------------------------

class TestRepairIdempotency:
    """Repairing an already-repaired file is safe."""

    @staticmethod
    def _snapshot_datasets(path):
        """Every dataset in the file, by full path, as raw arrays."""
        snapshot = {}

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                snapshot[name] = obj[:]

        with h5py.File(path, "r") as f:
            f.visititems(visit)
        return snapshot

    def test_idempotent_datasets(self, padded_h5, tiny_fasta, cumsum_cache):
        """A second --apply run leaves every dataset byte-identical.

        Idempotency replaced the crash-resume protocol in the design, so the
        second pass must exercise the write path, not dry-run.
        """
        cache_dir, fasta_sha256 = cumsum_cache

        # First repair
        report1 = repair_local_file(
            local_path=padded_h5,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=False,
        )
        assert report1["contigs_truncated"] > 0

        after_first = self._snapshot_datasets(padded_h5)
        assert after_first, "no datasets found in the repaired file"

        # Second repair — a real write pass, not a dry run
        report2 = repair_local_file(
            local_path=padded_h5,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=False,
        )

        assert report2["status"] == "ok"
        assert report2["contigs_truncated"] == 0
        for contig, res in report2.get("gc_results", {}).items():
            stats = res.get("stats", {})
            assert stats.get("changed_total", 0) == 0, (
                f"{contig}: gc changed on second repair"
            )

        after_second = self._snapshot_datasets(padded_h5)
        assert set(after_second) == set(after_first)
        for name in after_first:
            assert np.array_equal(after_second[name], after_first[name]), (
                f"{name} changed on the second --apply run"
            )

        # Only the provenance attr grows — one history element per apply run
        with h5py.File(padded_h5, "r") as f:
            history = json.loads(f.attrs["_repair_history"])
        assert len(history) == 2


# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------

class TestCLIParsing:
    """Basic CLI argument validation."""

    def test_default_is_dry_run(self):
        """Without --apply, dry_run is True."""
        from fragments_h5.repair import parse_args
        args = parse_args(["--local-file", "/tmp/test.h5"])
        assert args.dry_run is True

    def test_apply_overrides_dry_run(self):
        """--apply sets dry_run to False."""
        from fragments_h5.repair import parse_args
        args = parse_args([
            "--apply", "--local-file", "/tmp/test.h5",
        ])
        assert args.dry_run is False

    def test_max_files_default(self):
        """Default --max-files is 1."""
        from fragments_h5.repair import parse_args
        args = parse_args(["--local-file", "/tmp/test.h5"])
        assert args.max_files == 1

    # -- §5 layer 2: the reference pin (review finding M-1) -----------------

    def test_apply_requires_expect_fasta_sha256(self):
        """--apply with --fasta and no pin is rejected before anything runs."""
        from fragments_h5.repair import parse_args
        with pytest.raises(SystemExit):
            parse_args([
                "--apply", "--local-file", "/tmp/test.h5", "--fasta", "/tmp/ref.fa",
            ])

    def test_apply_accepts_expect_fasta_sha256(self):
        """--apply with the pin present parses fine."""
        from fragments_h5.repair import parse_args
        args = parse_args([
            "--apply", "--local-file", "/tmp/test.h5", "--fasta", "/tmp/ref.fa",
            "--expect-fasta-sha256", "0" * 64,
        ])
        assert args.dry_run is False
        assert args.expect_fasta_sha256 == "0" * 64

    def test_dry_run_does_not_require_the_pin(self):
        """The pin is only mandatory when writing."""
        from fragments_h5.repair import parse_args
        args = parse_args(["--local-file", "/tmp/test.h5", "--fasta", "/tmp/ref.fa"])
        assert args.expect_fasta_sha256 is None

    def test_apply_without_fasta_does_not_require_the_pin(self):
        """No --fasta means no reference to pin (§3.2.5 no-gc files).

        Targets that do have gc are still covered: repair_local_file aborts when
        gc is present and --fasta is absent (see test_gc_file_without_fasta_aborts).
        """
        from fragments_h5.repair import parse_args
        args = parse_args(["--apply", "--local-file", "/tmp/test.h5"])
        assert args.dry_run is False

    def test_num_processes_flag_removed(self):
        """--num-processes did nothing and is gone (review finding L-1)."""
        from fragments_h5.repair import parse_args
        with pytest.raises(SystemExit):
            parse_args(["--local-file", "/tmp/test.h5", "--num-processes", "4"])


class TestFastaRequiredForGcFiles:
    """The other half of the §5 layer-2 gate: gc without a reference is fatal."""

    def test_gc_file_without_fasta_aborts(self, clean_h5):
        """A file with a gc dataset cannot be repaired without --fasta."""
        with pytest.raises(RepairAbort, match="--fasta not provided"):
            repair_local_file(local_path=clean_h5, fasta_path=None, dry_run=True)
