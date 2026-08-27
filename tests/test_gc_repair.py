"""Tests for the GC repair tool — stage-0 checklist items (§7.4).

Covers:
  0.a — delete-and-recreate dataset test
  0.b — synthetic padding-row fixture with full check list
  0.c — 2-D mapq truncation along axis 0 only
  0.d — rounding agreement: exhaustive over the quantiser domain, plus random
        real-contig fragments, plus the engineered half-ulp set
  0.e — float32-accumulator-simulation test
  0.f — (alphabet histogram regeneration — tested via scan_fasta_alphabet)
  0.g — existing overflow regression test still passes (separate file)
  0.h — T23/T24 end-to-end integration over a real ~17.2 Mb saturating contig
        (marked slow)

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
    find_t23_crossing,
    hash_attrs_excluding,
    predict_index_contigs,
    recompute_gc_for_contig,
    rebuild_contig_index,
    should_index_contig,
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


# A contig that lives in the H5 but not in the FASTA. Real files carry ~165 of
# these (unplaced/unlocalised scaffolds) because the BAM header lists contigs
# the GC reference does not. Both dimensions are deliberately under the index
# thresholds — length < INDEX_BLOCK_SIZE and n < MIN_NUM_READS_FOR_INDEX — so
# the contig never carries an index and the index key-set check is unaffected.
SCAFFOLD_NAME = "chrUn_scaffold99"
SCAFFOLD_LEN = 3000
SCAFFOLD_N_FRAGS = 40


def _scaffold_contig(n_non_255=0):
    """Fragments on a FASTA-absent contig; gc is 255 ('no reference') by default."""
    starts = np.sort(
        np.random.RandomState(300).randint(0, SCAFFOLD_LEN - 300, SCAFFOLD_N_FRAGS)
    ).astype("int32")
    lengths = np.random.RandomState(301).randint(
        10, 200, SCAFFOLD_N_FRAGS
    ).astype("uint16")
    gc = np.full(SCAFFOLD_N_FRAGS, 255, dtype="uint8")
    gc[:n_non_255] = 128
    return {"starts": starts, "lengths": lengths, "gc": gc}


@pytest.fixture
def absent_contig_h5(tiny_fasta):
    """Factory building a fragment H5 whose FASTA lacks one of its contigs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        counter = [0]

        def build(n_non_255=0):
            contigs = _build_fixture_contigs(tiny_fasta, (42, 43, 44, 45))
            contigs[SCAFFOLD_NAME] = _scaffold_contig(n_non_255=n_non_255)
            counter[0] += 1
            path = os.path.join(tmpdir, f"absent_{counter[0]}.h5")
            _make_fragment_h5(
                path, contigs,
                contig_lengths=dict(
                    FIXTURE_CONTIG_LENGTHS, **{SCAFFOLD_NAME: SCAFFOLD_LEN}
                ),
            )
            return path

        yield build


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

    def test_rounding_exhaustive_over_quantiser_domain(self):
        """Every (numerator, length) pair with length <= 400, exhaustively.

        The quantiser's output is a pure function of (num, den). The numerator
        `cumsum[stop] - cumsum[start]` is a half-integer, and on any real contig
        it stays under ~1e8 — far inside the range where float64 represents
        half-integers exactly — so the subtraction is exact and two fragments
        sharing a (num, den) always produce the same byte.

        That makes random fragments a resampling of a bounded set of equivalence
        classes rather than a source of new cases: 10^5 draws from the tiny_fasta
        contig reach 22894 distinct classes and 10^6 reach only 32412, for ten
        times the runtime. This enumerates all 160800 classes with den <= 400,
        which strictly dominates any number of random draws over the same band.
        See design §3.4 for why this replaced the original ">= 10^7 fragments".
        """
        max_den = 400
        n_classes = 0
        for den in range(1, max_den + 1):
            # Every reachable numerator: 0, 0.5, 1.0, ... den (N contributes 0.5).
            numerators = np.arange(0, 2 * den + 1, dtype=np.float64) * 0.5
            k = len(numerators)

            # Lay the classes out back to back with stride `den`, so fragment i
            # is [i*den, (i+1)*den) and its numerator is numerators[i].
            cumsum = np.zeros(k * den + 1, dtype=np.float64)
            cumsum[den::den] = np.cumsum(numerators)
            starts = (np.arange(k) * den).astype("int32")
            lengths = np.full(k, den, dtype="uint16")
            assert np.array_equal(
                cumsum[starts.astype(np.int64) + den] - cumsum[starts], numerators
            ), f"den={den}: synthetic cumsum does not reproduce the numerators"

            actual = recompute_gc_for_contig(starts, lengths, cumsum)
            expected = np.array(
                [self._cpython_reference(num, den) for num in numerators],
                dtype=np.uint8,
            )
            n_bad = int(np.sum(actual != expected))
            assert n_bad == 0, (
                f"den={den}: {n_bad} of {k} numerators disagree; first at "
                f"num={numerators[int(np.argmax(actual != expected))]}"
            )
            n_classes += k

        assert n_classes == 160800, n_classes

    def test_rounding_random_fragments(self, tiny_fasta):
        """Rounding agreement over random fragments drawn from a real contig.

        This is the arm that exercises the real indexing path —
        `cumsum[stop] - cumsum[start]` against a cumsum built by
        get_g_or_c_cumsum — rather than a synthetic cumsum. Exhaustive coverage
        of the quantiser itself is
        test_rounding_exhaustive_over_quantiser_domain's job, which is why 10^5
        draws suffice here (design §3.4).

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

    def test_abort_on_zero_signature_without_sortedness_violation(self):
        """The other single-condition case: cond_b holds but cond_a does not.

        cond_a is `starts[-1] < starts[-2]`, and the phantom row's start is 0,
        so cond_a can only fail when the preceding start is also 0. That is what
        this builds: a real fragment at position 0 followed by an all-zero row.
        The zero signature fires, the sortedness check cannot, and the file must
        go to human review rather than being silently truncated.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test.h5")
            with h5py.File(h5_path, "w") as f:
                f.create_group("data")
                f.create_group("index")
                f.attrs["index_block_size"] = 5000
                f.attrs["max_fragment_length"] = 65535
                f.attrs["_contig_lengths_str"] = str({"chr1": 1000})
                f.attrs["_bam_header"] = ""
                f.attrs["_source_format"] = "BAM"

                grp = f.create_group("data/chr1")
                grp.create_dataset("starts", data=np.array([0, 0], dtype="int32"))
                grp.create_dataset("lengths", data=np.array([30, 0], dtype="uint16"))
                grp.create_dataset("mapq", data=np.array([[60, 60], [0, 0]], dtype="uint8"))
                grp.create_dataset("gc", data=np.array([100, 0], dtype="uint8"))
                grp.create_dataset("strand", data=np.array([b"+", b""], dtype="|S1"))
                f["fragment_length_counts"] = np.zeros(65536)

            with h5py.File(h5_path, "r") as f:
                with pytest.raises(RepairAbort, match="cond_a=False, cond_b=True"):
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

    regions: sequence of 'below' | 'at_T23' | 'middle' | 'at_T24' | 'above' —
        the float32 cumsum band that fragment i's stop coordinate lands in.
        The two 'at_' values sit exactly on a threshold; both 2**23 and 2**24
        are exactly representable in float32, so the comparison really is an
        equality and the choice of < vs <= is observable.
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
        "at_T23": np.float32(T23),
        "middle": np.float32(T23 + 1000.0),
        "at_T24": np.float32(T24),
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

    # -- exact threshold equality (audit finding F4) ------------------------
    #
    # Both bands are half-open: below is `< T23`, middle is `>= T23 & < T24`.
    # A fragment landing exactly on a threshold is the only input that can tell
    # `<` from `<=`, and neither threshold needs a 16 Mb contig to reach here
    # because _region_inputs writes the float32 cumsum directly.

    def test_exactly_t23_is_middle_band_not_below(self):
        """f32_at_stop == T23 belongs to the middle band, not the below band.

        With `below_t23` widened to `<=`, this fragment lands in both bands and
        its change is reported as a pre-saturation violation.
        """
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(
            ["at_T23"], n_flags=[True]
        )
        assert f32_cumsum[int(starts[0]) + FRAG_LEN] == np.float32(T23)

        old_gc = np.array([100], dtype=np.uint8)
        new_gc = np.array([101], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert violations == []
        assert stats["below_T23"] == 0
        assert stats["between_T23_T24"] == 1
        assert stats["changed_below_T23"] == 0
        assert stats["changed_between_T23_T24"] == 1

    def test_exactly_t24_is_above_band_not_middle(self):
        """f32_at_stop == T24 belongs to the above band, not the middle band.

        With the middle band widened to `<= T24`, this N-free fragment's change
        is reported as an illegal middle-band change on an N-free span.
        """
        starts, lengths, f32_cumsum, n_prefix = _region_inputs(
            ["at_T24"], n_flags=[False]
        )
        assert f32_cumsum[int(starts[0]) + FRAG_LEN] == np.float32(T24)

        old_gc = np.array([0], dtype=np.uint8)
        new_gc = np.array([200], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert violations == []
        assert stats["between_T23_T24"] == 0
        assert stats["changed_between_T23_T24"] == 0
        assert stats["above_T24"] == 1
        assert stats["changed_above_T24"] == 1

    # -- T23 half-integer crossing (pilot finding 2026-08-26) ----------------
    #
    # Below T23, float32 represents half-integers exactly. If the accumulator
    # arrives at T23 - 0.5 = 8,388,607.5 and the next base is G/C (+1.0), the
    # result 8,388,608.5 is not representable in [T23, T24) (only integers),
    # so it rounds to 8,388,608.0 — a -0.5 error on a non-N base. This is
    # the one position per contig where the middle-band error can change
    # without an N, and fragments spanning it legitimately change GC.

    def test_middle_band_n_free_change_spanning_t23_crossing_permitted(self):
        """An N-free fragment spanning the T23 half-integer crossing is allowed."""
        crossing_idx = 15
        size = 30
        f32_cumsum = np.zeros(size, dtype=np.float32)
        for i in range(size):
            if i < crossing_idx:
                f32_cumsum[i] = np.float32(T23 - 100 + i)
            else:
                f32_cumsum[i] = np.float32(T23 + (i - crossing_idx))
        # The value just before the crossing is a half-integer
        f32_cumsum[crossing_idx - 1] = np.float32(T23 - 0.5)

        # Fragment spans [10, 20): covers crossing at 15
        starts = np.array([10], dtype=np.int32)
        lengths = np.array([10], dtype=np.uint16)

        # No N in the span
        n_prefix = np.zeros(size, dtype=np.int64)

        old_gc = np.array([100], dtype=np.uint8)
        new_gc = np.array([101], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert violations == [], f"Expected no violations but got: {violations}"
        assert stats["changed_between_T23_T24"] == 1
        assert stats["changed_between_T23_T24_no_N"] == 0

    def test_middle_band_n_free_change_not_spanning_t23_crossing_blocked(self):
        """An N-free middle-band change NOT spanning the T23 crossing is still blocked."""
        crossing_idx = 15
        size = 30
        f32_cumsum = np.zeros(size, dtype=np.float32)
        for i in range(size):
            if i < crossing_idx:
                f32_cumsum[i] = np.float32(T23 - 100 + i)
            else:
                f32_cumsum[i] = np.float32(T23 + (i - crossing_idx))
        f32_cumsum[crossing_idx - 1] = np.float32(T23 - 0.5)

        # Fragment starts PAST the crossing: [20, 25)
        starts = np.array([20], dtype=np.int32)
        lengths = np.array([5], dtype=np.uint16)
        n_prefix = np.zeros(size, dtype=np.int64)

        old_gc = np.array([100], dtype=np.uint8)
        new_gc = np.array([101], dtype=np.uint8)

        stats, violations = classify_gc_changes(
            old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
        )
        assert len(violations) == 1
        assert "N-free" in violations[0]
        assert stats["changed_between_T23_T24_no_N"] == 1


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
# End-to-end integration across the T23 / T24 seam
#
# Everything above tests the three pieces of the saturation analysis in
# isolation: simulate_float32_cumsum against a directly accumulated float32
# array, build_n_prefix_sum against a 200 bp contig, and classify_gc_changes
# against hand-built f32_cumsum / n_prefix inputs. Nothing joined them, and the
# hand-built region inputs structurally cannot: reaching T24 needs ~16.8 Mb of
# accumulated G/C, and that band is where all 218 damaged files actually live.
#
# The fixture below is a real ~17.2 Mb contig written to a real FASTA. The chain
# under test is the production one end to end — get_g_or_c_cumsum ->
# simulate_float32_cumsum -> build_n_prefix_sum -> classify_gc_changes ->
# recompute_gc_for_contig — with no hand-constructed array anywhere in it.
#
# Cost: ~3 s, ~730 MB peak RSS, a 17.5 MB FASTA and a 137 MB cumsum cache in a
# temp dir. Marked slow, like tests/test_gc_cumsum_overflow.py (~676 MB).
# ---------------------------------------------------------------------------

# Per-base contribution as a 256-entry lookup, built from the same independent
# _GC_CONTRIBUTION table the oracle uses. Shares no code with repair.py or with
# fragments_h5.sequence.
_GC_LUT = np.zeros(256, dtype=np.float32)
for _ch, _v in _GC_CONTRIBUTION.items():
    _GC_LUT[ord(_ch)] = _v
    _GC_LUT[ord(_ch.lower())] = _v


def gc_contributions(sequence):
    """Per-base C+G contributions with the builder's leading pad element.

    get_g_or_c_cumsum prepends an 'a' before encoding (fragment.py:441), so the
    accumulated value at contig position i sits at index i. The 0.0 first
    element here reproduces that offset.
    """
    raw = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    out = np.empty(len(raw) + 1, dtype=np.float32)
    out[0] = np.float32(0.0)
    out[1:] = _GC_LUT[raw]
    return out


def naive_float32_gc_uint8(f32_cumsum, starts, lengths):
    """GC bytes exactly as the pre-fix builder produced them.

    fragment.py:496 took the numerator out of the float32 accumulator and
    divided in Python float; fragments_h5.py:845 quantised with
    int(round(gc * 254)). Feeding a float32 cumsum here reproduces the damage.
    """
    out = np.empty(len(starts), dtype=np.uint8)
    for i in range(len(starts)):
        start = int(starts[i])
        length = int(lengths[i])
        num = float(f32_cumsum[start + length] - f32_cumsum[start])
        out[i] = int(round(round(num / float(length), 5) * 254))
    return out


# Layout of the saturating contig. Each named block carries fragments; the
# filler between blocks is pure G, which advances the accumulator 1.0 per base.
SAT_BLOCK_LEN = 2500
SAT_FRAGS_PER_BLOCK = 250
SAT_DRIFT_N = 400_000
SAT_BLOCKS = (
    "below_free", "below_n", "mid_free", "mid_n", "boundary",
    "above_free", "above_n",
)


def _build_saturating_sequence():
    """A contig whose accumulated G/C crosses both 2**23 and 2**24.

    Returns (sequence, {block_name: (position, length)}). The layout is chosen
    so that every branch of the §5b region rule is exercised by real data:

      below_free / below_n  — under T23, where float32 is exact to half-integers
                              so even N-bearing spans must come back unchanged.
      mid_free / mid_n      — inside [T23, T24), where the ulp is 1.0. N-free
                              spans are still exact; N-bearing spans lose the
                              0.5 to round-half-to-even and must change.
      drift                 — a long N run inside the middle band. Adding 0.5 to
                              an even float32 integer is a tie that rounds back
                              to itself, so the float32 accumulator stalls
                              completely while float64 advances 0.5 per base.
                              This opens a ~200k-wide window in which the two
                              accumulators disagree about which band a fragment
                              is in. It carries no fragments; it exists to make
                              that disagreement observable.
      boundary              — placed in that window: float32 says middle band,
                              float64 says above T24. N-free, so its bytes are
                              unchanged either way and only the band counts move.
      above_free / above_n  — past T24, where float32 has flatlined, every
                              numerator is 0 and every stored byte is 0. This is
                              the production corruption signature.
    """
    rng = np.random.RandomState(2024)
    chunks = []
    blocks = {}
    acc = 0.0
    pos = 0

    def filler(target_acc):
        nonlocal acc, pos
        n = int(target_acc - acc)
        assert n > 0, f"filler target {target_acc} is behind accumulator {acc}"
        chunks.append("G" * n)
        acc += n
        pos += n

    def block(name, seq):
        nonlocal acc, pos
        blocks[name] = (pos, len(seq))
        chunks.append(seq)
        acc += float(_GC_LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)].sum())
        pos += len(seq)

    def varied(n, n_frac):
        weights = [(1 - n_frac) * w for w in (0.28, 0.25, 0.25, 0.22)] + [n_frac]
        return "".join(rng.choice(list("ACGTN"), size=n, p=weights))

    filler(T23 - 300_000)
    block("below_free", varied(SAT_BLOCK_LEN, 0.0))
    block("below_n", varied(SAT_BLOCK_LEN, 0.10))

    filler(T23 + 400_000)
    block("mid_free", varied(SAT_BLOCK_LEN, 0.0))
    block("mid_n", varied(SAT_BLOCK_LEN, 0.10))

    filler(T24 - 100_000)
    block("drift", "N" * SAT_DRIFT_N)
    block("boundary", varied(SAT_BLOCK_LEN, 0.0))

    # Push float32 the rest of the way to T24 and past the flatline onset.
    chunks.append("G" * 150_000)
    acc += 150_000
    pos += 150_000
    block("above_free", varied(SAT_BLOCK_LEN, 0.0))
    block("above_n", varied(SAT_BLOCK_LEN, 0.10))

    # Trailing margin so no fragment can run off the end of the contig.
    chunks.append("G" * 2000)
    acc += 2000
    pos += 2000
    return "".join(chunks), blocks


def _place_saturating_fragments(blocks):
    """Sorted starts / lengths, SAT_FRAGS_PER_BLOCK inside each named block."""
    rng = np.random.RandomState(99)
    all_starts, all_lengths = [], []
    for name in SAT_BLOCKS:
        pos, width = blocks[name]
        lengths = rng.randint(10, 250, SAT_FRAGS_PER_BLOCK).astype(np.int64)
        starts = pos + rng.randint(0, width - 250, SAT_FRAGS_PER_BLOCK).astype(np.int64)
        order = np.argsort(starts)
        all_starts.append(starts[order])
        all_lengths.append(lengths[order])
    starts = np.concatenate(all_starts).astype("int32")
    lengths = np.concatenate(all_lengths).astype("uint16")
    assert np.all(np.diff(starts) >= 0), "block layout must keep starts sorted"
    return starts, lengths


@pytest.fixture(scope="module")
def saturating_fixture():
    """A damaged H5 over a contig that really crosses T23 and T24.

    Yields paths plus the independently derived expectations. The big arrays are
    dropped before yielding so the per-test peak is dominated by repair itself.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        sequence, blocks = _build_saturating_sequence()
        starts, lengths = _place_saturating_fragments(blocks)
        stops = starts.astype(np.int64) + lengths.astype(np.int64)

        contributions = gc_contributions(sequence)
        f32_cumsum = np.cumsum(contributions)
        f64_cumsum = np.cumsum(contributions.astype(np.float64))

        damaged_gc = naive_float32_gc_uint8(f32_cumsum, starts, lengths)
        oracle_gc = naive_gc_uint8(sequence, starts, lengths)

        at_stop = f32_cumsum[stops]
        below = at_stop < T23
        middle = (at_stop >= T23) & (at_stop < T24)
        above = at_stop >= T24

        f64_at_stop = f64_cumsum[stops]
        f64_bands = (
            int(np.sum(f64_at_stop < T23)),
            int(np.sum((f64_at_stop >= T23) & (f64_at_stop < T24))),
            int(np.sum(f64_at_stop >= T24)),
        )

        raw = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        is_n = raw == ord("N")
        n_prefix = np.zeros(len(raw) + 1, dtype=np.int64)
        n_prefix[1:] = np.cumsum(is_n)
        has_n = (n_prefix[stops] - n_prefix[starts]) > 0

        changed = damaged_gc != oracle_gc
        expected = {
            "total_frags": len(starts),
            "below_T23": int(np.sum(below)),
            "between_T23_T24": int(np.sum(middle)),
            "above_T24": int(np.sum(above)),
            "changed_total": int(np.sum(changed)),
            "changed_below_T23": int(np.sum(changed & below)),
            "changed_between_T23_T24": int(np.sum(changed & middle)),
            "changed_between_T23_T24_no_N": int(np.sum(changed & middle & ~has_n)),
            "changed_above_T24": int(np.sum(changed & above)),
        }

        info = {
            "starts": starts,
            "lengths": lengths,
            "damaged_gc": damaged_gc,
            "oracle_gc": oracle_gc,
            "expected_stats": expected,
            "f32_bands": (expected["below_T23"], expected["between_T23_T24"],
                          expected["above_T24"]),
            "f64_bands": f64_bands,
            "f32_max": float(f32_cumsum.max()),
            "n_middle_n_free": int(np.sum(middle & ~has_n)),
            "n_below_with_n": int(np.sum(below & has_n)),
            "n_below_with_n_changed": int(np.sum(below & has_n & changed)),
            "n_rounding_sensitive": count_rounding_sensitive(sequence, starts, lengths),
        }

        fasta_path = os.path.join(tmpdir, "saturating.fa")
        _write_fasta(fasta_path, [("chrT", sequence)])

        h5_path = os.path.join(tmpdir, "damaged.h5")
        _make_fragment_h5(
            h5_path,
            {"chrT": {"starts": starts, "lengths": lengths, "gc": damaged_gc}},
            contig_lengths={"chrT": len(sequence)},
        )

        del sequence, contributions, f32_cumsum, f64_cumsum, raw, is_n, n_prefix

        cache_dir = os.path.join(tmpdir, "cumsum_cache")
        fasta_sha256 = compute_file_sha256(fasta_path)
        from fragments_h5.repair import build_cumsum_cache
        build_cumsum_cache(fasta_path, cache_dir, fasta_sha256)

        info.update(
            fasta=fasta_path, h5=h5_path,
            cache_dir=cache_dir, fasta_sha256=fasta_sha256,
        )
        yield info


@pytest.mark.slow
class TestSaturatingContigEndToEnd:
    """The T23/T24 seam, joined: real cumsum -> real simulation -> real rules."""

    def test_fixture_actually_crosses_both_thresholds(self, saturating_fixture):
        """Guard on the fixture: without these properties the tests below are vacuous."""
        info = saturating_fixture
        below, middle, above = info["f32_bands"]

        # The float32 accumulator really flatlines at exactly 2**24.
        assert info["f32_max"] == float(T24), info["f32_max"]

        # All three bands carry real fragments.
        assert below >= SAT_FRAGS_PER_BLOCK, below
        assert middle >= SAT_FRAGS_PER_BLOCK, middle
        assert above >= SAT_FRAGS_PER_BLOCK, above

        stats = info["expected_stats"]
        # Below T23 the ulp is <= 0.5, so half-integer N contributions survive and
        # nothing may change — including on spans that do contain N.
        assert info["n_below_with_n"] > 0
        assert info["n_below_with_n_changed"] == 0
        assert stats["changed_below_T23"] == 0

        # The middle band contains both kinds of span, and only the N-bearing
        # ones lose information.
        assert info["n_middle_n_free"] > 0
        assert stats["changed_between_T23_T24"] > 0
        assert stats["changed_between_T23_T24_no_N"] == 0

        # Above T24 the damage is total: every stored byte collapsed to 0.
        assert stats["changed_above_T24"] == above
        assert np.array_equal(
            np.unique(info["damaged_gc"][-2 * SAT_FRAGS_PER_BLOCK:]),
            np.array([0], dtype=np.uint8),
        )

        # The fixture can tell round(x, 5) from round(x, 4), same bar as the
        # small e2e fixtures (test_fixtures_are_rounding_sensitive).
        assert info["n_rounding_sensitive"] >= 3, info["n_rounding_sensitive"]

        # And it can tell a float32 accumulator from a float64 one: the drift
        # block moves a whole block of fragments across the T24 boundary.
        assert info["f32_bands"] != info["f64_bands"], (
            f"float32 and float64 agree on every band ({info['f32_bands']}); "
            f"the fixture cannot detect a simulation that skips the float32 cast"
        )

    def test_dry_run_band_classification_matches_independent_accumulator(
        self, saturating_fixture,
    ):
        """Every §5b statistic matches one computed from an independent float32 cumsum.

        This is the seam. The numbers on the left come out of the real chain
        (cached float64 cumsum -> simulate_float32_cumsum -> build_n_prefix_sum
        -> classify_gc_changes); the ones on the right come from accumulating
        the FASTA text in float32 with a lookup table written in this file.
        """
        info = saturating_fixture
        report = repair_local_file(
            local_path=info["h5"],
            fasta_path=info["fasta"],
            cumsum_cache_dir=info["cache_dir"],
            fasta_sha256=info["fasta_sha256"],
            dry_run=True,
        )
        assert report["status"] == "dry_run_ok"

        stats = report["gc_results"]["chrT"]["stats"]
        assert report["gc_results"]["chrT"]["violations"] == []
        for key, want in info["expected_stats"].items():
            assert stats[key] == want, (
                f"{key}: chain says {stats[key]}, independent float32 "
                f"accumulator says {want}"
            )

    def test_apply_restores_the_oracle_bytes(self, saturating_fixture):
        """Repairing the damaged file reproduces the naive character-count oracle."""
        import shutil

        info = saturating_fixture
        with tempfile.TemporaryDirectory() as tmpdir:
            target = os.path.join(tmpdir, "apply.h5")
            shutil.copy2(info["h5"], target)

            report = repair_local_file(
                local_path=target,
                fasta_path=info["fasta"],
                cumsum_cache_dir=info["cache_dir"],
                fasta_sha256=info["fasta_sha256"],
                dry_run=False,
            )
            assert report["status"] == "ok"
            assert report["contigs_truncated"] == 0

            with h5py.File(target, "r") as f:
                repaired = f["data/chrT/gc"][:]
                assert np.array_equal(f["data/chrT/starts"][:], info["starts"])
                assert np.array_equal(f["data/chrT/lengths"][:], info["lengths"])
                assert len(json.loads(f.attrs["_repair_history"])) == 1

            n_wrong = int(np.sum(repaired != info["oracle_gc"]))
            assert n_wrong == 0, (
                f"{n_wrong} of {len(repaired)} repaired bytes disagree with the "
                f"independent oracle; first at index "
                f"{int(np.argmax(repaired != info['oracle_gc']))}"
            )
            # And the repair really did something: the damaged input differed.
            assert not np.array_equal(info["damaged_gc"], info["oracle_gc"])


# ---------------------------------------------------------------------------
# T23 half-integer crossing end-to-end
#
# Below T23, float32 represents half-integers exactly, so the accumulator can
# arrive at 8,388,607.5. Adding 1.0 for a G/C gives 8,388,608.5, not
# representable in [T23, T24) (only integers), so it rounds to 8,388,608.0.
# This is the one non-N error-change position per contig.
#
# The fixture deterministically produces the crossing: 8,388,607 G's + 1 N +
# trailing G's. The float64 cumsum is built directly (not from a FASTA) so the
# test avoids writing an 8.4 MB file, but the float32 accumulation is real —
# simulate_float32_cumsum runs actual numpy float32 cumsum over the diffs.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def t23_crossing_fixture():
    """Float32 accumulator crossing T23 from a half-integer, exercised for real.

    Sequence layout: 8,388,607 G's + 1 N + 2000 G's.
    Cumsum: [..., 8388607.0, 8388607.5, 8388608.0(f32)/8388608.5(f64), ...]
    Crossing at cumsum index 8,388,609.
    """
    n_g_before = 8_388_607
    n_trailing = 2000
    cumsum_size = n_g_before + 1 + n_trailing + 1  # bases + leading pad

    # Build float64 cumsum: all-G contributes 1.0 each, one N contributes 0.5
    float64_cumsum = np.arange(cumsum_size, dtype=np.float64)
    float64_cumsum[n_g_before + 1:] -= 0.5  # N contributed 0.5, not 1.0

    # Real float32 simulation
    f32_cumsum = simulate_float32_cumsum(float64_cumsum)
    n_prefix = build_n_prefix_sum(float64_cumsum)

    crossing_idx = n_g_before + 2  # 8,388,609

    # Sequence for the oracle
    sequence = "G" * n_g_before + "N" + "G" * n_trailing

    # Fragments: (D) below T23, (C) spans N + crossing, (A) spans crossing no N,
    #            (B) past crossing no N
    starts = np.array([
        1000,                    # D: well below T23
        crossing_idx - 2,        # C: spans N + crossing
        crossing_idx - 1,        # A: spans crossing, no N
        crossing_idx,            # B: just past crossing
    ], dtype=np.int32)
    lengths = np.array([100, 3, 2, 100], dtype=np.uint16)

    damaged_gc = naive_float32_gc_uint8(f32_cumsum, starts, lengths)
    oracle_gc = naive_gc_uint8(sequence, starts, lengths)

    yield {
        "f32_cumsum": f32_cumsum,
        "float64_cumsum": float64_cumsum,
        "n_prefix": n_prefix,
        "starts": starts,
        "lengths": lengths,
        "damaged_gc": damaged_gc,
        "oracle_gc": oracle_gc,
        "crossing_idx": crossing_idx,
        "sequence": sequence,
    }


@pytest.mark.slow
class TestT23CrossingEndToEnd:
    """The T23 half-integer crossing, exercised with real float32 arithmetic."""

    def test_crossing_position_verified(self, t23_crossing_fixture):
        """Guard: the fixture really has a half-integer -> integer transition."""
        info = t23_crossing_fixture
        c = info["crossing_idx"]
        f32 = info["f32_cumsum"]

        # Value before the crossing is exactly T23 - 0.5
        assert float(f32[c - 1]) == float(np.float32(T23 - 0.5))
        # Value at the crossing is exactly T23
        assert float(f32[c]) == float(np.float32(T23))
        # find_t23_crossing returns the correct index
        assert find_t23_crossing(f32, info["n_prefix"]) == c

    def test_damaged_gc_differs_at_crossing_fragment(self, t23_crossing_fixture):
        """The float32 bug causes a GC change at the crossing fragment."""
        info = t23_crossing_fixture
        # Fragment A (index 2): spans crossing, no N — GC should differ
        assert info["damaged_gc"][2] != info["oracle_gc"][2], (
            f"damaged={info['damaged_gc'][2]}, oracle={info['oracle_gc'][2]}"
        )
        # Fragment B (index 3): past crossing, no N — GC should agree
        assert info["damaged_gc"][3] == info["oracle_gc"][3], (
            f"damaged={info['damaged_gc'][3]}, oracle={info['oracle_gc'][3]}"
        )
        # Fragment D (index 0): below T23 — GC should agree
        assert info["damaged_gc"][0] == info["oracle_gc"][0]

    def test_classify_allows_crossing_span(self, t23_crossing_fixture):
        """classify_gc_changes permits the N-free crossing change."""
        info = t23_crossing_fixture
        stats, violations = classify_gc_changes(
            info["damaged_gc"], info["oracle_gc"],
            info["starts"], info["lengths"],
            info["f32_cumsum"], info["n_prefix"],
        )
        assert violations == [], f"Unexpected violations: {violations}"
        assert stats["changed_between_T23_T24_t23_crossing"] >= 1
        assert stats["changed_between_T23_T24_no_N"] == 0
        assert stats["changed_below_T23"] == 0


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
# Index eligibility exactly AT the thresholds
#
# Both index guards are inequalities against a named constant, and the
# interesting behaviour of an inequality lives on the boundary, not near it.
# test_index_skip_few_reads uses 99 reads and test_index_skip_short_contig uses
# 4999 bp — each one step outside the threshold, so neither pins down which side
# of the boundary the constant sits on. A production file (RD-56670-Lib1,
# chrUn_KI270507v1) had exactly MIN_NUM_READS_FOR_INDEX rows *including the
# phantom padding row*; truncation took it to MIN-1, it correctly left the index,
# and the repair tool aborted mid-write because its index key-set invariant did
# not allow for that. These tests sit exactly on both boundaries.
# ---------------------------------------------------------------------------

class TestIndexThresholdBoundaries:

    @pytest.mark.parametrize("n_frags, expected", [
        (MIN_NUM_READS_FOR_INDEX - 1, False),
        (MIN_NUM_READS_FOR_INDEX, True),        # the boundary itself
        (MIN_NUM_READS_FOR_INDEX + 1, True),
    ])
    def test_read_count_guard_is_inclusive_at_the_threshold(self, n_frags, expected):
        """n >= MIN_NUM_READS_FOR_INDEX is indexed; the threshold value qualifies."""
        assert should_index_contig(n_frags, 100000, INDEX_BLOCK_SIZE) is expected

    @pytest.mark.parametrize("contig_length, expected", [
        (INDEX_BLOCK_SIZE - 1, False),
        (INDEX_BLOCK_SIZE, False),              # the boundary itself
        (INDEX_BLOCK_SIZE + 1, True),
    ])
    def test_length_guard_is_exclusive_at_the_threshold(self, contig_length, expected):
        """contig_length > index_block_size is indexed; the threshold value does not."""
        assert should_index_contig(200, contig_length, INDEX_BLOCK_SIZE) is expected

    def test_rebuild_at_exactly_the_read_threshold(self):
        """rebuild_contig_index agrees with the predicate on both sides of MIN."""
        with tempfile.TemporaryDirectory() as tmpdir:
            for n, expect_index in (
                (MIN_NUM_READS_FOR_INDEX, True),
                (MIN_NUM_READS_FOR_INDEX - 1, False),
            ):
                h5_path = os.path.join(tmpdir, f"n{n}.h5")
                _make_fragment_h5(h5_path, {
                    "chr1": {
                        "starts": np.arange(0, 10 * n, 10, dtype="int32"),
                        "lengths": np.full(n, 5, dtype="uint16"),
                    },
                }, contig_lengths={"chr1": 100000})

                with h5py.File(h5_path, "r+") as f:
                    rebuilt = rebuild_contig_index(f, "chr1", 100000, INDEX_BLOCK_SIZE)
                    assert rebuilt is expect_index, f"n={n}"
                    assert ("chr1" in f["index"]) is expect_index, f"n={n}"

    def test_predict_matches_rebuild_at_the_boundary(self):
        """The pre-write forecast and the actual rebuild cannot disagree.

        The forecast is what makes the key-set check possible before any
        mutation, so it must track rebuild_contig_index exactly — including at
        the threshold, where the two used to be written out separately.
        """
        cases = {
            "at_threshold": MIN_NUM_READS_FOR_INDEX,
            "below_threshold": MIN_NUM_READS_FOR_INDEX - 1,
        }
        contig_lengths = {c: 100000 for c in cases}

        predicted = predict_index_contigs(cases, contig_lengths, INDEX_BLOCK_SIZE)

        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "both.h5")
            _make_fragment_h5(h5_path, {
                c: {
                    "starts": np.arange(0, 10 * n, 10, dtype="int32"),
                    "lengths": np.full(n, 5, dtype="uint16"),
                }
                for c, n in cases.items()
            }, contig_lengths=contig_lengths)

            actual = set()
            with h5py.File(h5_path, "r+") as f:
                for c in cases:
                    if rebuild_contig_index(f, c, 100000, INDEX_BLOCK_SIZE):
                        actual.add(c)

        assert predicted == actual == {"at_threshold"}


# ---------------------------------------------------------------------------
# A contig sitting exactly ON the read threshold, end-to-end (§3.2.3)
# ---------------------------------------------------------------------------

# Absent from the FASTA (so gc is all-255 and needs no reference), but long
# enough to clear the length guard, and carrying exactly
# MIN_NUM_READS_FOR_INDEX rows once the phantom row is included. This is the
# shape of chrUn_KI270507v1 in RD-56670-Lib1.
THRESHOLD_NAME = "chrUn_threshold99"
THRESHOLD_LEN = 6000
THRESHOLD_N_REAL = MIN_NUM_READS_FOR_INDEX - 1


def _threshold_contig():
    starts = np.sort(
        np.random.RandomState(310).randint(1, THRESHOLD_LEN - 300, THRESHOLD_N_REAL)
    ).astype("int32")
    lengths = np.random.RandomState(311).randint(
        10, 200, THRESHOLD_N_REAL
    ).astype("uint16")
    return {
        "starts": starts,
        "lengths": lengths,
        "gc": np.full(THRESHOLD_N_REAL, 255, dtype="uint8"),
    }


class TestContigLeavesIndexAtThreshold:
    """A contig that drops below MIN_NUM_READS_FOR_INDEX purely because
    truncation removed its phantom row must not fail the repair."""

    @pytest.fixture
    def threshold_h5(self, tiny_fasta):
        with tempfile.TemporaryDirectory() as tmpdir:
            contigs = _build_fixture_contigs(tiny_fasta, (50, 51, 52, 53))
            contigs[THRESHOLD_NAME] = _threshold_contig()
            path = os.path.join(tmpdir, "threshold.h5")
            _make_fragment_h5(
                path, contigs, add_padding=True,
                contig_lengths=dict(
                    FIXTURE_CONTIG_LENGTHS, **{THRESHOLD_NAME: THRESHOLD_LEN}
                ),
            )
            yield path

    def test_fixture_sits_exactly_on_the_threshold(self, threshold_h5):
        """Guard the fixture itself: off by one either way and the test is vacuous."""
        with h5py.File(threshold_h5, "r") as f:
            n_with_phantom = len(f[f"data/{THRESHOLD_NAME}/starts"])
            assert n_with_phantom == MIN_NUM_READS_FOR_INDEX
            assert THRESHOLD_NAME in f["index"], (
                "contig must start out indexed, or there is nothing to lose"
            )
            assert f.attrs["index_block_size"] < THRESHOLD_LEN, (
                "contig must clear the length guard, or the drop has two causes"
            )

    def test_apply_succeeds_and_contig_leaves_the_index(
        self, threshold_h5, tiny_fasta, cumsum_cache
    ):
        """The whole repair completes; the contig is legitimately unindexed."""
        cache_dir, fasta_sha256 = cumsum_cache

        report = repair_local_file(
            local_path=threshold_h5,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=False,
        )

        assert report["status"] == "ok"
        assert report["index_contigs_dropped"] == [THRESHOLD_NAME]

        with h5py.File(threshold_h5, "r") as f:
            assert len(f[f"data/{THRESHOLD_NAME}/starts"]) == THRESHOLD_N_REAL
            assert THRESHOLD_NAME not in f["index"]
            # every other contig kept its index
            assert set(f["index"].keys()) == {"chr1"}

    def test_dry_run_reports_the_drop(self, threshold_h5, tiny_fasta, cumsum_cache):
        """The drop is visible before --apply, because it is decided pre-write."""
        cache_dir, fasta_sha256 = cumsum_cache

        report = repair_local_file(
            local_path=threshold_h5,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=True,
        )

        assert report["status"] == "dry_run_ok"
        assert report["index_contigs_dropped"] == [THRESHOLD_NAME]


class TestIndexKeySetStillGuarded:
    """The threshold exemption must not blunt the invariant."""

    @pytest.fixture
    def bogus_index_h5(self, tiny_fasta):
        """A padded file carrying an index on a contig too short to have one.

        Nothing about truncation explains losing this index, so the rebuild
        reshaping the key set is exactly the corruption the invariant exists to
        catch.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            contigs = _build_fixture_contigs(tiny_fasta, (50, 51, 52, 53))
            path = os.path.join(tmpdir, "bogus.h5")
            _make_fragment_h5(
                path, contigs, add_padding=True,
                contig_lengths=FIXTURE_CONTIG_LENGTHS,
            )
            with h5py.File(path, "r+") as f:
                assert FIXTURE_CONTIG_LENGTHS["chr2"] <= f.attrs["index_block_size"]
                assert "chr2" not in f["index"]
                f["index/chr2"] = np.array([0, len(f["data/chr2/starts"])])
            yield path

    def test_unexplained_index_loss_still_aborts(
        self, bogus_index_h5, tiny_fasta, cumsum_cache
    ):
        cache_dir, fasta_sha256 = cumsum_cache

        with pytest.raises(RepairAbort, match="Index key set changed"):
            repair_local_file(
                local_path=bogus_index_h5,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=False,
            )

    def test_unexplained_index_loss_aborts_before_writing(
        self, bogus_index_h5, tiny_fasta, cumsum_cache
    ):
        """The abort leaves the file byte-identical — no partial repair.

        On RD-56670-Lib1 this check fired *after* the datasets had been
        rewritten, growing the file from 64,941,883 to 74,417,833 bytes and
        leaving a half-repaired file behind. Validating the key set from
        in-memory state removes that window.
        """
        cache_dir, fasta_sha256 = cumsum_cache
        before = compute_file_sha256(bogus_index_h5)
        size_before = os.path.getsize(bogus_index_h5)

        with pytest.raises(RepairAbort):
            repair_local_file(
                local_path=bogus_index_h5,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=False,
            )

        assert compute_file_sha256(bogus_index_h5) == before
        assert os.path.getsize(bogus_index_h5) == size_before


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


# ---------------------------------------------------------------------------
# Contigs present in the H5 but absent from the FASTA (audit finding F2)
#
# This branch runs on ~165 scaffold contigs in every one of the 218 production
# files, and it is a §5 layer-3 wrong-reference defence: if the operator hands
# the tool a reference that is merely *different* rather than obviously wrong,
# the contigs it does not contain are exactly where that shows up. Both
# directions matter — skipping the genuinely unreferenced ones, and refusing
# to skip one that carries real GC data.
# ---------------------------------------------------------------------------

class TestContigAbsentFromFasta:
    """The `contig not in fasta_ref.references` branch, both directions."""

    def test_absent_contig_with_all_255_gc_is_skipped(
        self, absent_contig_h5, tiny_fasta, cumsum_cache,
    ):
        """gc == 255 everywhere means 'no reference for this contig' — skip it."""
        cache_dir, fasta_sha256 = cumsum_cache
        h5_path = absent_contig_h5()

        report = repair_local_file(
            local_path=h5_path,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=True,
        )

        assert report["status"] == "dry_run_ok"
        res = report["gc_results"][SCAFFOLD_NAME]
        assert res["action"] == "skipped_absent"
        assert res["violations"] == []
        assert res["stats"] == {"total_frags": SCAFFOLD_N_FRAGS, "all_255": True}

        # The FASTA-backed contigs still went down the recompute path, so the
        # skip is contig-scoped rather than a whole-file bail-out.
        for contig in ("chr1", "chr2"):
            assert report["gc_results"][contig]["action"] == "recomputed"

    def test_absent_contig_gc_untouched_by_apply(
        self, absent_contig_h5, tiny_fasta, cumsum_cache, expected_clean_gc,
    ):
        """--apply must not rewrite the skipped contig's gc, and must succeed."""
        cache_dir, fasta_sha256 = cumsum_cache
        h5_path = absent_contig_h5()

        with h5py.File(h5_path, "r") as f:
            before = f[f"data/{SCAFFOLD_NAME}/gc"][:]

        report = repair_local_file(
            local_path=h5_path,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=False,
        )
        assert report["status"] == "ok"

        with h5py.File(h5_path, "r") as f:
            after = f[f"data/{SCAFFOLD_NAME}/gc"][:]
            assert np.array_equal(after, before)
            assert np.all(after == 255)
            # ... while the referenced contigs still match the oracle
            for contig in ("chr1", "chr2"):
                assert np.array_equal(
                    f[f"data/{contig}/gc"][:], expected_clean_gc[contig]
                )

    def test_absent_contig_with_non_255_gc_aborts(
        self, absent_contig_h5, tiny_fasta, cumsum_cache,
    ):
        """Real GC bytes on a contig the FASTA lacks means the wrong reference."""
        cache_dir, fasta_sha256 = cumsum_cache
        h5_path = absent_contig_h5(n_non_255=3)

        with pytest.raises(RepairAbort, match="absent from FASTA but gc has 3 non-255"):
            repair_local_file(
                local_path=h5_path,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=True,
            )

    def test_absent_contig_abort_writes_nothing(
        self, absent_contig_h5, tiny_fasta, cumsum_cache,
    ):
        """The abort fires during analysis, before --apply touches the file."""
        cache_dir, fasta_sha256 = cumsum_cache
        h5_path = absent_contig_h5(n_non_255=1)

        sha_before = compute_file_sha256(h5_path)
        with pytest.raises(RepairAbort, match="absent from FASTA"):
            repair_local_file(
                local_path=h5_path,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=False,
            )
        assert compute_file_sha256(h5_path) == sha_before


# ---------------------------------------------------------------------------
# index_block_size comes from the file, not the module (audit finding F1)
# ---------------------------------------------------------------------------

class TestIndexBlockSizeFromFile:
    """§3.2.3: the index must be rebuilt with the block size the file was built with."""

    def test_rebuild_honours_non_default_block_size(self, tiny_fasta, cumsum_cache):
        """A file written with a non-default block size keeps that geometry.

        INDEX_BLOCK_SIZE has changed over the format's life, so a file's own
        `index_block_size` attr is the only authority. Rebuilding with the
        module constant instead would silently reshape the index of every file
        built under an older value.
        """
        cache_dir, fasta_sha256 = cumsum_cache
        block_size = 10000
        assert block_size != INDEX_BLOCK_SIZE, "fixture must not use the default"
        assert CHR1_LEN > block_size and CHR1_LEN > INDEX_BLOCK_SIZE, (
            "chr1 must get an index under both block sizes, so the test "
            "discriminates on the index contents rather than its existence"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "block_size.h5")
            _make_fragment_h5(
                h5_path,
                _build_fixture_contigs(tiny_fasta, (42, 43, 44, 45)),
                index_block_size=block_size,
                contig_lengths=FIXTURE_CONTIG_LENGTHS,
            )

            with h5py.File(h5_path, "r") as f:
                assert int(f.attrs["index_block_size"]) == block_size
                before = f["index/chr1"][:]

            report = repair_local_file(
                local_path=h5_path,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=False,
            )
            assert report["status"] == "ok"

            with h5py.File(h5_path, "r") as f:
                after = f["index/chr1"][:]
                starts = f["data/chr1/starts"][:]

            expected = np.append(
                np.searchsorted(
                    starts, np.arange(0, CHR1_LEN, block_size), side="left"
                ),
                len(starts),
            )
            assert np.array_equal(after, expected), (
                f"rebuilt index {after} != block-size-{block_size} index {expected}"
            )
            assert np.array_equal(after, before), "rebuild changed a correct index"

            # The module default would produce a differently shaped array, which
            # is what makes this test able to see the difference at all.
            n_default = len(np.arange(0, CHR1_LEN, INDEX_BLOCK_SIZE)) + 1
            assert len(after) != n_default, (
                f"both block sizes give {n_default} entries; fixture cannot discriminate"
            )


# ---------------------------------------------------------------------------
# Reference geometry gates (audit finding F3)
#
# Both are §5 layer-3 defences against being handed the wrong reference: a
# FASTA whose contig names match but whose coordinates do not.
# ---------------------------------------------------------------------------

class TestReferenceGeometryGates:
    """Contig-length agreement and in-bounds fragments are checked before any GC."""

    def test_contig_length_mismatch_aborts(self, tiny_fasta, cumsum_cache):
        """FASTA contig length != H5 contig_lengths — the references disagree."""
        cache_dir, fasta_sha256 = cumsum_cache

        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "wrong_length.h5")
            _make_fragment_h5(
                h5_path,
                _build_fixture_contigs(tiny_fasta, (42, 43, 44, 45)),
                # Every fragment is still in bounds and every gc byte is still
                # correct; only the recorded contig length is off by one, so
                # this gate is the only thing that can catch it.
                contig_lengths=dict(FIXTURE_CONTIG_LENGTHS, chr1=CHR1_LEN + 1),
            )

            with pytest.raises(RepairAbort, match="FASTA length 20000 != H5"):
                repair_local_file(
                    local_path=h5_path,
                    fasta_path=tiny_fasta,
                    cumsum_cache_dir=cache_dir,
                    fasta_sha256=fasta_sha256,
                    dry_run=True,
                )

    def test_fragment_past_contig_end_aborts(self, tiny_fasta, cumsum_cache):
        """A fragment whose stop runs past the FASTA contig — abort, do not index.

        Without this gate the overrunning stop is used as a cumsum index and the
        failure mode is a raw IndexError from deep inside the recompute rather
        than a RepairAbort the runner can classify.
        """
        cache_dir, fasta_sha256 = cumsum_cache

        contigs = _build_fixture_contigs(tiny_fasta, (42, 43, 44, 45))
        chr1 = contigs["chr1"]
        # Appended last and started past every existing start, so `starts`
        # stays non-decreasing and the padding-row detector still says 'clean'.
        chr1["starts"] = np.append(chr1["starts"], np.int32(CHR1_LEN - 100))
        chr1["lengths"] = np.append(chr1["lengths"], np.uint16(150))
        chr1["gc"] = np.append(chr1["gc"], np.uint8(100))
        assert np.all(np.diff(chr1["starts"]) >= 0)

        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "overrun.h5")
            _make_fragment_h5(
                h5_path, contigs, contig_lengths=FIXTURE_CONTIG_LENGTHS,
            )

            with pytest.raises(RepairAbort, match="extends past FASTA contig"):
                repair_local_file(
                    local_path=h5_path,
                    fasta_path=tiny_fasta,
                    cumsum_cache_dir=cache_dir,
                    fasta_sha256=fasta_sha256,
                    dry_run=True,
                )


# ---------------------------------------------------------------------------
# Mixed truncation verdicts across contigs
# ---------------------------------------------------------------------------

class TestMixedTruncationVerdicts:
    """A file where some contigs need truncating and others do not is blocking.

    The padding row was appended by a single builder pass, so it is present on
    every contig or none. A file that disagrees with itself is outside the
    understood failure mode and must go to human review.
    """

    def test_mixed_verdicts_abort(self, padded_h5, tiny_fasta, cumsum_cache):
        cache_dir, fasta_sha256 = cumsum_cache

        # Drop the phantom row from chr2 only, leaving chr1 padded.
        with h5py.File(padded_h5, "r+") as f:
            truncate_contig_datasets(f, "chr2")
            assert detect_padding_row(f, "chr1") == "truncate"
            assert detect_padding_row(f, "chr2") == "clean"

        sha_before = compute_file_sha256(padded_h5)
        with pytest.raises(RepairAbort, match="Mixed truncation verdicts"):
            repair_local_file(
                local_path=padded_h5,
                fasta_path=tiny_fasta,
                cumsum_cache_dir=cache_dir,
                fasta_sha256=fasta_sha256,
                dry_run=False,
            )
        assert compute_file_sha256(padded_h5) == sha_before


# ---------------------------------------------------------------------------
# FASTA handle lifetime (implementation review r2)
# ---------------------------------------------------------------------------

class TestFastaHandleLifetime:
    """The FASTA handle is released on every exit path, not just the happy one."""

    @staticmethod
    def _spy_on_fasta(monkeypatch):
        """Record every pysam.FastaFile repair.py opens."""
        import fragments_h5.repair as repair_mod

        opened = []
        real = repair_mod.pysam.FastaFile

        def spy(*args, **kwargs):
            handle = real(*args, **kwargs)
            opened.append(handle)
            return handle

        monkeypatch.setattr(repair_mod.pysam, "FastaFile", spy)
        return opened

    def test_handle_closed_on_success(
        self, clean_h5, tiny_fasta, cumsum_cache, monkeypatch,
    ):
        cache_dir, fasta_sha256 = cumsum_cache
        opened = self._spy_on_fasta(monkeypatch)

        repair_local_file(
            local_path=clean_h5,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=True,
        )

        assert opened, "repair_local_file never opened the FASTA"
        assert all(h.closed for h in opened)

    def test_handle_closed_on_abort(
        self, tiny_fasta, cumsum_cache, monkeypatch,
    ):
        """A RepairAbort raised between the open and the close must not leak it."""
        cache_dir, fasta_sha256 = cumsum_cache
        opened = self._spy_on_fasta(monkeypatch)

        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "wrong_length.h5")
            _make_fragment_h5(
                h5_path,
                _build_fixture_contigs(tiny_fasta, (42, 43, 44, 45)),
                contig_lengths=dict(FIXTURE_CONTIG_LENGTHS, chr1=CHR1_LEN + 1),
            )

            with pytest.raises(RepairAbort):
                repair_local_file(
                    local_path=h5_path,
                    fasta_path=tiny_fasta,
                    cumsum_cache_dir=cache_dir,
                    fasta_sha256=fasta_sha256,
                    dry_run=True,
                )

        assert opened, "the abort happened before the FASTA was opened"
        assert all(h.closed for h in opened), "FASTA handle leaked on abort"
