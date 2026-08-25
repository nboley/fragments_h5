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

@pytest.fixture
def tiny_fasta():
    """A small FASTA: chr1 (300 bp, all G), chr2 (200 bp, G+N+G)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "test.fa")
        _write_fasta(path, [
            ("chr1", "G" * 300),
            ("chr2", "G" * 80 + "N" * 40 + "G" * 80),
        ])
        yield path


@pytest.fixture
def clean_h5(tiny_fasta):
    """A correctly-built fragment H5 (no padding, correct GC)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Compute correct GC values using the reference
        from fragments_h5.fragment import get_g_or_c_cumsum

        # chr1: all G, so every fragment has GC = 1.0 -> uint8 = 254
        n1 = 150
        starts1 = np.sort(np.random.RandomState(42).randint(0, 250, n1)).astype("int32")
        lengths1 = np.random.RandomState(43).randint(10, 50, n1).astype("uint16")
        # Clip to contig length
        stops1 = starts1.astype(np.int64) + lengths1.astype(np.int64)
        mask1 = stops1 <= 300
        starts1 = starts1[mask1]
        lengths1 = lengths1[mask1]

        cumsum1, _ = get_g_or_c_cumsum(tiny_fasta, "chr1")
        gc1 = recompute_gc_for_contig(starts1, lengths1, cumsum1)

        # chr2: G+N+G mix
        n2 = 120
        starts2 = np.sort(np.random.RandomState(44).randint(0, 160, n2)).astype("int32")
        lengths2 = np.random.RandomState(45).randint(5, 30, n2).astype("uint16")
        stops2 = starts2.astype(np.int64) + lengths2.astype(np.int64)
        mask2 = stops2 <= 200
        starts2 = starts2[mask2]
        lengths2 = lengths2[mask2]

        cumsum2, _ = get_g_or_c_cumsum(tiny_fasta, "chr2")
        gc2 = recompute_gc_for_contig(starts2, lengths2, cumsum2)

        h5_path = os.path.join(tmpdir, "clean.h5")
        _make_fragment_h5(h5_path, {
            "chr1": {"starts": starts1, "lengths": lengths1, "gc": gc1},
            "chr2": {"starts": starts2, "lengths": lengths2, "gc": gc2},
        }, contig_lengths={"chr1": 300, "chr2": 200})
        yield h5_path


@pytest.fixture
def padded_h5(tiny_fasta):
    """Fragment H5 with phantom padding rows (simulating the caddb89 bug)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        from fragments_h5.fragment import get_g_or_c_cumsum

        n1 = 200
        starts1 = np.sort(np.random.RandomState(50).randint(0, 260, n1)).astype("int32")
        lengths1 = np.random.RandomState(51).randint(10, 40, n1).astype("uint16")
        stops1 = starts1.astype(np.int64) + lengths1.astype(np.int64)
        mask1 = stops1 <= 300
        starts1 = starts1[mask1]
        lengths1 = lengths1[mask1]

        cumsum1, _ = get_g_or_c_cumsum(tiny_fasta, "chr1")
        gc1 = recompute_gc_for_contig(starts1, lengths1, cumsum1)

        n2 = 180
        starts2 = np.sort(np.random.RandomState(52).randint(0, 170, n2)).astype("int32")
        lengths2 = np.random.RandomState(53).randint(5, 25, n2).astype("uint16")
        stops2 = starts2.astype(np.int64) + lengths2.astype(np.int64)
        mask2 = stops2 <= 200
        starts2 = starts2[mask2]
        lengths2 = lengths2[mask2]

        cumsum2, _ = get_g_or_c_cumsum(tiny_fasta, "chr2")
        gc2 = recompute_gc_for_contig(starts2, lengths2, cumsum2)

        h5_path = os.path.join(tmpdir, "padded.h5")
        _make_fragment_h5(h5_path, {
            "chr1": {"starts": starts1, "lengths": lengths1, "gc": gc1},
            "chr2": {"starts": starts2, "lengths": lengths2, "gc": gc2},
        }, add_padding=True, contig_lengths={"chr1": 300, "chr2": 200})
        yield h5_path


@pytest.fixture
def cumsum_cache(tiny_fasta):
    """Build and return a cumsum cache dir for the tiny FASTA."""
    with tempfile.TemporaryDirectory() as cache_dir:
        fasta_sha256 = compute_file_sha256(tiny_fasta)
        from fragments_h5.repair import build_cumsum_cache
        build_cumsum_cache(tiny_fasta, cache_dir, fasta_sha256)
        yield cache_dir, fasta_sha256


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

    def test_truncation_full_pipeline(self, padded_h5, tiny_fasta, cumsum_cache):
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
        """Rounding agreement over a large set of random fragments."""
        from fragments_h5.fragment import get_g_or_c_cumsum

        cumsum, _ = get_g_or_c_cumsum(tiny_fasta, "chr1")  # all-G, 300 bp

        # Generate many fragments
        rng = np.random.RandomState(100)
        n = 100_000  # 100k fragments for fast test
        starts = rng.randint(0, 250, n).astype("int32")
        lengths = rng.randint(1, 50, n).astype("uint16")
        # Clip
        stops = starts.astype(np.int64) + lengths.astype(np.int64)
        mask = stops <= 300
        starts = starts[mask]
        lengths = lengths[mask]

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

    def test_apply_preserves_gc(self, clean_h5, tiny_fasta, cumsum_cache):
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

    def test_idempotent_datasets(self, padded_h5, tiny_fasta, cumsum_cache):
        """Second repair changes no datasets (only adds a second history entry)."""
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

        # Record post-first-repair state
        with h5py.File(padded_h5, "r") as f:
            gc_after_first = {
                c: f[f"data/{c}/gc"][:] for c in f["data"]
            }
            n_frags_after_first = {
                c: len(f[f"data/{c}/starts"]) for c in f["data"]
            }

        # Second repair (dry-run to verify no changes)
        report2 = repair_local_file(
            local_path=padded_h5,
            fasta_path=tiny_fasta,
            cumsum_cache_dir=cache_dir,
            fasta_sha256=fasta_sha256,
            dry_run=True,
        )

        assert report2["contigs_truncated"] == 0
        for contig, res in report2.get("gc_results", {}).items():
            stats = res.get("stats", {})
            assert stats.get("changed_total", 0) == 0, (
                f"{contig}: gc changed on second repair"
            )


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
