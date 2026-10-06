import logging
import os
import pytest
import subprocess
import sys
import tempfile
import unittest.mock

import numpy
import pysam

from fragments_h5.fragments_h5 import build_fragments_h5, FragmentsH5, bam_to_fragments
import fragments_h5.fragments_h5 as fragments_h5_module

# TODO: Add test coverage for the following missing test cases:
# - MethylCounts and YM tag parsing (fragment.py)
# - Single-end BAM processing (single_end_bam_to_fragments)
# - set_mapq_255_to_none flag
# - allowed_contigs / --contigs filtering
# - Methylation (read_methyl=True)
# - Empty BAM files
# - BAM/FASTA contig length mismatches
# - sequence.pyx Cython module (one_hot_encode_sequences, reverse_complement, etc.)
# - _logging.py utilities

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "./data/")

# importing from datasets was giving me a circular import that i couldn't resolve, so copied this list here
GREATEST_HITS = [
    ("chr6", 99119615, 99119634),  # CTCF Constitutive
]


@pytest.fixture(scope="module")
def bam_path():
    return os.path.join(DATA_DIR, "./small.chr6.bam")


@pytest.fixture(scope="module")
def target_bam_path():
    return os.path.join(DATA_DIR, "./scATAC_breast_v1_chr6_99118615_99121634.hg38.bam")


@pytest.fixture(scope="module")
def fasta_file_path():
    return os.path.join(DATA_DIR, "./GRCh38.p12.genome.chr6_99110000_99130000.fa.gz")


@pytest.fixture(scope="module")
def duplicates_bam_path():
    return os.path.join(DATA_DIR, "./test_duplicates.bam")


@pytest.fixture(scope="module")
def small_h5_path(bam_path, fasta_file_path):
    with tempfile.TemporaryDirectory() as dirname:
        ofname = os.path.join(dirname, os.path.basename(bam_path) + ".frag.h5")
        build_fragments_h5(
            bam_path, ofname, fasta_filename=fasta_file_path
        )
        yield ofname


@pytest.fixture(scope="module")
def target_h5_path(target_bam_path, fasta_file_path):
    with tempfile.TemporaryDirectory() as dirname:
        ofname = os.path.join(
            dirname, "scATAC_breast_v1.chr6_99118615_99121634.fragments.h5"
        )

        cmd = [
            sys.executable, "-m", "fragments_h5.main",
            target_bam_path,
            ofname,
            "--contigs", "chr6",
            "--fasta", fasta_file_path,
            "--verbose",
        ]
        subprocess.run(cmd, check=True)
        yield ofname


def test_build_small_h5(small_h5_path):
    """Test that the fixture builds the fragment h5 successfully."""
    assert os.path.exists(small_h5_path), f"Small H5 file was not created at {small_h5_path}"
    with FragmentsH5(small_h5_path) as fh5:
        assert len(fh5.contig_lengths) > 0, "H5 file should contain at least one contig"


def test_build_target_h5(target_h5_path):
    """Test that the fixture builds the target fragment h5 successfully."""
    assert os.path.exists(target_h5_path), f"Target H5 file was not created at {target_h5_path}"
    with FragmentsH5(target_h5_path) as fh5:
        assert len(fh5.contig_lengths) > 0, "H5 file should contain at least one contig"


def fragments_eq(f1, f2):
    # we can't use the standard fragment __eq__ b/c the gc
    # storage is lossy, and so we can only check that they're close
    return (
        f1.chrom == f2.chrom
        and f1.start == f2.start
        and f1.stop == f2.stop
        and f1.mapq1 == f2.mapq1
        and f1.mapq2 == f2.mapq2
        and (
            (f1.gc is None and f2.gc is None)
            or abs(float(f1.gc) - float(f2.gc)) < 1e-2
        )
    )


def assert_fragments_identical(fs1, fs2):
    assert all(fragments_eq(f1, f2) for f1, f2 in zip(fs1, fs2))


def test_read_all(small_h5_path, bam_path, fasta_file_path):
    # Compare only fragments in the FASTA target region (99110000-99130000) where
    # real sequence data exists. Outside this region the FASTA is N-padded, and
    # float32 cumsum precision loss causes GC divergence between chunk-based
    # builds and full-chromosome computation.
    region_start, region_stop = 99110000, 99130000
    fmh5 = FragmentsH5(small_h5_path, "r")
    h5_fragments = list(fmh5.fetch('chr6', region_start, region_stop, return_gc=True))
    bam_fragments = list(
        bam_to_fragments(
            bam_path, 'chr6', start=region_start, stop=region_stop,
            max_tlen=fmh5.max_fragment_length, fasta_file=fasta_file_path
        )
    )
    assert len(h5_fragments) > 0, "Expected fragments in target region"
    assert_fragments_identical(h5_fragments, bam_fragments)


def test_fetch(small_h5_path, bam_path, fasta_file_path):
    fmh5 = FragmentsH5(small_h5_path, "r")
    for contig, start, stop in GREATEST_HITS:
        h5_fragments = list(fmh5.fetch(contig, start, stop, return_gc=True))
        bam_fragments = list(
            bam_to_fragments(
                bam_path,
                chrom=contig,
                start=start,
                stop=stop,
                max_tlen=fmh5.max_fragment_length,
                fasta_file=fasta_file_path,
            )
        )
        assert_fragments_identical(h5_fragments, bam_fragments)


def test_fetches_over_target_region(target_h5_path, target_bam_path, fasta_file_path):
    fmh5 = FragmentsH5(target_h5_path, "r")

    contig, region_start, region_stop = "chr6", 99119615, 99119634
    for offset in range(0, 501, 100):
        h5_fragments = list(
            fmh5.fetch(contig, region_start - offset, region_stop + offset, return_gc=True)
        )
        bam_fragments = []
        for fragment in bam_to_fragments(
            target_bam_path,
            chrom=contig,
            start=region_start - offset,
            stop=region_stop + offset,
            max_tlen=fmh5.max_fragment_length,
            fasta_file=fasta_file_path,
        ):
            bam_fragments.append(fragment)
        assert_fragments_identical(h5_fragments, bam_fragments)
        print(offset, len(h5_fragments), len(bam_fragments))


def test_fetch_counts(target_h5_path):
    frag_h5 = FragmentsH5(target_h5_path)
    counts = frag_h5.fetch_counts("chr6", 99118615, 99121634)
    assert counts == 184


def test_fetch_cache_v_no_cache(small_h5_path):
    assert list(
        FragmentsH5(small_h5_path, "r", cache_pointers=True).fetch(*GREATEST_HITS[0])
    ) == list(
        FragmentsH5(small_h5_path, "r", cache_pointers=False).fetch(*GREATEST_HITS[0])
    )


def test_context_manager(small_h5_path):
    """Test that FragmentsH5 works as a context manager."""
    with FragmentsH5(small_h5_path) as fh5:
        starts, stops, _ = fh5.fetch_array("chr6")
        assert len(starts) > 0


def test_properties(small_h5_path):
    """Test various properties of FragmentsH5."""
    fh5 = FragmentsH5(small_h5_path)
    
    # Test filename/name properties
    assert fh5.filename == fh5.name
    assert small_h5_path in fh5.filename
    
    # Test n_fragments
    assert fh5.n_fragments > 0
    assert fh5.n_frags == fh5.n_fragments
    
    # Test fragment_length_counts
    assert len(fh5.fragment_length_counts) == fh5.max_fragment_length + 1
    assert fh5.fragment_length_counts.sum() == fh5.n_fragments
    
    fh5.close()


def test_fetch_empty_region(small_h5_path):
    """Test fetching from a region with no fragments."""
    fh5 = FragmentsH5(small_h5_path)
    
    # Query a region far from any data (beginning of chr6)
    starts, stops, supp_data = fh5.fetch_array("chr6", 0, 100)
    assert len(starts) == 0
    assert len(stops) == 0
    
    fh5.close()


def test_fetch_missing_contig(small_h5_path):
    """Test fetching from a contig that doesn't exist raises KeyError."""
    fh5 = FragmentsH5(small_h5_path)
    
    # Query a contig not in the data - should raise KeyError
    with pytest.raises(KeyError):
        fh5.fetch_array("chr_nonexistent", 0, 1000)
    
    fh5.close()


def test_fetch_with_max_frag_len(target_h5_path):
    """Test that max_frag_len filtering works."""
    fh5 = FragmentsH5(target_h5_path)
    
    contig, start, stop = "chr6", 99118615, 99121634
    
    # Fetch all fragments
    starts_all, stops_all, _ = fh5.fetch_array(contig, start, stop)
    
    # Fetch with max_frag_len=200
    starts_filtered, stops_filtered, _ = fh5.fetch_array(
        contig, start, stop, max_frag_len=200
    )
    
    # Filtered should have fewer or equal fragments
    assert len(starts_filtered) <= len(starts_all)
    
    # All filtered fragments should be <= 200 bp
    lengths = stops_filtered - starts_filtered
    assert all(lengths <= 200)
    
    fh5.close()


def test_fetch_with_midpoint_filter(target_h5_path):
    """Test filter_to_midpoint_frags option."""
    fh5 = FragmentsH5(target_h5_path)
    
    contig, start, stop = "chr6", 99119615, 99119634
    
    # Fetch overlapping fragments
    starts_overlap, stops_overlap, _ = fh5.fetch_array(
        contig, start, stop, filter_to_midpoint_frags=False
    )
    
    # Fetch midpoint-filtered fragments
    starts_midpoint, stops_midpoint, _ = fh5.fetch_array(
        contig, start, stop, filter_to_midpoint_frags=True
    )
    
    # Midpoint-filtered should have fewer or equal fragments
    assert len(starts_midpoint) <= len(starts_overlap)
    
    # All midpoint-filtered fragments should have midpoints in [start, stop)
    midpoints = starts_midpoint + (stops_midpoint - starts_midpoint) // 2
    assert all((midpoints >= start) & (midpoints < stop))
    
    fh5.close()


def test_fetch_supplementary_data(target_h5_path):
    """Test fetching with supplementary data options."""
    fh5 = FragmentsH5(target_h5_path)
    
    contig, start, stop = "chr6", 99119615, 99119634
    
    # Test return_mapqs
    starts, stops, supp = fh5.fetch_array(
        contig, start, stop, return_mapqs=True
    )
    assert "mapq" in supp
    assert supp["mapq"].shape == (len(starts), 2)
    
    # Test return_gc
    starts, stops, supp = fh5.fetch_array(
        contig, start, stop, return_gc=True
    )
    assert "gc" in supp
    assert len(supp["gc"]) == len(starts)
    
    # Test return_strand
    starts, stops, supp = fh5.fetch_array(
        contig, start, stop, return_strand=True
    )
    assert "strand" in supp
    assert len(supp["strand"]) == len(starts)
    
    fh5.close()


def test_region_beyond_contig_raises(small_h5_path):
    """Test that querying beyond contig end raises ValueError."""
    fh5 = FragmentsH5(small_h5_path)
    
    # Get the contig length
    contig_len = fh5.contig_lengths["chr6"]
    
    # Query starting beyond the contig should raise
    with pytest.raises(ValueError, match="beyond the contig end"):
        fh5.fetch_array("chr6", contig_len + 1000, contig_len + 2000)
    
    fh5.close()


def test_pickle_support(small_h5_path):
    """Test that FragmentsH5 can be pickled (for multiprocessing)."""
    import pickle
    
    fh5 = FragmentsH5(small_h5_path)
    
    # Get some data before pickling
    starts_before, stops_before, _ = fh5.fetch_array("chr6")
    
    # Pickle and unpickle
    pickled = pickle.dumps(fh5)
    fh5_restored = pickle.loads(pickled)
    
    # Verify data is the same after unpickling
    starts_after, stops_after, _ = fh5_restored.fetch_array("chr6")
    
    assert list(starts_before) == list(starts_after)
    assert list(stops_before) == list(stops_after)
    
    fh5.close()
    fh5_restored.close()


def _reference_has_flags(h5_path):
    """Recompute the four structural flags straight from raw h5py.

    Deliberately a transcription of the pre-cache `any(...)` expressions rather than a
    call into FragmentsH5, so that these tests cannot pass by agreeing with the code
    they are meant to check.
    """
    import h5py

    with h5py.File(h5_path, "r") as f:
        data = f["data"]
        return {
            "has_methyl": any("num_cpgs" in data[c] for c in data.keys()),
            "has_strand": any(
                ("strand" in data[c]) and (len(data[c]["strand"].shape) == 1)
                for c in data.keys()
            ),
            "has_gc": any("gc" in data[c] for c in data.keys()),
            "has_fragment_end_clipped": any(
                "fragment_end_clipped" in data[c] for c in data.keys()
            ),
        }


def _actual_has_flags(fh5):
    return {
        "has_methyl": fh5.has_methyl,
        "has_strand": fh5.has_strand,
        "has_gc": fh5.has_gc,
        "has_fragment_end_clipped": fh5.has_fragment_end_clipped,
    }


@pytest.fixture(scope="module")
def no_gc_h5_path(bam_path):
    """A build with no --fasta, so `gc` is absent -> exercises the False branch."""
    with tempfile.TemporaryDirectory() as dirname:
        ofname = os.path.join(dirname, "no_gc.frag.h5")
        build_fragments_h5(bam_path, ofname)
        yield ofname


@pytest.fixture(scope="module")
def non_uniform_h5_path(many_contig_h5_path):
    """A 30-contig h5 where each optional dataset exists on only SOME contigs.

    The builder cannot produce this: it applies read_gc / read_strand / read_methyl /
    store_fragment_end_clipped globally to every contig via SubBuildArgs, so every
    builder-made file is uniform and `any()` and `all()` agree on it. That uniformity let
    an `any()` -> `all()` mutant survive the whole suite, so the non-uniform file has to be
    forged with raw h5py (the same tactic as two_bit_strand_h5_path).

    This is the exact shape AGENT_CONTEXT.md 7.3 flags as a live, unfixed hazard: a file
    with `gc` on some contigs but not others reports has_gc == True, and a consumer that
    then iterates all contigs gets a KeyError. The tests below pin the CURRENT `any()`
    semantics because this change is behaviour-preserving -- they do not assert that
    `any()` is the right answer.
    """
    import h5py
    import shutil

    with tempfile.TemporaryDirectory() as dirname:
        ofname = os.path.join(dirname, "non_uniform.frag.h5")
        shutil.copy(many_contig_h5_path, ofname)
        with h5py.File(ofname, "r+") as f:
            contigs = sorted(f["data"].keys())
            assert len(contigs) >= 3, "need several contigs for non-uniformity to exist"
            first, second = contigs[0], contigs[1]
            n_first = f["data"][first]["starts"].shape[0]

            # present on exactly ONE contig -> any() True, all() False
            f["data"][first].create_dataset(
                "gc", data=numpy.zeros(n_first, dtype="uint8"), dtype="uint8"
            )
            f["data"][first].create_dataset(
                "num_cpgs", data=numpy.zeros(n_first, dtype="uint8"), dtype="uint8"
            )
            f["data"][first].create_dataset(
                "fragment_end_clipped",
                data=numpy.zeros(n_first, dtype="uint8"),
                dtype="uint8",
            )
            # strand exists on every contig from the build; remove it from one
            # -> still any() True, but all() False
            if "strand" in f["data"][second]:
                del f["data"][second]["strand"]
        yield ofname


def test_has_properties_use_any_not_all_semantics(non_uniform_h5_path):
    """On a non-uniform file all four must report True, i.e. any() not all().

    This is the assertion that kills an `any()` -> `all()` mutant, which otherwise survives
    the entire suite because every builder-produced fixture is uniform.

    Pins current behaviour, which AGENT_CONTEXT.md 7.3 records as a known hazard rather
    than a desirable design: a consumer trusting has_gc == True and then iterating every
    contig will KeyError on the contigs that lack it. Switching to all() would be a real
    behaviour change and needs its own decision -- it is not in scope for a caching change.
    """
    import h5py

    with h5py.File(non_uniform_h5_path, "r") as f:
        contigs = sorted(f["data"].keys())
        n = len(contigs)
        with_gc = sum("gc" in f["data"][c] for c in contigs)
        with_strand = sum("strand" in f["data"][c] for c in contigs)

    # the fixture is only meaningful if it is genuinely non-uniform
    assert 0 < with_gc < n, f"gc on {with_gc}/{n} contigs; fixture is not non-uniform"
    assert 0 < with_strand < n, f"strand on {with_strand}/{n} contigs; not non-uniform"

    fh5 = FragmentsH5(non_uniform_h5_path)
    try:
        assert _actual_has_flags(fh5) == {
            "has_methyl": True,
            "has_strand": True,
            "has_gc": True,
            "has_fragment_end_clipped": True,
        }
    finally:
        fh5.close()


@pytest.fixture(scope="module")
def methyl_h5_path():
    """A build with read_methyl=True, so `num_cpgs` exists -> has_methyl is True.

    Without this, every fixture has has_methyl == False, and a mutant that walks the
    contigs but ignores `num_cpgs` (`any(False for contig in ...)`) passes the whole
    suite. Methylation comes from the YM tag, whose format is fixed by
    MethylCounts.init_from_ym_tag in fragment.py.
    """
    ym = (
        "unconverted_cytosines:3;converted_cytosines:7;"
        "unconverted_cpgs:2;converted_cpgs:5"
    )
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 200_000}],
    }
    with tempfile.TemporaryDirectory() as dirname:
        bam = os.path.join(dirname, "methyl.bam")
        seq_len, tlen = 10, 120
        with pysam.AlignmentFile(bam, "wb", header=header) as outf:
            for i in range(5):
                pos1 = 1_000 + i * 1_000
                pos2 = pos1 + tlen - seq_len
                for flag, start, mate, tl in (
                    (0x1 | 0x2 | 0x20 | 0x40, pos1, pos2, tlen),
                    (0x1 | 0x2 | 0x10 | 0x80, pos2, pos1, -tlen),
                ):
                    a = pysam.AlignedSegment()
                    a.query_name = f"methyl_read_{i}"
                    a.reference_id = 0
                    a.reference_start = start
                    a.cigarstring = f"{seq_len}M"
                    a.mapping_quality = 60
                    a.query_sequence = "A" * seq_len
                    a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
                    a.flag = flag
                    a.next_reference_id = 0
                    a.next_reference_start = mate
                    a.template_length = tl
                    a.set_tag("YM", ym)
                    outf.write(a)
        pysam.sort("-o", bam, bam)
        pysam.index(bam)

        ofname = os.path.join(dirname, "methyl.frag.h5")
        build_fragments_h5(
            bam, ofname, num_processes=1, read_strand=True, read_methyl=True,
            store_fragment_end_clipped=False,
        )
        yield ofname


@pytest.fixture(scope="module")
def no_strand_h5_path(bam_path):
    """A build with read_strand=False, so `strand` is absent -> has_strand is False.

    Without this, every fixture has has_strand == True and an `any(True for ...)` mutant
    survives the suite.
    """
    with tempfile.TemporaryDirectory() as dirname:
        ofname = os.path.join(dirname, "no_strand.frag.h5")
        build_fragments_h5(bam_path, ofname, read_strand=False)
        yield ofname


@pytest.fixture(scope="module")
def two_bit_strand_h5_path(bam_path):
    """A build whose `strand` dataset is 2-D, which must read as has_strand == False.

    This is the only non-trivial logic in has_strand: the `len(shape) == 1` guard exists
    because some very old "small frag" h5s stored strand in two bits. No fixture
    exercised that guard, so a mutant ignoring it survived. Built normally, then the
    strand dataset is replaced with a 2-column one via raw h5py -- the builder cannot
    produce this shape any more, which is exactly why it has to be forged here.
    """
    import h5py

    with tempfile.TemporaryDirectory() as dirname:
        ofname = os.path.join(dirname, "two_bit_strand.frag.h5")
        build_fragments_h5(bam_path, ofname, read_strand=True)
        with h5py.File(ofname, "r+") as f:
            for contig in f["data"]:
                grp = f["data"][contig]
                if "strand" not in grp:
                    continue
                n = grp["strand"].shape[0]
                del grp["strand"]
                grp.create_dataset(
                    "strand", data=numpy.zeros((n, 2), dtype="uint8"), dtype="uint8"
                )
        yield ofname


@pytest.fixture(scope="module")
def many_contig_h5_path():
    """A fragment h5 spanning many contigs, so a per-contig scan is clearly visible.

    Real production files have ~195 contigs (hg38 with alts/randoms); 30 is enough to
    make a per-contig-per-access scan unmistakable while staying fast to build.
    """
    n_contigs = 30
    contigs = [(f"chr{i}", 200_000) for i in range(1, n_contigs + 1)]
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    }
    with tempfile.TemporaryDirectory() as dirname:
        bam = os.path.join(dirname, "many_contig.bam")
        seq_len, tlen = 10, 120
        with pysam.AlignmentFile(bam, "wb", header=header) as outf:
            for contig, _ in contigs:
                tid = outf.get_tid(contig)
                for i in range(5):
                    pos1 = 1_000 + i * 1_000
                    pos2 = pos1 + tlen - seq_len
                    for flag, start, mate in (
                        (0x1 | 0x2 | 0x20 | 0x40, pos1, pos2),
                        (0x1 | 0x2 | 0x10 | 0x80, pos2, pos1),
                    ):
                        a = pysam.AlignedSegment()
                        a.query_name = f"{contig}_read_{i}"
                        a.reference_id = tid
                        a.reference_start = start
                        a.cigarstring = f"{seq_len}M"
                        a.mapping_quality = 60
                        a.query_sequence = "A" * seq_len
                        a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
                        a.flag = flag
                        a.next_reference_id = tid
                        a.next_reference_start = mate
                        a.template_length = tlen if start == pos1 else -tlen
                        outf.write(a)
        pysam.sort("-o", bam, bam)
        pysam.index(bam)

        ofname = os.path.join(dirname, "many_contig.frag.h5")
        build_fragments_h5(
            bam, ofname, num_processes=1, read_strand=True,
            store_fragment_end_clipped=False,
        )
        yield ofname


@pytest.mark.parametrize(
    "fixture_name",
    [
        "small_h5_path",
        "target_h5_path",
        "no_gc_h5_path",
        "many_contig_h5_path",
        "methyl_h5_path",
        "no_strand_h5_path",
        "two_bit_strand_h5_path",
    ],
)
def test_has_properties_match_independent_reference(request, fixture_name):
    """Caching must not change what any has_* property returns."""
    h5_path = request.getfixturevalue(fixture_name)
    expected = _reference_has_flags(h5_path)

    fh5 = FragmentsH5(h5_path)
    try:
        assert _actual_has_flags(fh5) == expected
    finally:
        fh5.close()


# The three guards below exist because a value test is only as strong as its fixtures'
# disagreement. Each flag needs at least one True and one False fixture, or a mutant that
# hardcodes the majority answer passes every assertion. Two such mutants were found by
# executing them: `has_methyl -> any(False for ...)` and `has_strand -> any(True for ...)`.


def test_has_gc_differs_between_fixtures(small_h5_path, no_gc_h5_path):
    """Guard that the parametrized test above spans both a True and a False gc case."""
    assert _reference_has_flags(small_h5_path)["has_gc"] is True
    assert _reference_has_flags(no_gc_h5_path)["has_gc"] is False


def test_has_methyl_differs_between_fixtures(methyl_h5_path, small_h5_path):
    """Pin that some fixture has has_methyl True, else an always-False mutant survives."""
    assert _reference_has_flags(methyl_h5_path)["has_methyl"] is True
    assert _reference_has_flags(small_h5_path)["has_methyl"] is False


def test_has_fragment_end_clipped_differs_between_fixtures(
    small_h5_path, many_contig_h5_path
):
    """Pin that some fixture has has_fragment_end_clipped False.

    The spread was implicit (many_contig_h5_path builds with
    store_fragment_end_clipped=False) but nothing asserted it, so a fixture change could
    have silently removed the only False case and let a hardcoded mutant through.
    """
    assert _reference_has_flags(small_h5_path)["has_fragment_end_clipped"] is True
    assert _reference_has_flags(many_contig_h5_path)["has_fragment_end_clipped"] is False


def test_has_strand_differs_between_fixtures(
    small_h5_path, no_strand_h5_path, two_bit_strand_h5_path
):
    """Pin that some fixture has has_strand False, else an always-True mutant survives.

    Covers both False routes: strand absent entirely, and strand present but 2-D (the
    legacy two-bit layout the `len(shape) == 1` guard exists for).
    """
    assert _reference_has_flags(small_h5_path)["has_strand"] is True
    assert _reference_has_flags(no_strand_h5_path)["has_strand"] is False
    assert _reference_has_flags(two_bit_strand_h5_path)["has_strand"] is False


def _count_h5py_group_access(monkeypatch):
    """Count h5py Group lookups. Returns a mutable counter.

    `get` and `keys` are counted as well as `__contains__`/`__getitem__`: patching only
    the dunders leaves `Group.get()` as a silent route to the file, so an implementation
    using `.get()` would read as zero accesses. The counter should be hard to evade, not
    merely sufficient for the current implementation.
    """
    import h5py

    counts = {"contains": 0, "getitem": 0, "get": 0, "keys": 0}
    orig_contains = h5py.Group.__contains__
    orig_getitem = h5py.Group.__getitem__
    orig_get = h5py.Group.get
    orig_keys = h5py.Group.keys

    def counting_contains(self, key):
        counts["contains"] += 1
        return orig_contains(self, key)

    def counting_getitem(self, key):
        counts["getitem"] += 1
        return orig_getitem(self, key)

    def counting_get(self, key, *args, **kwargs):
        counts["get"] += 1
        return orig_get(self, key, *args, **kwargs)

    def counting_keys(self):
        counts["keys"] += 1
        return orig_keys(self)

    monkeypatch.setattr(h5py.Group, "__contains__", counting_contains)
    monkeypatch.setattr(h5py.Group, "__getitem__", counting_getitem)
    monkeypatch.setattr(h5py.Group, "get", counting_get)
    monkeypatch.setattr(h5py.Group, "keys", counting_keys)
    return counts


def test_has_properties_scan_once(small_h5_path, monkeypatch):
    """The structural scan happens once per open handle, not once per access.

    Before this was cached, each has_* access walked every contig, issuing one
    Group.__contains__ and one Group.__getitem__ per contig per access. Now repeated
    access must touch the HDF5 file zero additional times.
    """
    counts = _count_h5py_group_access(monkeypatch)

    fh5 = FragmentsH5(small_h5_path)
    try:
        # first access resolves each answer; that one scan is the cost we are keeping
        _actual_has_flags(fh5)
        counts.update(dict.fromkeys(counts, 0))

        for _ in range(25):
            fh5.has_methyl
            fh5.has_strand
            fh5.has_gc
            fh5.has_fragment_end_clipped

        assert all(v == 0 for v in counts.values()), (
            f"has_* properties re-scanned the h5 on access: {counts}"
        )
    finally:
        fh5.close()


def test_scan_cost_is_independent_of_access_count(many_contig_h5_path, monkeypatch):
    """The scan cost must not scale with the number of accesses.

    Uses a many-contig file, where the pre-cache cost was one Group.__contains__ plus
    one Group.__getitem__ *per contig per access*. has_methyl is the worst case: it is
    False here, so `any()` could not short-circuit and walked every contig every time.
    """
    import h5py

    with h5py.File(many_contig_h5_path, "r") as f:
        n_contigs = len(f["data"].keys())
    assert n_contigs >= 20, f"fixture only has {n_contigs} contigs; test would be weak"

    fh5 = FragmentsH5(many_contig_h5_path)
    try:
        assert fh5.has_methyl is False, "has_methyl must be the non-short-circuiting case"

        counts = _count_h5py_group_access(monkeypatch)
        for _ in range(40):
            fh5.has_methyl
        assert all(v == 0 for v in counts.values()), (
            f"40 has_methyl accesses on a {n_contigs}-contig file cost {counts}; "
            f"the pre-cache code would have cost ~{40 * n_contigs} of each"
        )
    finally:
        fh5.close()


def test_first_scan_short_circuits_per_flag(many_contig_h5_path, monkeypatch):
    """Each flag must stop at the first contig that settles it.

    This bounds the cost of the FIRST scan. Every other caching test installs its
    counters *after* the answer is already cached, so none of them can tell a
    short-circuiting scan from a full walk. That gap was real: wrapping the generator
    in a list -- `any([...])` instead of `any(...)` -- returns the identical value for
    every file while visiting every contig, and passed all 20 other tests.

    The per-flag short-circuit is the whole reason this is `cached_property` per flag
    rather than one eager scan in `__init__`. An eager combined pass has to resolve the
    worst-case flag, which measured 20 -> 256 ms at open on a 195-contig file. A mutant
    that silently removes the short-circuit re-creates that regression with no failing
    test, so the cost asymmetry needs pinning directly.

    On this fixture `strand` is present on every contig, so `has_strand` settles True at
    contig #1. `num_cpgs` is absent everywhere, so `has_methyl` can only settle False by
    visiting all of them. The gap between those two costs IS the short-circuit.
    """
    import h5py

    with h5py.File(many_contig_h5_path, "r") as f:
        n_contigs = len(f["data"].keys())
    assert n_contigs >= 20, f"fixture has only {n_contigs} contigs; test would be weak"

    counts = _count_h5py_group_access(monkeypatch)

    # has_strand is True at the first contig, so it must not walk the file.
    fh5 = FragmentsH5(many_contig_h5_path)
    try:
        counts["contains"] = 0
        assert fh5.has_strand is True
        short_circuit_cost = counts["contains"]
    finally:
        fh5.close()

    # has_methyl is False everywhere, so it has no choice but to visit every contig.
    fh5 = FragmentsH5(many_contig_h5_path)
    try:
        counts["contains"] = 0
        assert fh5.has_methyl is False
        full_walk_cost = counts["contains"]
    finally:
        fh5.close()

    assert short_circuit_cost <= 2, (
        f"has_strand cost {short_circuit_cost} Group.__contains__ calls on a "
        f"{n_contigs}-contig file, but it is True at contig #1 and must stop there. "
        f"`any([...])` in place of `any(...)` produces exactly this symptom: same "
        f"return value, every contig visited."
    )
    assert full_walk_cost >= n_contigs, (
        f"has_methyl cost only {full_walk_cost} calls; it is False on every contig "
        f"and must therefore visit all {n_contigs} of them. A cost below that means "
        f"it is not really scanning."
    )


def test_has_properties_survive_pickle(many_contig_h5_path, monkeypatch):
    """Unpickling reopens the same filename, so the cached answers stay valid.

    They are carried in __dict__ by __getstate__ on purpose: a forked worker should
    inherit them rather than redo the scan. Uses the many-contig fixture so that a
    rescan (~one Group.__contains__ per contig) is distinguishable from __setstate__'s
    two unavoidable reopen lookups ( _f["index"] and _f["data"] ).
    """
    import pickle
    import h5py

    with h5py.File(many_contig_h5_path, "r") as f:
        n_contigs = len(f["data"].keys())

    expected = _reference_has_flags(many_contig_h5_path)

    fh5 = FragmentsH5(many_contig_h5_path)
    _actual_has_flags(fh5)  # resolve before pickling, so the answers are in __dict__
    pickled = pickle.dumps(fh5)
    fh5.close()

    counts = _count_h5py_group_access(monkeypatch)
    restored = pickle.loads(pickled)
    try:
        assert _actual_has_flags(restored) == expected
        # a structural rescan would cost one __contains__ per contig; reopening costs none
        assert counts["contains"] == 0, (
            f"unpickling re-scanned the h5 structure: {counts} "
            f"(a rescan of {n_contigs} contigs is what this rules out)"
        )
        # __setstate__ costs exactly two lookups: _f["index"] and _f["data"]. Bounding at
        # n_contigs instead would still admit a partial rescan of a short-circuiting flag.
        assert counts["getitem"] <= 2, (
            f"unpickling issued {counts['getitem']} group lookups; __setstate__ needs "
            f'exactly 2 (_f["index"], _f["data"]). Full counts: {counts}'
        )
        assert counts["get"] == 0 and counts["keys"] == 0, (
            f"unpickling reached the file by get()/keys(): {counts}"
        )
    finally:
        restored.close()


def test_resolved_has_properties_outlive_close(small_h5_path):
    """An answer computed while open stays readable after close().

    It describes the file's structure, not the liveness of the handle -- the same
    reason contig_lengths / max_fragment_length / n_fragments are readable after
    close(). This is the one behavioural change: pre-cache, every post-close access
    raised ValueError("Invalid group (or file) id").
    """
    expected = _reference_has_flags(small_h5_path)

    fh5 = FragmentsH5(small_h5_path)
    resolved = _actual_has_flags(fh5)
    fh5.close()

    assert _actual_has_flags(fh5) == resolved == expected
    # the metadata whose post-close behaviour this now matches
    assert fh5.contig_lengths["chr6"] > 0
    assert fh5.n_fragments > 0


@pytest.mark.parametrize(
    "prop_name", ["has_methyl", "has_strand", "has_gc", "has_fragment_end_clipped"]
)
def test_unresolved_has_property_after_close_still_raises(small_h5_path, prop_name):
    """A never-computed answer still raises on a dead handle, as it always did.

    Pins the limit of the change: caching only relaxes post-close behaviour for
    answers already known. It never makes a previously-working access fail, and it
    does not resurrect a dead handle.

    This also happens to be the assertion that kills a naive `return False` mutant --
    such a mutant never touches the file, so it answers on a dead handle instead of
    raising. Parametrized over all four so that holds for each of them.
    """
    fh5 = FragmentsH5(small_h5_path)
    fh5.close()

    with pytest.raises(ValueError, match="Invalid group"):
        getattr(fh5, prop_name)


def test_include_duplicates(duplicates_bam_path, fasta_file_path):
    """Test that include_duplicates parameter correctly includes/excludes duplicate-marked fragments."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Build h5 excluding duplicates (default)
        h5_no_dups = os.path.join(tmpdir, "no_dups.h5")
        build_fragments_h5(
            duplicates_bam_path, h5_no_dups, fasta_filename=fasta_file_path,
            include_duplicates=False, num_processes=1
        )
        
        # Build h5 including duplicates
        h5_with_dups = os.path.join(tmpdir, "with_dups.h5")
        build_fragments_h5(
            duplicates_bam_path, h5_with_dups, fasta_filename=fasta_file_path,
            include_duplicates=True, num_processes=1
        )
        
        # Compare fragment counts
        fh5_no_dups = FragmentsH5(h5_no_dups)
        fh5_with_dups = FragmentsH5(h5_with_dups)
        
        # Without duplicates: should have 1 fragment (the non-duplicate)
        assert fh5_no_dups.n_fragments == 1
        
        # With duplicates: should have 2 fragments (both the original and duplicate)
        assert fh5_with_dups.n_fragments == 2
        
        fh5_no_dups.close()
        fh5_with_dups.close()


def test_fragment_end_clipped_storage_and_read(
    bam_path, fasta_file_path
):
    """Test that fragment_end_clipped is stored/omitted based on args and read/counted correctly.

    - Build with store_fragment_end_clipped=True (default): file has the dataset; counts from
      fetch_array(return_fragment_end_clipped=True) match counts from bam_to_fragments.
    - Build with store_fragment_end_clipped=False: file does not have the dataset; requesting
      return_fragment_end_clipped=True raises.
    """
    import numpy as np

    with tempfile.TemporaryDirectory() as tmpdir:
        # --- Ground truth: count fragment_end_clipped from BAM ---
        bam_clipped_true = 0
        bam_clipped_false = 0
        bam_clipped_unknown = 0
        bam_total = 0
        for frag in bam_to_fragments(
            bam_path,
            "chr6",
            max_tlen=65535,
            fasta_file=fasta_file_path,
        ):
            bam_total += 1
            if frag.fragment_end_clipped is True:
                bam_clipped_true += 1
            elif frag.fragment_end_clipped is False:
                bam_clipped_false += 1
            else:
                bam_clipped_unknown += 1

        assert bam_total > 0, "test BAM should have at least one fragment on chr6"

        # --- Build H5 with store_fragment_end_clipped=True (default) ---
        h5_with_clipped = os.path.join(tmpdir, "with_clipped.h5")
        build_fragments_h5(
            bam_path,
            h5_with_clipped,
            fasta_filename=fasta_file_path,
            store_fragment_end_clipped=True,
            num_processes=1,
        )

        fh5_with = FragmentsH5(h5_with_clipped)
        assert fh5_with.has_fragment_end_clipped is True

        starts, stops, supp = fh5_with.fetch_array(
            "chr6", return_fragment_end_clipped=True
        )
        assert "fragment_end_clipped" in supp
        fec = supp["fragment_end_clipped"]
        assert fec.dtype in (np.uint8, np.dtype("uint8"))
        assert len(fec) == len(starts) == fh5_with.n_fragments

        read_clipped_true = int((fec == 1).sum())
        read_clipped_false = int((fec == 0).sum())
        read_clipped_unknown = int((fec == 255).sum())

        assert read_clipped_true == bam_clipped_true
        assert read_clipped_false == bam_clipped_false
        assert read_clipped_unknown == bam_clipped_unknown
        assert read_clipped_true + read_clipped_false + read_clipped_unknown == bam_total

        # --- Build H5 with store_fragment_end_clipped=False ---
        h5_without_clipped = os.path.join(tmpdir, "without_clipped.h5")
        build_fragments_h5(
            bam_path,
            h5_without_clipped,
            fasta_filename=fasta_file_path,
            store_fragment_end_clipped=False,
            num_processes=1,
        )

        fh5_without = FragmentsH5(h5_without_clipped)
        assert fh5_without.has_fragment_end_clipped is False

        with pytest.raises(ValueError, match="does not contain fragment_end_clipped"):
            fh5_without.fetch_array("chr6", return_fragment_end_clipped=True)

        # --- fetch() returns correct fragment_end_clipped when present ---
        with_clipped_frags = list(
            fh5_with.fetch(return_fragment_end_clipped=True)
        )
        assert len(with_clipped_frags) == bam_total
        assert sum(1 for f in with_clipped_frags if f.fragment_end_clipped is True) == bam_clipped_true
        assert sum(1 for f in with_clipped_frags if f.fragment_end_clipped is False) == bam_clipped_false
        assert sum(1 for f in with_clipped_frags if f.fragment_end_clipped is None) == bam_clipped_unknown

        fh5_with.close()
        fh5_without.close()


def test_build_fails_if_output_exists(bam_path, fasta_file_path):
    """Test that build_fragments_h5 raises when the output file already exists."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ofname = os.path.join(tmpdir, "existing.h5")
        build_fragments_h5(
            bam_path, ofname, fasta_filename=fasta_file_path, num_processes=1
        )
        assert os.path.exists(ofname)

        # Building again to the same path should fail (h5py opens with mode 'x')
        with pytest.raises(FileExistsError):
            build_fragments_h5(
                bam_path, ofname, fasta_filename=fasta_file_path, num_processes=1
            )


@pytest.mark.timeout(30)
def test_multiprocessing_with_small_bam(duplicates_bam_path, fasta_file_path):
    """Test that multiprocessing works correctly with small BAMs.

    This test specifically validates that using multiple workers with a tiny BAM
    (only 4 reads, 1 contig) doesn't cause hangs or deadlocks. This scenario
    previously caused issues with forkserver due to race conditions between
    server initialization and quick job completion.

    We use more workers than contigs (8 workers, 1 contig) to stress-test
    the edge case where some workers never get jobs.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        output_h5 = os.path.join(tmpdir, "multiproc_test.h5")

        # Use 8 workers with a BAM that has only 1 contig
        # This is the worst-case scenario: many idle workers
        build_fragments_h5(
            duplicates_bam_path,
            output_h5,
            fasta_filename=fasta_file_path,
            include_duplicates=True,
            num_processes=8  # Much more than needed
        )

        # Verify output is correct
        fh5 = FragmentsH5(output_h5)
        assert fh5.n_fragments == 2
        fh5.close()


@pytest.mark.timeout(60)
def test_multiprocessing_stress_test(bam_path, fasta_file_path):
    """Stress test multiprocessing by running multiple builds in sequence.

    This simulates running multiple tests back-to-back (as pytest does)
    to ensure there are no state/cleanup issues between runs.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Run 5 builds in a row with different worker counts
        for i, num_procs in enumerate([2, 4, 8, 1, 4]):
            output_h5 = os.path.join(tmpdir, f"stress_test_{i}.h5")

            build_fragments_h5(
                bam_path,
                output_h5,
                fasta_filename=fasta_file_path,
                num_processes=num_procs
            )

            # Quick validation
            fh5 = FragmentsH5(output_h5)
            assert fh5.n_fragments > 0
            fh5.close()


@pytest.mark.timeout(30)
def test_multiprocessing_pool_is_actually_used(bam_path, fasta_file_path, caplog):
    """Verify that num_processes > 1 actually takes the multiprocessing code path."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_h5 = os.path.join(tmpdir, "pool_check.h5")

        with caplog.at_level(logging.INFO, logger="fragments_h5.fragments_h5"):
            build_fragments_h5(
                bam_path,
                output_h5,
                fasta_filename=fasta_file_path,
                num_processes=4,
            )

        assert any(
            "Using multiprocessing with 4 workers" in msg for msg in caplog.messages
        ), (
            f"Expected multiprocessing log message but got: {caplog.messages}"
        )

        fh5 = FragmentsH5(output_h5)
        assert fh5.n_fragments > 0
        fh5.close()


@pytest.mark.timeout(30)
def test_single_process_path_is_used(bam_path, fasta_file_path, caplog):
    """Verify that num_processes=1 takes the single-process code path."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_h5 = os.path.join(tmpdir, "single_check.h5")

        with caplog.at_level(logging.INFO, logger="fragments_h5.fragments_h5"):
            build_fragments_h5(
                bam_path,
                output_h5,
                fasta_filename=fasta_file_path,
                num_processes=1,
            )

        assert any(
            "Using single-process path" in msg for msg in caplog.messages
        ), (
            f"Expected single-process log message but got: {caplog.messages}"
        )

        fh5 = FragmentsH5(output_h5)
        assert fh5.n_fragments > 0
        fh5.close()


def _build_and_compare_chunk_sizes(bam_path, fasta_file_path, chunk_size, num_processes=1):
    """Build two H5 files — one with default chunk size, one with a small chunk size —
    and verify all datasets are identical."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Build reference with large chunk size (entire contig = single chunk)
        ref_h5_path = os.path.join(tmpdir, "reference.h5")
        build_fragments_h5(
            bam_path, ref_h5_path, fasta_filename=fasta_file_path,
            num_processes=1,
        )

        # Build with small chunk size to force multiple chunks per contig
        chunked_h5_path = os.path.join(tmpdir, "chunked.h5")
        with unittest.mock.patch.object(fragments_h5_module, 'GENOMIC_CHUNK_SIZE', chunk_size):
            build_fragments_h5(
                bam_path, chunked_h5_path, fasta_filename=fasta_file_path,
                num_processes=num_processes,
            )

        # Compare all datasets
        with FragmentsH5(ref_h5_path) as ref, FragmentsH5(chunked_h5_path) as chunked:
            assert ref.n_fragments == chunked.n_fragments, (
                f"Fragment count mismatch: {ref.n_fragments} vs {chunked.n_fragments}"
            )
            # Compare contigs that have data (not all contigs in contig_lengths
            # will have data groups — only those with mapped reads)
            ref_data_contigs = set(ref._f['data'].keys())
            chunked_data_contigs = set(chunked._f['data'].keys())
            assert ref_data_contigs == chunked_data_contigs, (
                f"Data contig mismatch: {ref_data_contigs} vs {chunked_data_contigs}"
            )

            for contig in ref_data_contigs:
                ref_data = ref._f[f'data/{contig}']
                chunked_data = chunked._f[f'data/{contig}']

                assert set(ref_data.keys()) == set(chunked_data.keys()), (
                    f"Dataset keys mismatch for {contig}"
                )

                for ds_name in ref_data.keys():
                    ref_arr = ref_data[ds_name][:]
                    chunked_arr = chunked_data[ds_name][:]
                    numpy.testing.assert_array_equal(
                        ref_arr, chunked_arr,
                        err_msg=f"Data mismatch for {contig}/{ds_name}"
                    )

                # Verify index if it exists
                if contig in ref._f.get('index', {}):
                    assert contig in chunked._f['index'], (
                        f"Index missing for {contig} in chunked build"
                    )
                    numpy.testing.assert_array_equal(
                        ref._f[f'index/{contig}'][:],
                        chunked._f[f'index/{contig}'][:],
                        err_msg=f"Index mismatch for {contig}"
                    )


@pytest.mark.timeout(120)
def test_chunk_merge_correctness_single_process(bam_path, fasta_file_path):
    """Verify that splitting a contig into small chunks and merging produces
    identical results to processing the whole contig at once.

    Uses 5M chunk size: chr6 (170M bases) splits into ~34 chunks.
    Reads are in a ~20kb window around position 99.1M, so 1-2 chunks have data.
    Tests the merge path including concatenation of chunk arrays."""
    _build_and_compare_chunk_sizes(bam_path, fasta_file_path, chunk_size=10_000_000)


@pytest.mark.timeout(120)
def test_chunk_merge_correctness_multiprocess(bam_path, fasta_file_path):
    """Same as above but with multiprocessing to test concurrent chunk processing."""
    _build_and_compare_chunk_sizes(bam_path, fasta_file_path, chunk_size=5_000_000, num_processes=4)


@pytest.mark.timeout(120)
def test_chunk_merge_small_chunks(bam_path, fasta_file_path):
    """Test with 1M chunks to create more chunks and increase boundary crossings.
    The ~20kb data window at position ~99.1M will be split across 1-2 chunks."""
    _build_and_compare_chunk_sizes(bam_path, fasta_file_path, chunk_size=10_000_000)


@pytest.mark.timeout(120)
def test_chunk_merge_with_target_bam(target_bam_path, fasta_file_path):
    """Test chunk merge with the larger scATAC BAM file."""
    _build_and_compare_chunk_sizes(target_bam_path, fasta_file_path, chunk_size=10_000_000)
