import os
import pytest
import subprocess
import tempfile
import pysam

from fragments_h5.fragments_h5 import build_fragments_h5, FragmentsH5, bam_to_fragments

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

        cmd = f"""
        build-fragments-h5 \
            {target_bam_path} \
            {ofname} \
            --contigs chr6 \
            --fasta {fasta_file_path} \
            --verbose \
        """.strip()
        subprocess.run(cmd, shell=True, check=True)
        yield ofname


def test_build_small_h5(small_h5_path):
    # just test that the fixture builds the fragment h5
    # this also helps to make the timings more interpretable
    pass


def test_build_target_h5(target_h5_path):
    # just test that the fixture builds the fragment h5
    # this also helps to make the timings more interpretable
    pass


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
            or (f1.gc - f2.gc < 1e-2)
        )
    )


def assert_fragments_identical(fs1, fs2):
    assert all(fragments_eq(f1, f2) for f1, f2 in zip(fs1, fs2))


def test_read_all(small_h5_path, bam_path, fasta_file_path):
    fmh5 = FragmentsH5(small_h5_path, "r")
    h5_fragments = list(fmh5.fetch(return_gc=True))
    bam_fragments = list(
        bam_to_fragments(
            bam_path, 'chr6', max_tlen=fmh5.max_fragment_length, fasta_file=fasta_file_path
        )
    )
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


def test_include_duplicates(duplicates_bam_path, fasta_file_path):
    """Test that include_duplicates parameter correctly includes/excludes duplicate-marked fragments.
    
    TODO: Investigate why this test hangs when num_processes > 1. The multiprocessing
    may be causing deadlocks or issues with the small test BAM file. Currently using
    num_processes=1 as a workaround.
    """
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
