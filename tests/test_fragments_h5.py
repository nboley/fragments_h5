import os
import pytest
import subprocess
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
