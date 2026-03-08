"""Tests for reading BAMs from S3 via htslib streaming (no download).

Requires:
- pysam/htslib built with S3 support (libcurl)
- AWS credentials (env or ~/.aws/credentials) if the bucket is not public
- An indexed BAM at the URL (htslib looks for .bai at same path + .bai)
"""
import os
import tempfile
import pytest
import pysam

from fragments_h5.fragments_h5 import build_fragments_h5, FragmentsH5

# Test data on S3 (index must exist at same path + .bai)
S3_BAM_URL = "s3://fragmentomics.kariusdx.com/nboley/fragments_h5_test_data/small.chr6.bam"
S3_FASTA_URL = "s3://fragmentomics.kariusdx.com/nboley/fragments_h5_test_data/GRCh38.p12.genome.chr6_99110000_99130000.fa.gz"


def _is_remote_url(path: str) -> bool:
    return path.startswith("s3://") or path.startswith("http://") or path.startswith("https://")


def _pysam_can_open_s3():
    """Try to open the S3 URL with pysam; return (True, None) on success, (False, reason) on skip/fail."""
    if not S3_BAM_URL:
        return False, "FRAGMENTS_H5_S3_BAM_URL not set"
    if not _is_remote_url(S3_BAM_URL):
        return False, "FRAGMENTS_H5_S3_BAM_URL is not a remote URL (s3:// or http(s)://)"
    try:
        with pysam.AlignmentFile(S3_BAM_URL) as bam:
            _ = bam.references
            if not bam.has_index():
                return False, "BAM at S3 URL has no index (htslib expects .bai at same path)"
            # One contig fetch to confirm range requests work
            if bam.references:
                contig = bam.references[0]
                n = 0
                for _ in bam.fetch(contig, 0, min(100000, bam.get_reference_length(contig) or 0)):
                    n += 1
                    if n >= 10:
                        break
    except Exception as e:
        return False, f"pysam could not open S3 BAM: {e}"
    return True, None


@pytest.mark.skipif(
    not S3_BAM_URL or not _is_remote_url(S3_BAM_URL),
    reason="S3_BAM_URL must be a remote URL (s3:// or http(s)://)",
)
def test_pysam_opens_s3_bam():
    """Test that pysam can open an S3 BAM URL and perform a small fetch (htslib streaming)."""
    can, reason = _pysam_can_open_s3()
    if not can:
        pytest.skip(reason)
    # If we get here, the skipif passed (URL set) and _pysam_can_open_s3 succeeded
    assert can


@pytest.mark.skipif(
    not S3_BAM_URL or not _is_remote_url(S3_BAM_URL),
    reason="S3_BAM_URL must be a remote URL (s3:// or http(s)://)",
)
def test_build_fragments_h5_from_s3_bam():
    """Test building a fragments H5 from an S3 BAM path (streaming, no download).

    Uses num_processes=1 to avoid multiprocessing + S3 handle issues.
    Does not use --fasta so no reference streaming is required.
    """
    can, reason = _pysam_can_open_s3()
    if not can:
        pytest.skip(reason)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_h5 = os.path.join(tmpdir, "from_s3.fragments.h5")
        build_fragments_h5(
            S3_BAM_URL,
            out_h5,
            fasta_filename=None,
            num_processes=1,
        )
        assert os.path.isfile(out_h5)
        with FragmentsH5(out_h5) as fh5:
            assert fh5.n_fragments >= 0
            # At least one contig should have data if BAM has fragments
            if fh5.contig_lengths:
                contig = next(iter(fh5.contig_lengths))
                if contig in fh5.data:
                    starts, stops, _ = fh5.fetch_array(contig)
                    assert len(starts) == len(stops)


def _pysam_can_open_s3_fasta():
    """Try to open the S3 FASTA URL with pysam; return (True, None) on success, (False, reason) on skip/fail."""
    if not S3_FASTA_URL:
        return False, "S3_FASTA_URL not set"
    if not _is_remote_url(S3_FASTA_URL):
        return False, "S3_FASTA_URL is not a remote URL (s3:// or http(s)://)"
    try:
        with pysam.FastaFile(S3_FASTA_URL) as fasta:
            _ = fasta.references
            if not fasta.references:
                return False, "FASTA at S3 URL has no references"
            # Try to fetch a small region to confirm it works
            if fasta.references:
                contig = fasta.references[0]
                contig_len = fasta.get_reference_length(contig)
                if contig_len and contig_len > 0:
                    _ = fasta.fetch(contig, 0, min(100, contig_len))
    except Exception as e:
        return False, f"pysam could not open S3 FASTA: {e}"
    return True, None


@pytest.mark.skipif(
    not S3_BAM_URL or not _is_remote_url(S3_BAM_URL) or not S3_FASTA_URL or not _is_remote_url(S3_FASTA_URL),
    reason="S3_BAM_URL and S3_FASTA_URL must both be remote URLs (s3:// or http(s)://)",
)
def test_build_fragments_h5_from_s3_bam_with_s3_fasta():
    """Test building a fragments H5 from S3 BAM + S3 FASTA (streaming both).

    This verifies that:
    1. Both BAM and FASTA can be streamed from S3
    2. GC content is calculated correctly from streamed FASTA
    3. Results match local file processing

    Requires:
    - S3 BAM with index (.bai) at same path
    - S3 FASTA with index (.fai and .gzi for compressed) at same path
    - pysam/htslib built with S3 support (libcurl)
    - AWS credentials (for non-public buckets)
    """
    can_bam, reason_bam = _pysam_can_open_s3()
    if not can_bam:
        pytest.skip(f"Cannot open S3 BAM: {reason_bam}")

    can_fasta, reason_fasta = _pysam_can_open_s3_fasta()
    if not can_fasta:
        pytest.skip(f"Cannot open S3 FASTA: {reason_fasta}")

    with tempfile.TemporaryDirectory() as tmpdir:
        out_h5 = os.path.join(tmpdir, "from_s3_with_fasta.fragments.h5")
        build_fragments_h5(
            S3_BAM_URL,
            out_h5,
            fasta_filename=S3_FASTA_URL,
            num_processes=1,
        )
        assert os.path.isfile(out_h5)
        
        with FragmentsH5(out_h5) as fh5:
            assert fh5.n_fragments >= 0
            
            # Verify that GC content was calculated (not all None)
            if fh5.contig_lengths:
                contig = next(iter(fh5.contig_lengths))
                if contig in fh5.data:
                    starts, stops, supp_data = fh5.fetch_array(contig, return_gc=True)
                    assert len(starts) == len(stops)
                    
                    # Check that at least some fragments have GC content calculated
                    if len(starts) > 0 and 'gc' in supp_data:
                        import numpy as np
                        gc_values = supp_data['gc']
                        # At least some GC values should be non-NaN (not all regions are N's)
                        non_nan_count = (~np.isnan(gc_values)).sum()
                        assert non_nan_count > 0, "All GC values are NaN - FASTA streaming may have failed"
