"""Tests for CLI argument validation and library-level TSV flag neutralization."""

import logging
import os
import subprocess
import sys
import tempfile

import pytest

from fragments_h5.fragments_h5 import build_fragments_h5, FragmentsH5
from fragments_h5.fragment import is_fragment_file

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")
# pytest already runs under an interpreter that has fragments_h5 importable,
# so reuse it rather than hardcoding a path into one developer's conda env.
PYTHON = sys.executable


def _run_cli(args, expect_fail=False):
    """Run build-fragments-h5 via subprocess and return (returncode, stderr)."""
    cmd = [PYTHON, "-m", "fragments_h5.main"] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    if expect_fail:
        assert result.returncode != 0, f"Expected failure but got exit 0; stderr: {result.stderr}"
    return result.returncode, result.stderr


# ── Range validation ──

class TestRangeValidation:
    def test_se_max_fragment_length_zero_rejected(self):
        rc, stderr = _run_cli(
            ["dummy.bam", "out.h5", "--single-end", "--se-max-fragment-length", "0"],
            expect_fail=True,
        )
        assert rc == 2  # parser.error gives exit 2
        assert "--se-max-fragment-length must be between 1 and 65535" in stderr

    def test_se_max_fragment_length_negative_rejected(self):
        rc, stderr = _run_cli(
            ["dummy.bam", "out.h5", "--single-end", "--se-max-fragment-length", "-1"],
            expect_fail=True,
        )
        assert rc == 2
        assert "--se-max-fragment-length must be between 1 and 65535" in stderr

    def test_se_max_fragment_length_too_large_rejected(self):
        rc, stderr = _run_cli(
            ["dummy.bam", "out.h5", "--single-end", "--se-max-fragment-length", "65536"],
            expect_fail=True,
        )
        assert rc == 2
        assert "--se-max-fragment-length must be between 1 and 65535" in stderr

    def test_min_mapq_negative_rejected(self):
        rc, stderr = _run_cli(
            ["dummy.bam", "out.h5", "--min-mapq", "-1"],
            expect_fail=True,
        )
        assert rc == 2
        assert "--min-mapq must be >= 0" in stderr


# ── Argument consistency (BAM paths) ──

class TestBamArgConsistency:
    def test_single_end_requires_se_max_fragment_length(self):
        rc, stderr = _run_cli(
            ["dummy.bam", "out.h5", "--single-end"],
            expect_fail=True,
        )
        assert rc == 2
        assert "--se-max-fragment-length is required" in stderr

    def test_se_max_fragment_length_requires_single_end(self):
        rc, stderr = _run_cli(
            ["dummy.bam", "out.h5", "--se-max-fragment-length", "120"],
            expect_fail=True,
        )
        assert rc == 2
        assert "--se-max-fragment-length can only be used with --single-end" in stderr


# ── TSV input skips BAM-only checks ──

class TestTsvSkipsRequiredCheck:
    """For TSV input, --single-end and --se-max-fragment-length are ignored by
    the library, so the CLI must not force callers to supply them."""

    def test_tsv_single_end_without_se_max_accepted(self):
        """--single-end on TSV input should not require --se-max-fragment-length.
        The command will fail later (missing --fasta), but NOT with the
        required-check error."""
        rc, stderr = _run_cli(
            ["dummy.tsv.gz", "out.h5", "--single-end"],
            expect_fail=True,
        )
        # Should NOT fail with the --se-max-fragment-length required error
        assert "--se-max-fragment-length is required" not in stderr

    def test_tsv_se_max_without_single_end_accepted(self):
        """--se-max-fragment-length on TSV input should not require --single-end."""
        rc, stderr = _run_cli(
            ["dummy.tsv.gz", "out.h5", "--se-max-fragment-length", "120"],
            expect_fail=True,
        )
        assert "--se-max-fragment-length can only be used with --single-end" not in stderr


# ── Library-level TSV neutralization ──

@pytest.fixture(scope="module")
def tsv_fragment_file():
    """Create a minimal bgzipped, tabix-indexed TSV fragment file for testing."""
    import pysam

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write plain TSV first, then bgzip via pysam
        plain_path = os.path.join(tmpdir, "test_frags.tsv")
        with open(plain_path, "w") as f:
            # BED-like: chrom, start, stop, name, score, strand
            # Coordinates within the FASTA test region (chr6:99110000-99130000)
            f.write("chr6\t99115000\t99115200\tfrag1\t0\t+\n")
            f.write("chr6\t99115300\t99115500\tfrag2\t0\t+\n")
            f.write("chr6\t99115600\t99115700\tfrag3\t0\t-\n")

        # bgzip (produces test_frags.tsv.gz and removes plain file)
        pysam.tabix_compress(plain_path, plain_path + ".gz", force=True)
        bgz_path = plain_path + ".gz"

        # tabix index
        pysam.tabix_index(bgz_path, preset="bed", force=True)

        yield bgz_path


@pytest.fixture(scope="module")
def fasta_file_path():
    return os.path.join(DATA_DIR, "GRCh38.p12.genome.chr6_99110000_99130000.fa.gz")


class TestTsvNeutralization:
    """Test that build_fragments_h5 warns and neutralizes BAM-only flags for TSV input."""

    def test_single_end_neutralized_for_tsv(self, tsv_fragment_file, fasta_file_path, caplog):
        """--single-end + --se-max-fragment-length should be warned and neutralized.
        No fragments should be dropped by the SE length filter."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "out.h5")
            with caplog.at_level(logging.WARNING, logger="fragments_h5.fragments_h5"):
                build_fragments_h5(
                    tsv_fragment_file,
                    output,
                    fasta_filename=fasta_file_path,
                    single_end=True,
                    se_max_fragment_length=50,  # smaller than any fragment
                    num_processes=1,
                )
            # Should have warned about --single-end
            assert any("--single-end flag is meaningless" in m for m in caplog.messages)

            # All 3 fragments should be present (the SE filter was neutralized,
            # not applied — without neutralization, all would be dropped since
            # every fragment is >50bp)
            with FragmentsH5(output) as fh5:
                assert fh5.n_fragments == 3

    def test_min_mapq_neutralized_for_tsv(self, tsv_fragment_file, fasta_file_path, caplog):
        """--min-mapq should be warned and neutralized for TSV input."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "out.h5")
            with caplog.at_level(logging.WARNING, logger="fragments_h5.fragments_h5"):
                build_fragments_h5(
                    tsv_fragment_file,
                    output,
                    fasta_filename=fasta_file_path,
                    min_mapq=30,
                    num_processes=1,
                )
            assert any("--min-mapq flag is meaningless" in m for m in caplog.messages)

            # All fragments should still be present
            with FragmentsH5(output) as fh5:
                assert fh5.n_fragments == 3

    def test_tsv_without_flags_builds_normally(self, tsv_fragment_file, fasta_file_path):
        """TSV build without any BAM-only flags should work without warnings."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "out.h5")
            build_fragments_h5(
                tsv_fragment_file,
                output,
                fasta_filename=fasta_file_path,
                num_processes=1,
            )
            with FragmentsH5(output) as fh5:
                assert fh5.n_fragments == 3


# ── SE filter gating ──

class TestSeFilterGating:
    """Verify the SE length filter at build_sub_fragments_h5 is gated on single_end."""

    def test_se_filter_not_applied_when_not_single_end(self, tsv_fragment_file, fasta_file_path):
        """Even if se_max_fragment_length leaks through somehow, the filter
        should not apply unless single_end is True. The TSV neutralization
        already sets single_end=False, so this is a defense-in-depth test."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "out.h5")
            # Call with single_end=True + tiny se_max_fragment_length on TSV.
            # The library neutralizes these for TSV, so all fragments survive.
            build_fragments_h5(
                tsv_fragment_file,
                output,
                fasta_filename=fasta_file_path,
                single_end=True,
                se_max_fragment_length=1,  # Would drop everything if not neutralized
                num_processes=1,
            )
            with FragmentsH5(output) as fh5:
                assert fh5.n_fragments == 3


# ── is_fragment_file ──

class TestIsFragmentFile:
    def test_tsv_gz(self):
        assert is_fragment_file("frags.tsv.gz") is True

    def test_bed_gz(self):
        assert is_fragment_file("frags.bed.gz") is True

    def test_bam(self):
        assert is_fragment_file("input.bam") is False

    def test_plain_tsv(self):
        assert is_fragment_file("frags.tsv") is False
