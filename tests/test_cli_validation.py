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


@pytest.fixture(scope="module")
def bam_path():
    """Real BAM fixture (also used by tests/test_fragments_h5.py). The SE
    length-filter gate in build_sub_fragments_h5 (`if single_end and
    se_max_fragment_length is not None and ...`) is only reachable with BAM
    input: TSV input always neutralizes single_end and se_max_fragment_length
    before the gate is ever consulted (see TestTsvNeutralization)."""
    return os.path.join(DATA_DIR, "small.chr6.bam")


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
        """--min-mapq should be warned and neutralized for TSV input.

        Note: the n_fragments == 3 assertion below cannot fail regardless of
        neutralization. tsv_to_fragments accepts **kwargs and never reads
        min_mapq, so TSV input is structurally incapable of MAPQ filtering
        whether or not the neutralizer runs — it's a baseline sanity check
        (the build still produces output), not evidence of neutralization.
        The caplog warning assertion below is the only part of this test with
        teeth: it is the sole check that the neutralizer's warning path
        actually fires.
        """
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

    def test_single_end_and_length_neutralized_together_for_tsv(
        self, tsv_fragment_file, fasta_file_path
    ):
        """This is a duplicate of test_single_end_neutralized_for_tsv with a more
        extreme se_max_fragment_length value (1 instead of 50). It exercises TSV
        NEUTRALIZATION, not the SE length-filter gate in build_sub_fragments_h5
        (`if single_end and se_max_fragment_length is not None and ...`): for
        TSV input the neutralizer sets single_end=False and
        se_max_fragment_length=None before the gate is ever consulted, so this
        test would pass even if that gate were deleted entirely. Renamed from
        the former (misleadingly named)
        TestSeFilterGating.test_se_filter_not_applied_when_not_single_end, which
        claimed to test the gate but actually only exercised this same
        neutralization path. See TestSeFilterGating below for real gate tests."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "out.h5")
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


# ── SE filter gating (BAM input only — see bam_path fixture) ──

class TestSeFilterGating:
    """Verify the SE length-filter gate in build_sub_fragments_h5:

        if single_end and se_max_fragment_length is not None and (fragment.stop - fragment.start) > se_max_fragment_length:

    actually requires single_end=True before applying the SE length filter.
    This is only testable with BAM input: for TSV input the neutralizer in
    build_fragments_h5 always clears single_end and se_max_fragment_length
    before the gate is reached (see TestTsvNeutralization), so a TSV-based test
    can never exercise this line. Prior to this test class the gate had zero
    test coverage and could be deleted with no test failures.
    """

    def test_gate_blocks_filter_when_not_single_end(self, bam_path, fasta_file_path):
        """single_end=False + se_max_fragment_length=1 must NOT filter anything.

        Every real fragment in the BAM fixture is far longer than 1bp, so if the
        `single_end and` conjunct were removed from the SE length-filter gate
        in build_sub_fragments_h5, the length filter would apply regardless of
        single_end and n_fragments would collapse to 0. Detects: deleting
        `single_end and` from that gate.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            baseline = os.path.join(tmpdir, "baseline.h5")
            build_fragments_h5(
                bam_path, baseline, fasta_filename=fasta_file_path, num_processes=1,
            )
            with FragmentsH5(baseline) as fh5:
                baseline_count = fh5.n_fragments
            assert baseline_count > 0, "test BAM should have at least one fragment"

            output = os.path.join(tmpdir, "out.h5")
            build_fragments_h5(
                bam_path,
                output,
                fasta_filename=fasta_file_path,
                single_end=False,
                se_max_fragment_length=1,
                num_processes=1,
            )
            with FragmentsH5(output) as fh5:
                assert fh5.n_fragments == baseline_count

    def test_gate_applies_filter_when_single_end(self, bam_path, fasta_file_path):
        """single_end=True + se_max_fragment_length=1 must filter out every
        fragment, since every real fragment exceeds 1bp. This is the first test
        anywhere that the SE length filter actually removes fragments on its
        primary (BAM, single_end=True) use case. Removing the `single_end and`
        conjunct from the SE length-filter gate in build_sub_fragments_h5 does
        not change this test's outcome (the filter still applies either way),
        so this test alone does not detect that mutation — it is the
        companion test above that does.

        build_fragments_h5 raises ValueError rather than writing an empty h5
        when zero fragments survive filtering (see the "No fragments were
        extracted" check in build_fragments_h5); that ValueError is the
        observable signal here that every fragment was dropped.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "out.h5")
            with pytest.raises(ValueError, match="No fragments were extracted"):
                build_fragments_h5(
                    bam_path,
                    output,
                    fasta_filename=fasta_file_path,
                    single_end=True,
                    se_max_fragment_length=1,
                    num_processes=1,
                )


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
