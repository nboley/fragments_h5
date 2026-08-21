"""Tests for build provenance (item 5): _build_argv and _build_version attributes."""
import importlib.metadata
import json
import os
import tempfile

import h5py
import pysam
import pytest

from fragments_h5.fragments_h5 import build_fragments_h5, FragmentsH5


def _make_simple_bam(path):
    """Create a minimal single-end BAM with a few reads for provenance testing."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 10000}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as outf:
        for i in range(3):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i}"
            a.reference_id = 0
            a.reference_start = i * 100
            a.cigarstring = "10M"
            a.mapping_quality = 60
            a.query_sequence = "A" * 10
            a.query_qualities = pysam.qualitystring_to_array("I" * 10)
            a.flag = 0
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            outf.write(a)
    pysam.sort("-o", path, path)
    pysam.index(path)


def test_provenance_absent_for_library_caller():
    """Library caller (no build_argv) should have _build_version but not _build_argv."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "test.bam")
        _make_simple_bam(bam)
        h5_path = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5_path,
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
        )
        with h5py.File(h5_path, "r") as f:
            assert "_build_argv" not in f.attrs
            assert f.attrs["_build_version"] == importlib.metadata.version("fragments-h5")


def test_provenance_recorded_when_argv_passed():
    """When build_argv is passed, it round-trips as JSON."""
    fake_argv = ["build-fragments-h5", "in.bam", "out.h5", "--single-end", "--se-max-fragment-length", "500"]
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "test.bam")
        _make_simple_bam(bam)
        h5_path = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5_path,
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
            build_argv=fake_argv,
        )
        with h5py.File(h5_path, "r") as f:
            assert "_build_argv" in f.attrs
            assert json.loads(f.attrs["_build_argv"]) == fake_argv


def test_provenance_accessors():
    """FragmentsH5.build_argv and .build_version return parsed values."""
    fake_argv = ["build-fragments-h5", "input.bam", "output.h5"]
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "test.bam")
        _make_simple_bam(bam)
        h5_path = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5_path,
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
            build_argv=fake_argv,
        )
        with FragmentsH5(h5_path) as fh5:
            assert fh5.build_argv == fake_argv
            assert fh5.build_version == importlib.metadata.version("fragments-h5")


def test_file_without_provenance_opens():
    """Files without provenance attrs (pre-2.12.0) should open with both accessors as None."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "test.bam")
        _make_simple_bam(bam)
        h5_path = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5_path,
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
        )
        # Delete provenance attrs to simulate a pre-2.12.0 file
        with h5py.File(h5_path, "r+") as f:
            if "_build_argv" in f.attrs:
                del f.attrs["_build_argv"]
            if "_build_version" in f.attrs:
                del f.attrs["_build_version"]

        with FragmentsH5(h5_path) as fh5:
            assert fh5.build_argv is None
            assert fh5.build_version is None
