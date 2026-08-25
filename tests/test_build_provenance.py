"""Tests for build provenance (item 5): _build_argv and _build_version attributes."""
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile
import types
import unittest.mock

import h5py
import pysam
import pytest

from fragments_h5.fragments_h5 import (
    build_fragments_h5, FragmentsH5, _resolve_build_code_revision,
)
import fragments_h5.fragments_h5 as fh5_module
import fragments_h5.main as main_module


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
        # Delete provenance attrs to simulate a pre-2.12.0 file.
        # _build_argv is genuinely conditional here: this file was built without
        # build_argv (see call above), so it was never written. _build_version is
        # always present for an installed build, so it is deleted unconditionally.
        with h5py.File(h5_path, "r+") as f:
            if "_build_argv" in f.attrs:
                del f.attrs["_build_argv"]
            del f.attrs["_build_version"]
            if "_build_code_revision" in f.attrs:
                del f.attrs["_build_code_revision"]

        with FragmentsH5(h5_path) as fh5:
            assert fh5.build_argv is None
            assert fh5.build_version is None
            assert fh5.build_code_revision is None


def test_cli_records_build_argv():
    """End-to-end: build-fragments-h5's main() should record the real
    sys.argv it was invoked with into `_build_argv`.

    Previously this pass-through (main.py:177, `build_argv=sys.argv`) was
    only verified by inference from library-level tests that hand-construct
    a fake_argv and pass it directly to build_fragments_h5 -- main() itself
    was never driven through a full build.

    Mutation this detects: removing `build_argv=sys.argv` (or passing None)
    from the build_fragments_h5 call at main.py:177 makes `_build_argv`
    absent from the output file, and the `"_build_argv" in f.attrs`
    assertion below goes red.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "test.bam")
        _make_simple_bam(bam)
        h5_path = os.path.join(tmpdir, "out.h5")
        argv = [
            "build-fragments-h5", bam, h5_path,
            "--single-end", "--se-max-fragment-length", "1000",
            "--num-processes", "1",
        ]
        with unittest.mock.patch.object(sys, "argv", argv):
            main_module.main()

        with h5py.File(h5_path, "r") as f:
            assert "_build_argv" in f.attrs
            assert json.loads(f.attrs["_build_argv"]) == argv


def test_package_not_found_omits_build_version(monkeypatch):
    """When importlib.metadata.version raises PackageNotFoundError (e.g. an
    uninstalled tree), `_build_version` must be omitted from the output
    rather than written as an error string or crashing the build.

    Previously only the *result* of this branch was simulated -- by
    deleting `_build_version` after a successful build in
    test_file_without_provenance_opens -- but the `except
    PackageNotFoundError: pass` branch (fragments_h5.py:1184-1185) itself
    was never executed by any test.

    Mutation this detects: removing the try/except around
    `importlib.metadata.version(...)` at fragments_h5.py:1182-1185 makes
    PackageNotFoundError propagate out of build_fragments_h5 instead of
    being swallowed; this test's build_fragments_h5() call would raise
    instead of completing.
    """
    def _raise(*args, **kwargs):
        raise importlib.metadata.PackageNotFoundError("fragments-h5")

    monkeypatch.setattr(importlib.metadata, "version", _raise)

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
            assert "_build_version" not in f.attrs


# --- _build_code_revision tests ---


def test_revision_prefix_contract():
    """Build a file; assert _build_code_revision is either absent or a str
    starting with one of the four valid prefixes.

    Mutation this detects: changing a prefix string (e.g. "git:" -> "GIT:")
    in the resolver makes the assertion fail.
    """
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
            rev = f.attrs.get("_build_code_revision")
            if rev is not None:
                valid_prefixes = ("git:", "baked:", "dist:", "dist-editable:")
                assert isinstance(rev, str)
                assert any(rev.startswith(p) for p in valid_prefixes), \
                    f"unexpected prefix in {rev!r}"


def test_revision_git_branch():
    """When running from a git-tracked file, the resolver returns 'git:...'.

    Skip on the gate itself (ls-files --error-unmatch exiting nonzero),
    not on .git presence — they diverge where the gate matters.

    Mutation this detects: removing the "git:" prefix from the return
    statement makes startswith("git:") fail.
    """
    pkg_dir = os.path.dirname(fh5_module.__file__)
    gate = subprocess.run(
        ["git", "-C", pkg_dir, "ls-files", "--error-unmatch",
         os.path.basename(fh5_module.__file__)],
        capture_output=True, timeout=5, text=True,
    )
    if gate.returncode != 0:
        pytest.skip("package file not tracked by git")

    result = _resolve_build_code_revision()
    assert result is not None
    assert result.startswith("git:"), f"expected 'git:' prefix, got {result!r}"


def test_revision_gate_rejects_untracked_copy(monkeypatch):
    """A package file copied into a scratch repo but NOT git-added must
    NOT get a 'git:' prefix. Tests the ls-files --error-unmatch gate.

    Mutation this detects: removing the gate (the ls-files --error-unmatch
    call and its returncode check) would let the untracked copy inherit
    a 'git:' label from the enclosing repo's describe output.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run(["git", "init", tmpdir], capture_output=True, check=True)
        subprocess.run(
            ["git", "-C", tmpdir, "commit", "--allow-empty", "-m", "init"],
            capture_output=True, check=True,
        )
        pkg_src = fh5_module.__file__
        dest = os.path.join(tmpdir, os.path.basename(pkg_src))
        shutil.copy2(pkg_src, dest)

        monkeypatch.setattr(fh5_module, "__file__", dest)
        result = _resolve_build_code_revision()

    assert result is None or not result.startswith("git:"), \
        f"untracked copy should not get git: prefix, got {result!r}"


def test_revision_baked_branch(monkeypatch):
    """When git is unavailable but _build_revision.py exists, the resolver
    returns 'baked:<value>'.

    This is the container path: .dockerignore excludes .git, the runtime
    image has no git binary, and the Dockerfile bakes a revision at build
    time.

    Mutation this detects: removing the `from ._build_revision import
    REVISION` try/except block makes the result fall through to
    dist-editable (or dist) instead of returning 'baked:vFAKE'.
    """
    def _no_git(*args, **kwargs):
        raise FileNotFoundError("git")
    monkeypatch.setattr(fh5_module.subprocess, "run", _no_git)

    fake_mod = types.ModuleType("fragments_h5._build_revision")
    fake_mod.REVISION = "vFAKE"
    monkeypatch.setitem(sys.modules, "fragments_h5._build_revision", fake_mod)

    result = _resolve_build_code_revision()
    assert result == "baked:vFAKE"


def test_revision_dist_branch(monkeypatch):
    """When git is unavailable, no baked revision, and the install is NOT
    editable, the resolver returns 'dist:<version>'.

    Must assert startswith("dist:") and NOT startswith("dist-editable:") —
    a loose startswith("dist") would conflate the two branches.

    Mutation this detects: removing the importlib.metadata.version fallback
    in step 3 makes the result None instead of 'dist:...'.
    """
    def _no_git(*args, **kwargs):
        raise FileNotFoundError("git")
    monkeypatch.setattr(fh5_module.subprocess, "run", _no_git)
    monkeypatch.delitem(sys.modules, "fragments_h5._build_revision", raising=False)

    # Force the editable check to return non-editable.
    # importlib.metadata.version() calls distribution().version internally,
    # so the wrapper must delegate all attribute access except read_text.
    real_distribution = importlib.metadata.distribution

    class _NonEditableDist:
        def __init__(self, dist):
            self._dist = dist
        def read_text(self, name):
            if name == "direct_url.json":
                return None
            return self._dist.read_text(name)
        def __getattr__(self, name):
            return getattr(self._dist, name)

    def _fake_dist(name):
        return _NonEditableDist(real_distribution(name))
    monkeypatch.setattr(importlib.metadata, "distribution", _fake_dist)

    result = _resolve_build_code_revision()
    assert result is not None
    assert result.startswith("dist:"), f"expected 'dist:' prefix, got {result!r}"
    assert not result.startswith("dist-editable:"), \
        f"should be 'dist:', not 'dist-editable:', got {result!r}"


def test_revision_dist_editable_branch(monkeypatch):
    """When git is unavailable and the install IS editable, the resolver
    returns 'dist-editable:<version>'.

    This machine has an editable install, so with git disabled the natural
    fallback is dist-editable.

    Mutation this detects: removing the editable check (the direct_url.json
    probe) makes the result 'dist:...' instead of 'dist-editable:...'.
    """
    def _no_git(*args, **kwargs):
        raise FileNotFoundError("git")
    monkeypatch.setattr(fh5_module.subprocess, "run", _no_git)
    monkeypatch.delitem(sys.modules, "fragments_h5._build_revision", raising=False)

    result = _resolve_build_code_revision()
    assert result is not None
    assert result.startswith("dist-editable:"), \
        f"expected 'dist-editable:' prefix, got {result!r}"


def test_revision_total_failure(monkeypatch):
    """When git, baked, and importlib all fail, the resolver returns None
    and the attribute is omitted.

    Mutation this detects: returning a sentinel like "unknown" instead of
    None would fail the `is None` assertion.
    """
    def _no_git(*args, **kwargs):
        raise FileNotFoundError("git")
    monkeypatch.setattr(fh5_module.subprocess, "run", _no_git)
    monkeypatch.delitem(sys.modules, "fragments_h5._build_revision", raising=False)

    def _raise(*args, **kwargs):
        raise importlib.metadata.PackageNotFoundError("fragments-h5")
    monkeypatch.setattr(importlib.metadata, "version", _raise)

    result = _resolve_build_code_revision()
    assert result is None


def test_revision_total_failure_omits_attribute(monkeypatch):
    """End-to-end: when the resolver returns None, the h5 file must NOT
    contain _build_code_revision, and the accessor must be None.

    Mutation this detects: writing a sentinel value when the resolver
    returns None (e.g. f.attrs["_build_code_revision"] = "unknown")
    would fail the `not in f.attrs` assertion.
    """
    def _no_git(*args, **kwargs):
        raise FileNotFoundError("git")
    monkeypatch.setattr(fh5_module.subprocess, "run", _no_git)
    monkeypatch.delitem(sys.modules, "fragments_h5._build_revision", raising=False)

    def _raise(*args, **kwargs):
        raise importlib.metadata.PackageNotFoundError("fragments-h5")
    monkeypatch.setattr(importlib.metadata, "version", _raise)

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
            assert "_build_code_revision" not in f.attrs
        with FragmentsH5(h5_path) as fh5:
            assert fh5.build_code_revision is None


def test_revision_old_file_degradation():
    """Files predating _build_code_revision open with build_code_revision
    == None.

    Mutation this detects: using attrs["_build_code_revision"] instead of
    attrs.get("_build_code_revision") in __init__ would KeyError on files
    lacking the attribute.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "test.bam")
        _make_simple_bam(bam)
        h5_path = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5_path,
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
        )
        with h5py.File(h5_path, "r+") as f:
            if "_build_code_revision" in f.attrs:
                del f.attrs["_build_code_revision"]

        with FragmentsH5(h5_path) as fh5:
            assert fh5.build_code_revision is None
