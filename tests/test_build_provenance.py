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
    """Library caller (no build_argv) should have neither _build_argv nor
    _build_version. _build_version is no longer written for newly built
    files (2026-08-25 decision) -- _build_code_revision is the authoritative
    code-identity field now; this asserts the new-file contract directly.

    Mutation this detects: reinstating the removed
    `f.attrs["_build_version"] = importlib.metadata.version(...)` write
    makes the `"_build_version" not in f.attrs` assertion go red.
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
            assert "_build_argv" not in f.attrs
            assert "_build_version" not in f.attrs


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
    """FragmentsH5.build_argv returns the parsed value. build_version is
    None for a newly built file, since it is no longer written
    (2026-08-25 decision) -- see test_old_file_build_version_still_exposed
    for the backward-compatibility accessor contract on files that do have
    the attribute.
    """
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
            assert fh5.build_version is None


def test_old_file_build_version_still_exposed():
    """A file carrying _build_version (as written by 2.12.0/2.12.1, before
    the 2026-08-25 decision to stop writing it) must still expose that value
    through the build_version accessor -- new files simply never have the
    attribute to begin with, but old ones that do must keep working.

    Mutation this detects: removing (or hardcoding to None) the
    `.attrs.get("_build_version")` read in FragmentsH5.__init__ makes the
    equality assertion below fail.
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
        # Simulate a file written by 2.12.0/2.12.1, which wrote this
        # attribute unconditionally (when determinable).
        with h5py.File(h5_path, "r+") as f:
            f.attrs["_build_version"] = "2.11.0"
        with FragmentsH5(h5_path) as fh5:
            assert fh5.build_version == "2.11.0"


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
        # also never written for a new file (2026-08-25 decision), so both
        # deletions below are conditional/no-ops in practice; kept unconditional-
        # style guards so this test still passes if that decision is ever
        # reversed.
        with h5py.File(h5_path, "r+") as f:
            if "_build_argv" in f.attrs:
                del f.attrs["_build_argv"]
            if "_build_version" in f.attrs:
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


# NOTE (2026-08-25): test_package_not_found_omits_build_version previously
# lived here. It existed solely to cover the `except PackageNotFoundError:
# pass` branch around the `_build_version` write, which has been deleted
# entirely (not written at all, unconditionally) -- there is no longer a
# try/except at that site for PackageNotFoundError to exercise. It was not
# replaced 1:1: the broader property it also touched on -- that
# PackageNotFoundError anywhere in build_fragments_h5 must not crash the
# build -- is already independently covered by
# test_revision_total_failure_omits_attribute below, which monkeypatches the
# same importlib.metadata.version and drives a full build_fragments_h5()
# call to completion. test_provenance_absent_for_library_caller above covers
# the new contract (newly built files never have _build_version, regardless
# of import state).


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


def test_revision_untracked_suffix(monkeypatch):
    """An untracked file anywhere in the repo (the tracked package file
    itself unmodified) makes the result carry the '-untracked' suffix.

    Mutation this detects: removing the ls-files --others --exclude-standard
    check (or the `revision += "-untracked"` line) makes '-untracked' absent
    from the result, and this test's `in result` assertion goes red.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run(["git", "init", tmpdir], capture_output=True, check=True)
        pkg_src = fh5_module.__file__
        dest = os.path.join(tmpdir, os.path.basename(pkg_src))
        shutil.copy2(pkg_src, dest)
        subprocess.run(["git", "-C", tmpdir, "add", os.path.basename(dest)],
                        capture_output=True, check=True)
        subprocess.run(["git", "-C", tmpdir, "commit", "-m", "init"],
                        capture_output=True, check=True)

        # An unrelated untracked file, unrelated to the package file itself.
        with open(os.path.join(tmpdir, "scratch.txt"), "w") as f:
            f.write("untracked\n")

        monkeypatch.setattr(fh5_module, "__file__", dest)
        result = _resolve_build_code_revision()

    assert result is not None
    assert result.startswith("git:")
    assert "-untracked" in result, f"expected '-untracked' in {result!r}"


def test_revision_dirty_without_untracked_suffix(monkeypatch):
    """A tracked-file modification with no untracked files present yields
    '-dirty' (from `git describe --dirty`) but must NOT also carry
    '-untracked' -- the two suffixes are genuinely distinct signals.

    Mutation this detects: keying the '-untracked' suffix off a signal that
    also fires for ordinary dirtiness (e.g. `git status --porcelain`, which
    reports untracked files with the same leading marker as modified ones)
    would make '-untracked' incorrectly appear here.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run(["git", "init", tmpdir], capture_output=True, check=True)
        pkg_src = fh5_module.__file__
        dest = os.path.join(tmpdir, os.path.basename(pkg_src))
        shutil.copy2(pkg_src, dest)
        subprocess.run(["git", "-C", tmpdir, "add", os.path.basename(dest)],
                        capture_output=True, check=True)
        subprocess.run(["git", "-C", tmpdir, "commit", "-m", "init"],
                        capture_output=True, check=True)

        # Modify the tracked package file itself -- no untracked files exist.
        with open(dest, "a") as f:
            f.write("\n# scratch modification\n")

        monkeypatch.setattr(fh5_module, "__file__", dest)
        result = _resolve_build_code_revision()

    assert result is not None
    assert result.startswith("git:")
    assert "-dirty" in result, f"expected '-dirty' in {result!r}"
    assert "-untracked" not in result, f"unexpected '-untracked' in {result!r}"


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

    The editable install is FAKED here, mirroring test_revision_dist_branch.
    An earlier version of this test relied on the real environment being an
    editable install, which asserts a property of the machine rather than of
    the code: importlib.metadata returns whichever *.dist-info / *.egg-info
    appears first on sys.path, and an egg-info carries no direct_url.json, so
    it reads as NON-editable. That made this test fail on an ordinary editable
    checkout of this repo.

    Mutation this detects: removing the editable check (the direct_url.json
    probe) makes the result 'dist:...' instead of 'dist-editable:...'.
    """
    def _no_git(*args, **kwargs):
        raise FileNotFoundError("git")
    monkeypatch.setattr(fh5_module.subprocess, "run", _no_git)
    monkeypatch.delitem(sys.modules, "fragments_h5._build_revision", raising=False)

    # Force the editable check to return editable.
    # importlib.metadata.version() calls distribution().version internally,
    # so the wrapper must delegate all attribute access except read_text.
    real_distribution = importlib.metadata.distribution

    class _EditableDist:
        def __init__(self, dist):
            self._dist = dist
        def read_text(self, name):
            if name == "direct_url.json":
                return json.dumps({
                    "url": "file:///fake/checkout",
                    "dir_info": {"editable": True},
                })
            return self._dist.read_text(name)
        def __getattr__(self, name):
            return getattr(self._dist, name)

    def _fake_dist(name):
        return _EditableDist(real_distribution(name))
    monkeypatch.setattr(importlib.metadata, "distribution", _fake_dist)

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


def test_revision_written_into_real_build(monkeypatch):
    """A real build_fragments_h5() call must actually write
    _build_code_revision into the output file's attrs -- not merely have
    an accessor that would return None for both a legitimate absence and a
    broken write.

    The resolver itself (_resolve_build_code_revision) is exercised in
    isolation by the test_revision_*_branch tests above; this test does not
    re-verify resolver logic. It monkeypatches the resolver to a fixed,
    non-rotting sentinel value instead of asserting on whatever the ambient
    git/dist environment happens to produce, because a real 'git:...'
    string changes at every commit (this project has been bitten by rotting
    hardcoded references three times) -- pinning the resolver's return value
    lets this test assert byte-for-byte equality against the write site
    only, which is the actual gap: everything upstream of the `f.attrs[...]
    = _code_rev` line is already covered elsewhere.

    Mutation this detects: replacing
    `f.attrs["_build_code_revision"] = _code_rev` with `pass` (leaving the
    `if _code_rev is not None:` guard intact) makes the attribute absent
    from a build where the resolver is known to have returned a non-None
    value, so `"_build_code_revision" in f.attrs` goes red. Before this
    test, that exact mutation left all 17 tests in this file green.
    """
    sentinel = "git:sentinel-does-not-rot-1234567"
    monkeypatch.setattr(fh5_module, "_resolve_build_code_revision", lambda: sentinel)

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
            assert "_build_code_revision" in f.attrs
            assert f.attrs["_build_code_revision"] == sentinel
        with FragmentsH5(h5_path) as fh5:
            assert fh5.build_code_revision == sentinel


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
