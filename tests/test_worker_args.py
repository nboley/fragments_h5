"""Tests for the SubBuildArgs dataclass and its derived accessor in
build_fragments_h5 (see docs/architecture/worker_args_refactor.md).

Both tests here target the code paths implicated in the released v2.10.1
defect: a stale positional read of the worker-args tuple,
`sum(a[3] - a[2] for a in args)`, computed `chunk_start - output_contig`
(int - str) and raised TypeError on every build. That accessor is now
`total_bases = sum(a.chunk_stop - a.chunk_start for a in args)`.
"""
import os
import tempfile

import pysam
import pytest

import fragments_h5.fragments_h5 as fh5_module
from fragments_h5.fragments_h5 import SubBuildArgs, build_fragments_h5


def _write_se_bam(path, contigs, reads):
    """contigs: list of (name, length). reads: list of dicts with name/contig/pos/cigar."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as outf:
        for r in reads:
            a = pysam.AlignedSegment()
            a.query_name = r["name"]
            a.reference_id = outf.get_tid(r["contig"])
            a.reference_start = r["pos"]
            a.cigarstring = r.get("cigar", "10M")
            a.mapping_quality = r.get("mapq", 60)
            seq_len = r.get("seq_len", 10)
            a.query_sequence = "A" * seq_len
            a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
            a.flag = r.get("flag", 0)
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            outf.write(a)
    pysam.sort("-o", path, path)
    pysam.index(path)


# ── Item 3: total_bases must pin a VALUE, not just avoid an exception ──

def test_total_bases_progress_total_matches_actual_chunk_sum(monkeypatch):
    """`total_bases = sum(a.chunk_stop - a.chunk_start for a in args)` (the
    direct successor of the positional expression that shipped a released
    TypeError on every build at v2.10.1) is pinned two ways:

    1. The value handed to tqdm's `total=` kwarg is intercepted, so the
       progress bar not crashing is no longer sufficient -- the actual
       number must be right.
    2. GENOMIC_CHUNK_SIZE is patched down so a contig splits into multiple
       chunks with chunk_start != 0 for later chunks. With chunk_start == 0
       (the untouched-default case, where a contig fits in one chunk),
       `contig_length - chunk_start` and `chunk_stop - chunk_start` are both
       equal to contig_length -- numerically indistinguishable, exactly the
       blind spot that let a wrong-field expression through undetected in
       the original tuple-index defect. Forcing chunk_start != 0 for some
       chunks closes that blind spot.

    A contig with zero mapped fragments (chrZ) is included and must
    contribute nothing to the total, since it never becomes a SubBuildArgs
    task at all (fragments_h5.py's skip-zero-mapped-contigs check).
    """
    captured_total = {}
    orig_tqdm = fh5_module.tqdm

    class _CapturingTqdm(orig_tqdm):
        def __init__(self, *args, **kwargs):
            captured_total["total"] = kwargs.get("total")
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(fh5_module, "tqdm", _CapturingTqdm)
    monkeypatch.setattr(fh5_module, "GENOMIC_CHUNK_SIZE", 700)

    captured_args = []
    original_worker = fh5_module.build_sub_fragments_h5

    def _spy(args):
        captured_args.append(args)
        return original_worker(args)

    monkeypatch.setattr(fh5_module, "build_sub_fragments_h5", _spy)

    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "in.bam")
        len_a, len_z = 2000, 500
        reads = [
            {"name": f"r{i}", "contig": "chrA", "pos": i * 100, "cigar": "10M"}
            for i in range(3)
        ]
        # chrZ has zero mapped reads and must be entirely excluded from the total.
        _write_se_bam(bam, [("chrA", len_a), ("chrZ", len_z)], reads)

        h5 = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5,
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
        )

    assert len(captured_args) >= 2, "expected chrA to split into multiple chunks"
    assert any(a.chunk_start != 0 for a in captured_args), (
        "expected at least one chunk with chunk_start != 0 -- otherwise "
        "chunk_stop - chunk_start and contig_length - chunk_start are "
        "numerically identical and the test cannot distinguish them"
    )
    assert all(a.bam_contig == "chrA" for a in captured_args), (
        "chrZ (zero mapped fragments) must not produce any worker task"
    )

    expected_total = sum(a.chunk_stop - a.chunk_start for a in captured_args)
    assert expected_total == len_a
    assert captured_total["total"] == expected_total


# ── Item 4: the design's central claim -- SubBuildArgs is not subscriptable ──

def test_sub_build_args_not_subscriptable():
    """The whole point of replacing the 17-element positional tuple with a
    dataclass (docs/architecture/worker_args_refactor.md §3.1) is that `a[4] -
    a[3]`-style access must be impossible: a NamedTuple is a tuple subclass
    and would keep positional access working, silently preserving the exact
    accessor pattern that shipped the v2.10.1 TypeError. Pin that claim
    directly.
    """
    args = SubBuildArgs(
        input_fname="in.bam",
        bam_contig="chr1",
        output_contig="chr1",
        chunk_start=0,
        chunk_stop=100,
        contig_length=100,
        fasta_filename=None,
        single_end=True,
        se_max_fragment_length=None,
        read_gc=False,
        read_strand=True,
        read_methyl=False,
        set_mapq_255_to_none=False,
        include_duplicates=False,
        store_fragment_end_clipped=True,
        tmp_dir_name="/tmp",
        min_mapq=None,
    )
    with pytest.raises(TypeError):
        args[0]
