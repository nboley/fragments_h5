"""Tests for --contig-name-map / contig_name_map (fragments_h5.py:948, 1062-1065,
main.py:145-159).

Coverage before this file: ZERO (confirmed by grep). This is the flag that makes
output_contig != bam_contig, which is exactly the field whose insertion into the
17-element worker-args tuple shipped a released TypeError (v2.10.1, see
docs/pending/worker_args_refactor.md and CLAUDE memory). Every assertion here is
written to distinguish a MAPPED output name/value from an UNMAPPED one -- a test
that would pass whether or not the rename happened is worthless for that purpose.
"""
import os
import subprocess
import sys
import tempfile

import numpy
import pysam
import pytest

from fragments_h5.fragments_h5 import build_fragments_h5, FragmentsH5

PYTHON = sys.executable


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


def _write_fasta(path, contigs):
    """contigs: dict of name -> sequence string."""
    with open(path, "w") as f:
        for name, seq in contigs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    pysam.faidx(path)


def _acgt(length):
    return ("ACGT" * (length // 4 + 1))[:length]


def _run_cli(args):
    cmd = [PYTHON, "-m", "fragments_h5.main"] + args
    return subprocess.run(cmd, capture_output=True, text=True)


# ── Central case: full rename ──

def test_contig_name_map_renames_output_group_and_fragments():
    """A mapped BAM contig must appear in the output h5 under the OUTPUT name,
    with fragment positions intact, and the FASTA lookup for GC must also use
    the output name (fragments_h5.py:785 `fasta_chrom = output_contig if
    output_contig != bam_contig else None`).

    The FASTA below contains ONLY the output name. If the worker's FASTA
    lookup fell back to the (unmapped) BAM name, get_g_or_c_cumsum would
    silently return None (fragment.py:424 "chrom not in fasta_file" -> None,
    no exception), and every gc value would be written as the 255 sentinel
    (fragments_h5.py:845), and fetch_array decodes that sentinel back to NaN
    on read (fragments_h5.py:570-574 `gc_none_mask = _tmp_gc == 255`). So an
    all-real (non-NaN) gc array is direct evidence the output name was used
    for the FASTA lookup, not just for the group name.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        contig_len = 2000
        bam = os.path.join(tmpdir, "in.bam")
        reads = [
            {"name": f"r{i}", "contig": "chrA", "pos": i * 100, "cigar": "10M"}
            for i in range(3)
        ]
        _write_se_bam(bam, [("chrA", contig_len)], reads)

        fasta = os.path.join(tmpdir, "ref.fa")
        _write_fasta(fasta, {"chrB": _acgt(contig_len)})

        h5 = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5, fasta_filename=fasta,
            contig_name_map={"chrA": "chrB"},
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
        )

        with FragmentsH5(h5) as fh5:
            assert set(fh5.data.keys()) == {"chrB"}
            assert "chrA" not in fh5.data
            assert fh5.contig_lengths == {"chrB": contig_len}

            starts, stops, supp = fh5.fetch_array("chrB", return_gc=True)
            assert sorted(starts.tolist()) == [0, 100, 200]
            assert sorted(stops.tolist()) == [10, 110, 210]
            assert not numpy.isnan(supp["gc"]).any(), (
                "gc came back NaN -- FASTA lookup used the wrong (unmapped) "
                "contig name (fetch_array decodes the 255 sentinel to NaN)"
            )

            with pytest.raises(KeyError):
                fh5.fetch_array("chrA")


# ── Rename survives multiprocessing (SubBuildArgs pickled across a fork) ──

def test_contig_name_map_renames_output_group_and_fragments_multiprocess():
    """All other renaming tests use num_processes=1, so the rename has never
    been exercised through pool.imap_unordered, where SubBuildArgs (holding
    both bam_contig and output_contig) is pickled to a worker process. This
    is exactly the code path the worker-args-refactor changed. Two contigs
    are used so two real tasks are dispatched across 2 worker processes.

    As in test_contig_name_map_renames_output_group_and_fragments, the FASTA
    contains ONLY the output names, so a non-NaN gc array is direct evidence
    that fasta_chrom = output_contig (fragments_h5.py:824) survived the
    pickle round-trip and was used inside the worker, not just the group
    name written by the main process.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        len_a, len_c = 2000, 3000
        bam = os.path.join(tmpdir, "in.bam")
        reads = [
            {"name": "a0", "contig": "chrA", "pos": 100, "cigar": "10M"},
            {"name": "a1", "contig": "chrA", "pos": 500, "cigar": "10M"},
            {"name": "c0", "contig": "chrC", "pos": 200, "cigar": "10M"},
            {"name": "c1", "contig": "chrC", "pos": 900, "cigar": "10M"},
        ]
        _write_se_bam(bam, [("chrA", len_a), ("chrC", len_c)], reads)

        fasta = os.path.join(tmpdir, "ref.fa")
        # Only the OUTPUT names are present -- a fallback to the BAM name
        # inside the worker would find nothing in the FASTA and silently
        # write the 255 "no gc" sentinel instead of raising.
        _write_fasta(fasta, {"chrB": _acgt(len_a), "chrD": _acgt(len_c)})

        h5 = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5, fasta_filename=fasta,
            contig_name_map={"chrA": "chrB", "chrC": "chrD"},
            single_end=True, se_max_fragment_length=1000,
            num_processes=2, store_fragment_end_clipped=False,
        )

        with FragmentsH5(h5) as fh5:
            assert set(fh5.data.keys()) == {"chrB", "chrD"}
            assert "chrA" not in fh5.data
            assert "chrC" not in fh5.data
            assert fh5.contig_lengths == {"chrB": len_a, "chrD": len_c}

            starts_b, _, supp_b = fh5.fetch_array("chrB", return_gc=True)
            assert sorted(starts_b.tolist()) == [100, 500]
            assert not numpy.isnan(supp_b["gc"]).any(), (
                "gc for renamed contig chrB came back NaN under "
                "num_processes=2 -- the FASTA lookup inside the worker "
                "used the wrong (unmapped/pre-pickle) contig name"
            )

            starts_d, _, supp_d = fh5.fetch_array("chrD", return_gc=True)
            assert sorted(starts_d.tolist()) == [200, 900]
            assert not numpy.isnan(supp_d["gc"]).any(), (
                "gc for renamed contig chrD came back NaN under "
                "num_processes=2 -- the FASTA lookup inside the worker "
                "used the wrong (unmapped/pre-pickle) contig name"
            )

            with pytest.raises(KeyError):
                fh5.fetch_array("chrA")
            with pytest.raises(KeyError):
                fh5.fetch_array("chrC")


# ── Partial map: mapped and unmapped contigs coexist ──

def test_contig_name_map_partial_mapped_and_unmapped_coexist():
    """A map covering only some contigs must rename exactly those, leaving the
    rest under their original BAM name in the SAME output file
    (fragments_h5.py:1065 `_map_name = lambda c: contig_name_map.get(c, c)`).

    The FASTA carries both chrB (chrA's renamed target) and chrC (unchanged).
    Fetching gc for both confirms the per-contig FASTA lookup uses the right
    name in each case -- for chrB that means `fasta_chrom = output_contig if
    output_contig != bam_contig else None` (fragments_h5.py:824) must resolve
    to the renamed name, not silently fall back to "chrA" (which isn't even
    in the FASTA). See test_contig_name_map_renames_output_group_and_fragments
    for why NaN-vs-real gc is the signal that catches a wrong-name lookup."""
    with tempfile.TemporaryDirectory() as tmpdir:
        len_a, len_c = 2000, 3000
        bam = os.path.join(tmpdir, "in.bam")
        reads = [
            {"name": "a0", "contig": "chrA", "pos": 100, "cigar": "10M"},
            {"name": "a1", "contig": "chrA", "pos": 500, "cigar": "10M"},
            {"name": "c0", "contig": "chrC", "pos": 200, "cigar": "10M"},
            {"name": "c1", "contig": "chrC", "pos": 900, "cigar": "10M"},
        ]
        _write_se_bam(bam, [("chrA", len_a), ("chrC", len_c)], reads)

        fasta = os.path.join(tmpdir, "ref.fa")
        # chrA's renamed target ("chrB") plus chrC's own (unchanged) name.
        _write_fasta(fasta, {"chrB": _acgt(len_a), "chrC": _acgt(len_c)})

        h5 = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bam, h5, fasta_filename=fasta,
            contig_name_map={"chrA": "chrB"},  # chrC intentionally absent from the map
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
        )

        with FragmentsH5(h5) as fh5:
            assert set(fh5.data.keys()) == {"chrB", "chrC"}
            assert "chrA" not in fh5.data
            assert fh5.contig_lengths == {"chrB": len_a, "chrC": len_c}

            starts_b, _, supp_b = fh5.fetch_array("chrB", return_gc=True)
            assert sorted(starts_b.tolist()) == [100, 500]
            assert not numpy.isnan(supp_b["gc"]).any(), (
                "gc for renamed contig chrB came back NaN -- FASTA lookup "
                "used the wrong (unmapped) contig name"
            )

            starts_c, _, supp_c = fh5.fetch_array("chrC", return_gc=True)
            assert sorted(starts_c.tolist()) == [200, 900]
            assert not numpy.isnan(supp_c["gc"]).any(), (
                "gc for unmapped contig chrC came back NaN -- FASTA lookup "
                "used the wrong contig name"
            )


# ── --contigs / allowed_contigs interaction ──

def test_allowed_contigs_uses_bam_names_not_output_names():
    """allowed_contigs (--contigs) is resolved against BAM names BEFORE the
    rename is applied (fragments_h5.py:1044-1046 builds contig_lengths from
    bam_fp.references and indexes it with allowed_contigs; contig_name_map is
    applied afterward at :1062-1068). An output name is never a valid entry."""
    with tempfile.TemporaryDirectory() as tmpdir:
        contig_len = 2000
        bam = os.path.join(tmpdir, "in.bam")
        reads = [
            {"name": "a0", "contig": "chrA", "pos": 100, "cigar": "10M"},
            {"name": "d0", "contig": "chrD", "pos": 100, "cigar": "10M"},
        ]
        _write_se_bam(bam, [("chrA", contig_len), ("chrD", contig_len)], reads)

        fasta = os.path.join(tmpdir, "ref.fa")
        _write_fasta(fasta, {"chrE": _acgt(contig_len), "chrD": _acgt(contig_len)})

        # Restricting with the BAM (input) name succeeds and still renames.
        h5_ok = os.path.join(tmpdir, "ok.h5")
        build_fragments_h5(
            bam, h5_ok, fasta_filename=fasta,
            contig_name_map={"chrA": "chrE"},
            allowed_contigs=["chrA"],
            single_end=True, se_max_fragment_length=1000,
            num_processes=1, store_fragment_end_clipped=False,
        )
        with FragmentsH5(h5_ok) as fh5:
            assert set(fh5.data.keys()) == {"chrE"}

        # Restricting with the OUTPUT name fails -- "chrE" is not a BAM contig.
        # This must be caught by the PRIMARY check -- the
        # `{c: all_contig_lengths[c] for c in allowed_contigs}` comprehension
        # (fragments_h5.py:1085) that resolves allowed_contigs against
        # BAM-name-space lengths -- not by an incidental downstream lookup
        # like `num_mapped_fragments[bam_contig]` (fragments_h5.py:1128),
        # which would also raise KeyError('chrE') by coincidence (chrE never
        # appears in bam_fp.get_index_statistics() either) even if the
        # primary check were broken and let "chrE" slip through as a BAM
        # name. Asserting only `pytest.raises(KeyError)` cannot distinguish
        # these two origins, so we inspect the traceback frame where the
        # KeyError was actually raised.
        h5_bad = os.path.join(tmpdir, "bad.h5")
        with pytest.raises(KeyError) as exc_info:
            build_fragments_h5(
                bam, h5_bad, fasta_filename=fasta,
                contig_name_map={"chrA": "chrE"},
                allowed_contigs=["chrE"],
                single_end=True, se_max_fragment_length=1000,
                num_processes=1, store_fragment_end_clipped=False,
            )
        offending_line = str(exc_info.traceback[-1].statement)
        assert "all_contig_lengths[c]" in offending_line, (
            "KeyError was raised from the wrong line -- expected the "
            "allowed_contigs/all_contig_lengths comprehension (the primary "
            f"BAM-name-space check), got: {offending_line!r}"
        )


# ── CLI surface: main.py --contig-name-map file parsing ──

def test_cli_contig_name_map_renames_output():
    """The CLI's TSV map-file parser (main.py:145-159) must produce the same
    rename as passing contig_name_map directly to the library."""
    with tempfile.TemporaryDirectory() as tmpdir:
        contig_len = 2000
        bam = os.path.join(tmpdir, "in.bam")
        _write_se_bam(
            bam, [("chrA", contig_len)],
            [{"name": "r0", "contig": "chrA", "pos": 100, "cigar": "10M"}],
        )

        fasta = os.path.join(tmpdir, "ref.fa")
        _write_fasta(fasta, {"chrB": _acgt(contig_len)})

        map_path = os.path.join(tmpdir, "map.tsv")
        with open(map_path, "w") as f:
            f.write("chrA\tchrB\n")

        h5 = os.path.join(tmpdir, "out.h5")
        result = _run_cli([
            bam, h5, "--fasta", fasta, "--contig-name-map", map_path,
            "--single-end", "--se-max-fragment-length", "1000",
        ])
        assert result.returncode == 0, result.stderr

        with FragmentsH5(h5) as fh5:
            assert set(fh5.data.keys()) == {"chrB"}
            assert "chrA" not in fh5.data


def test_cli_contig_name_map_skips_comments_and_blank_lines():
    """Blank lines and '#'-prefixed comment lines in the map file must be
    skipped (main.py:151-153 `if not line or line.startswith('#'): continue`),
    not treated as malformed entries, and the real mapping line must still
    take effect."""
    with tempfile.TemporaryDirectory() as tmpdir:
        contig_len = 2000
        bam = os.path.join(tmpdir, "in.bam")
        _write_se_bam(
            bam, [("chrA", contig_len)],
            [{"name": "r0", "contig": "chrA", "pos": 100, "cigar": "10M"}],
        )

        fasta = os.path.join(tmpdir, "ref.fa")
        _write_fasta(fasta, {"chrB": _acgt(contig_len)})

        map_path = os.path.join(tmpdir, "map.tsv")
        with open(map_path, "w") as f:
            f.write("# comment line\n\nchrA\tchrB\n\n# trailing comment\n")

        h5 = os.path.join(tmpdir, "out.h5")
        result = _run_cli([
            bam, h5, "--fasta", fasta, "--contig-name-map", map_path,
            "--single-end", "--se-max-fragment-length", "1000",
        ])
        assert result.returncode == 0, result.stderr

        with FragmentsH5(h5) as fh5:
            assert set(fh5.data.keys()) == {"chrB"}
            assert "chrA" not in fh5.data


def test_cli_contig_name_map_malformed_line_rejected():
    """A map line without exactly 2 tab-separated columns must abort the
    build with a clear error and non-zero exit (main.py:154-158), before any
    output file is written."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "in.bam")
        _write_se_bam(
            bam, [("chrA", 2000)],
            [{"name": "r0", "contig": "chrA", "pos": 100, "cigar": "10M"}],
        )

        map_path = os.path.join(tmpdir, "map.tsv")
        with open(map_path, "w") as f:
            f.write("chrA\tchrB\textra\n")

        h5 = os.path.join(tmpdir, "out.h5")
        result = _run_cli([bam, h5, "--contig-name-map", map_path])

        assert result.returncode == 1
        assert "Invalid contig name map line" in result.stderr
        assert not os.path.exists(h5)
