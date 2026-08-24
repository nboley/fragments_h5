"""Tests for fragment selection correctness changes (items 1-3).

Tests secondary alignment exclusion, SE oversized span raising, and
contig-with-single-mapped-read inclusion.
"""
import logging
import os
import tempfile

import numpy
import pysam
import pytest

from fragments_h5.fragments_h5 import build_fragments_h5, FragmentsH5


def _write_se_bam(path, contigs, reads):
    """Write a single-end BAM with the given contigs and reads.

    contigs: list of (name, length)
    reads: list of dicts with keys: name, contig, pos, cigar, mapq, flag_bits (optional)
    """
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
            a.cigarstring = r["cigar"]
            a.mapping_quality = r.get("mapq", 60)
            a.query_sequence = "A" * r.get("seq_len", 10)
            a.query_qualities = pysam.qualitystring_to_array("I" * r.get("seq_len", 10))
            a.flag = r.get("flag", 0)
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            outf.write(a)
    pysam.sort("-o", path, path)
    pysam.index(path)


def _write_pe_bam(path, contigs, pairs):
    """Write a paired-end BAM with proper pairs.

    pairs: list of dicts with keys: name, contig, pos1, pos2, tlen, cigar, mapq,
           extra_flags (optional, added to read1/read2)
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as outf:
        for p in pairs:
            seq_len = p.get("seq_len", 10)
            extra = p.get("extra_flags", 0)
            tlen = p["tlen"]
            # Read 1 (forward)
            a1 = pysam.AlignedSegment()
            a1.query_name = p["name"]
            a1.reference_id = outf.get_tid(p["contig"])
            a1.reference_start = p["pos1"]
            a1.cigarstring = p.get("cigar", f"{seq_len}M")
            a1.mapping_quality = p.get("mapq", 60)
            a1.query_sequence = "A" * seq_len
            a1.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
            a1.flag = 0x1 | 0x2 | 0x20 | 0x40 | extra  # paired, proper, mate_reverse, read1
            a1.next_reference_id = outf.get_tid(p["contig"])
            a1.next_reference_start = p["pos2"]
            a1.template_length = tlen
            outf.write(a1)
            # Read 2 (reverse)
            a2 = pysam.AlignedSegment()
            a2.query_name = p["name"]
            a2.reference_id = outf.get_tid(p["contig"])
            a2.reference_start = p["pos2"]
            a2.cigarstring = p.get("cigar", f"{seq_len}M")
            a2.mapping_quality = p.get("mapq", 60)
            a2.query_sequence = "A" * seq_len
            a2.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
            a2.flag = 0x1 | 0x2 | 0x10 | 0x80 | extra  # paired, proper, reverse, read2
            a2.next_reference_id = outf.get_tid(p["contig"])
            a2.next_reference_start = p["pos1"]
            a2.template_length = -tlen
            outf.write(a2)
    pysam.sort("-o", path, path)
    pysam.index(path)


def _build_se_h5(bam_path, h5_path, se_max_fragment_length=1000, **kwargs):
    build_fragments_h5(
        bam_path, h5_path,
        single_end=True,
        se_max_fragment_length=se_max_fragment_length,
        num_processes=1,
        read_strand=True,
        store_fragment_end_clipped=False,
        **kwargs,
    )


def _build_pe_h5(bam_path, h5_path, **kwargs):
    build_fragments_h5(
        bam_path, h5_path,
        num_processes=1,
        read_strand=True,
        store_fragment_end_clipped=False,
        **kwargs,
    )


# --- Item 1: Secondary alignment exclusion ---

def test_secondary_alignments_excluded_single_end():
    """5 primary reads + 3 secondary -> should get 5 fragments, not 8."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "se_secondary.bam")
        # 5 primary reads
        reads = [
            {"name": f"read_{i}", "contig": "chr1", "pos": i * 100, "cigar": "10M", "flag": 0}
            for i in range(5)
        ]
        # 3 secondary alignments at positions 5000, 5100, 5200
        for i, pos in enumerate([5000, 5100, 5200]):
            reads.append({
                "name": f"sec_{i}", "contig": "chr1", "pos": pos,
                "cigar": "10M", "flag": 0x100,  # is_secondary
            })
        _write_se_bam(bam, [("chr1", 10000)], reads)

        h5 = os.path.join(tmpdir, "out.h5")
        _build_se_h5(bam, h5)

        with FragmentsH5(h5) as fh5:
            starts, stops, _ = fh5.fetch_array("chr1")
            assert len(starts) == 5, f"Expected 5 fragments, got {len(starts)}"
            # None of the secondary positions should appear
            for sec_pos in [5000, 5100, 5200]:
                assert sec_pos not in starts, f"Secondary alignment at {sec_pos} should be excluded"


def test_secondary_alignments_excluded_paired_end():
    """6 proper pairs + 1 secondary with positive TLEN -> should get 6, not 7."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "pe_secondary.bam")

        pairs = [
            {"name": f"pair_{i}", "contig": "chr1", "pos1": i * 200, "pos2": i * 200 + 100,
             "tlen": 110, "cigar": "10M", "mapq": 60}
            for i in range(6)
        ]

        # Write pairs + one secondary read manually
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 100000}],
        }
        with pysam.AlignmentFile(bam, "wb", header=header) as outf:
            for p in pairs:
                seq_len = 10
                # Read 1
                a1 = pysam.AlignedSegment()
                a1.query_name = p["name"]
                a1.reference_id = 0
                a1.reference_start = p["pos1"]
                a1.cigarstring = "10M"
                a1.mapping_quality = 60
                a1.query_sequence = "A" * seq_len
                a1.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
                a1.flag = 0x1 | 0x2 | 0x20 | 0x40  # paired, proper, mate_reverse, read1
                a1.next_reference_id = 0
                a1.next_reference_start = p["pos2"]
                a1.template_length = p["tlen"]
                outf.write(a1)
                # Read 2
                a2 = pysam.AlignedSegment()
                a2.query_name = p["name"]
                a2.reference_id = 0
                a2.reference_start = p["pos2"]
                a2.cigarstring = "10M"
                a2.mapping_quality = 60
                a2.query_sequence = "A" * seq_len
                a2.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
                a2.flag = 0x1 | 0x2 | 0x10 | 0x80  # paired, proper, reverse, read2
                a2.next_reference_id = 0
                a2.next_reference_start = p["pos1"]
                a2.template_length = -p["tlen"]
                outf.write(a2)

            # Secondary alignment with positive in-range TLEN at pos 70000
            sec = pysam.AlignedSegment()
            sec.query_name = "sec_read"
            sec.reference_id = 0
            sec.reference_start = 70000
            sec.cigarstring = "10M"
            sec.mapping_quality = 60
            sec.query_sequence = "A" * 10
            sec.query_qualities = pysam.qualitystring_to_array("I" * 10)
            sec.flag = 0x1 | 0x2 | 0x20 | 0x40 | 0x100  # paired, proper, mate_reverse, read1, SECONDARY
            sec.next_reference_id = 0
            sec.next_reference_start = 70100
            sec.template_length = 110
            outf.write(sec)

        pysam.sort("-o", bam, bam)
        pysam.index(bam)

        h5 = os.path.join(tmpdir, "out.h5")
        _build_pe_h5(bam, h5)

        with FragmentsH5(h5) as fh5:
            starts, stops, _ = fh5.fetch_array("chr1")
            assert len(starts) == 6, f"Expected 6 fragments, got {len(starts)}"
            assert 70000 not in starts, "Secondary alignment at 70000 should be excluded"


# --- Item 2: SE oversized span ---

def test_se_oversized_span_raises():
    """A read with reference span 70020 should raise ValueError with useful info."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "oversized.bam")
        reads = [
            {"name": "normal", "contig": "chr1", "pos": 100, "cigar": "10M", "seq_len": 10, "flag": 0},
            {"name": "oversized_read", "contig": "chr1", "pos": 1000, "cigar": "10M70000N10M", "seq_len": 20, "flag": 0},
        ]
        _write_se_bam(bam, [("chr1", 100000)], reads)

        h5 = os.path.join(tmpdir, "out.h5")
        with pytest.raises(ValueError, match="70020") as exc_info:
            _build_se_h5(bam, h5, se_max_fragment_length=None)
        msg = str(exc_info.value)
        assert "chr1" in msg
        assert "oversized_read" in msg
        assert "1000" in msg


def test_se_span_boundary():
    """Span of 65535 should succeed; 65536 should raise ValueError."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Span = 65535: 10M + (65535-20)N + 10M
        n_ok = 65535 - 20
        bam_ok = os.path.join(tmpdir, "boundary_ok.bam")
        reads_ok = [
            {"name": "read_ok", "contig": "chr1", "pos": 100,
             "cigar": f"10M{n_ok}N10M", "seq_len": 20, "flag": 0},
        ]
        _write_se_bam(bam_ok, [("chr1", 200000)], reads_ok)

        h5_ok = os.path.join(tmpdir, "ok.h5")
        _build_se_h5(bam_ok, h5_ok, se_max_fragment_length=None)
        with FragmentsH5(h5_ok) as fh5:
            starts, stops, _ = fh5.fetch_array("chr1")
            assert len(starts) == 1
            lengths = stops - starts
            assert lengths[0] == 65535

        # Span = 65536: 10M + (65536-20)N + 10M
        n_bad = 65536 - 20
        bam_bad = os.path.join(tmpdir, "boundary_bad.bam")
        reads_bad = [
            {"name": "read_bad", "contig": "chr1", "pos": 100,
             "cigar": f"10M{n_bad}N10M", "seq_len": 20, "flag": 0},
        ]
        _write_se_bam(bam_bad, [("chr1", 200000)], reads_bad)

        h5_bad = os.path.join(tmpdir, "bad.h5")
        with pytest.raises(ValueError, match="65536"):
            _build_se_h5(bam_bad, h5_bad, se_max_fragment_length=None)


def test_se_max_fragment_length_still_skips():
    """With se_max_fragment_length=1000, oversized reads are skipped, not raised."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "skip.bam")
        reads = [
            {"name": "normal", "contig": "chr1", "pos": 100, "cigar": "10M", "seq_len": 10, "flag": 0},
            {"name": "oversized", "contig": "chr1", "pos": 1000, "cigar": "10M70000N10M", "seq_len": 20, "flag": 0},
        ]
        _write_se_bam(bam, [("chr1", 100000)], reads)

        h5 = os.path.join(tmpdir, "out.h5")
        _build_se_h5(bam, h5, se_max_fragment_length=1000)

        with FragmentsH5(h5) as fh5:
            starts, stops, _ = fh5.fetch_array("chr1")
            assert len(starts) == 1, "Only the 10bp fragment should survive"
            lengths = stops - starts
            assert lengths[0] == 10


# --- Item 3: Contig with single mapped read ---

def test_contig_with_single_mapped_read_included():
    """A contig with exactly 1 mapped read should be included, not skipped."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "single_read.bam")
        reads = [
            # contigA: 5 reads
            *[{"name": f"a_{i}", "contig": "contigA", "pos": i * 100, "cigar": "10M", "flag": 0} for i in range(5)],
            # contigB: 1 read
            {"name": "b_0", "contig": "contigB", "pos": 5000, "cigar": "10M", "flag": 0},
        ]
        _write_se_bam(bam, [("contigA", 10000), ("contigB", 10000)], reads)

        h5 = os.path.join(tmpdir, "out.h5")
        _build_se_h5(bam, h5)

        with FragmentsH5(h5) as fh5:
            assert "contigA" in fh5.data
            assert "contigB" in fh5.data
            starts_a, _, _ = fh5.fetch_array("contigA")
            starts_b, _, _ = fh5.fetch_array("contigB")
            assert len(starts_a) == 5
            assert len(starts_b) == 1
            assert starts_b[0] == 5000


def test_contig_with_zero_mapped_reads_skipped(caplog, monkeypatch):
    """A contig with 0 mapped reads should be skipped without dispatching a
    worker task for it, and an info message logged.

    The skip predicate (fragments_h5.py:1083-1085) doesn't change what ends
    up in the output -- build_sub_fragments_h5 already returns None for a
    zero-fragment chunk, so `"contigB" not in fh5.data` holds either way.
    What the skip predicate actually changes is whether a worker task is
    dispatched for contigB at all. So the `"contigB" not in fh5.data`
    assertion below is a baseline sanity check (data absence, expected
    regardless of the skip logic), not evidence the skip ran. The two
    assertions with teeth are: (1) the log message, and (2) the monkeypatch
    spy confirming build_sub_fragments_h5 is never invoked for contigB --
    that's the direct, observable difference the skip logic produces.
    """
    import fragments_h5.fragments_h5 as fh5_module

    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "zero_reads.bam")
        reads = [
            *[{"name": f"a_{i}", "contig": "contigA", "pos": i * 100, "cigar": "10M", "flag": 0} for i in range(5)],
        ]
        # Include contigB in header but no reads for it
        _write_se_bam(bam, [("contigA", 10000), ("contigB", 10000)], reads)

        h5 = os.path.join(tmpdir, "out.h5")

        called_contigs = []
        original_worker = fh5_module.build_sub_fragments_h5

        def _spy(args):
            called_contigs.append(args[1])  # bam_contig
            return original_worker(args)

        monkeypatch.setattr(fh5_module, "build_sub_fragments_h5", _spy)

        with caplog.at_level(logging.INFO):
            _build_se_h5(bam, h5)

        assert "contigB" not in called_contigs, "No worker task should be dispatched for contigB"
        assert "contigA" in called_contigs

        assert "contigB" not in FragmentsH5(h5).data
        assert any("skipping" in r.message and "zero mapped" in r.message for r in caplog.records)


def test_pe_contig_with_single_mapped_read_skipped(caplog, monkeypatch):
    """A PE contig with exactly 1 mapped read (mate elsewhere/unmapped) should be
    skipped without dispatching a worker task for it, and an info message logged.

    A lone mapped alignment on a contig cannot form a fragment in paired-end mode:
    if its mate were on the same contig, get_index_statistics().mapped would be
    >= 2. As in test_contig_with_zero_mapped_reads_skipped, asserting
    "contigB" not in fh5.data would be vacuous here -- a single unpaired read
    produces no fragment regardless of whether the skip logic runs. The
    assertions with teeth are the worker-dispatch spy and the log message.
    """
    import fragments_h5.fragments_h5 as fh5_module

    with tempfile.TemporaryDirectory() as tmpdir:
        bam = os.path.join(tmpdir, "pe_single_read.bam")

        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "contigA", "LN": 100000}, {"SN": "contigB", "LN": 10000}],
        }
        with pysam.AlignmentFile(bam, "wb", header=header) as outf:
            # contigA: 5 normal proper pairs
            for i in range(5):
                seq_len = 10
                pos1 = i * 200
                pos2 = pos1 + 100
                tlen = 110
                a1 = pysam.AlignedSegment()
                a1.query_name = f"pair_{i}"
                a1.reference_id = outf.get_tid("contigA")
                a1.reference_start = pos1
                a1.cigarstring = "10M"
                a1.mapping_quality = 60
                a1.query_sequence = "A" * seq_len
                a1.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
                a1.flag = 0x1 | 0x2 | 0x20 | 0x40  # paired, proper, mate_reverse, read1
                a1.next_reference_id = outf.get_tid("contigA")
                a1.next_reference_start = pos2
                a1.template_length = tlen
                outf.write(a1)
                a2 = pysam.AlignedSegment()
                a2.query_name = f"pair_{i}"
                a2.reference_id = outf.get_tid("contigA")
                a2.reference_start = pos2
                a2.cigarstring = "10M"
                a2.mapping_quality = 60
                a2.query_sequence = "A" * seq_len
                a2.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
                a2.flag = 0x1 | 0x2 | 0x10 | 0x80  # paired, proper, reverse, read2
                a2.next_reference_id = outf.get_tid("contigA")
                a2.next_reference_start = pos1
                a2.template_length = -tlen
                outf.write(a2)

            # contigB: a single mapped read whose mate is unmapped -- exactly one
            # mapped alignment on the contig, no possible partner for a fragment.
            lone = pysam.AlignedSegment()
            lone.query_name = "lone_read"
            lone.reference_id = outf.get_tid("contigB")
            lone.reference_start = 5000
            lone.cigarstring = "10M"
            lone.mapping_quality = 60
            lone.query_sequence = "A" * 10
            lone.query_qualities = pysam.qualitystring_to_array("I" * 10)
            lone.flag = 0x1 | 0x8 | 0x40  # paired, mate_unmapped, read1
            lone.next_reference_id = -1
            lone.next_reference_start = -1
            lone.template_length = 0
            outf.write(lone)

        pysam.sort("-o", bam, bam)
        pysam.index(bam)

        # Anchor the test to the actual precondition rather than an assumption
        # about how pysam counted it.
        with pysam.AlignmentFile(bam) as bam_fp:
            stats = {s.contig: s.mapped for s in bam_fp.get_index_statistics()}
        assert stats["contigB"] == 1

        h5 = os.path.join(tmpdir, "out.h5")

        called_contigs = []
        original_worker = fh5_module.build_sub_fragments_h5

        def _spy(args):
            called_contigs.append(args[1])  # bam_contig
            return original_worker(args)

        monkeypatch.setattr(fh5_module, "build_sub_fragments_h5", _spy)

        with caplog.at_level(logging.INFO):
            _build_pe_h5(bam, h5)

        assert "contigB" not in called_contigs, "No worker task should be dispatched for contigB"
        assert "contigA" in called_contigs

        assert "contigB" not in FragmentsH5(h5).data
        assert any("skipping" in r.message and "zero mapped" in r.message for r in caplog.records)
