"""Empty genomic chunks must not load reference sequence.

get_g_or_c_cumsum fetches, one-hot encodes and cumsums ~10 Mbp of FASTA. Before
the lazy loader, every chunk of a non-empty contig paid that cost even when it
held no fragment, so build time scaled with genome covered, not fragment count.
"""

import os
import tempfile

import h5py
import numpy as np
import pysam
import pytest

import fragments_h5.fragment as fragment_mod
import fragments_h5.fragments_h5 as fragments_h5_mod
from fragments_h5.fragment import (
    MAX_FRAG_LENGTH,
    bam_to_fragments,
    single_end_bam_to_fragments,
    tsv_to_fragments,
)
from fragments_h5.fragments_h5 import FragmentsH5, build_fragments_h5

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")
BAM = os.path.join(DATA_DIR, "small.chr6.bam")
CHR6_FASTA = os.path.join(DATA_DIR, "GRCh38.p12.genome.chr6_99110000_99130000.fa.gz")

CONTIG = "chrS"
CONTIG_LEN = 10_000
CHUNK = 1_000


@pytest.fixture
def cumsum_calls(monkeypatch):
    calls = []
    real = fragment_mod.get_g_or_c_cumsum

    def counting(*args, **kwargs):
        calls.append((args, kwargs))
        return real(*args, **kwargs)

    monkeypatch.setattr(fragment_mod, "get_g_or_c_cumsum", counting)
    return calls


@pytest.fixture(scope="module")
def synthetic_seq():
    rng = np.random.default_rng(0)
    return "".join(rng.choice(list("ACGT"), size=CONTIG_LEN))


@pytest.fixture(scope="module")
def synthetic_fasta(synthetic_seq):
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "synthetic.fa")
        with open(path, "w") as f:
            f.write(f">{CONTIG}\n")
            for i in range(0, CONTIG_LEN, 80):
                f.write(synthetic_seq[i : i + 80] + "\n")
        pysam.faidx(path)
        yield path


def _write_bed(tmpdir, intervals):
    plain = os.path.join(tmpdir, "frags.bed")
    with open(plain, "w") as f:
        for i, (start, stop) in enumerate(intervals):
            f.write(f"{CONTIG}\t{start}\t{stop}\tfrag{i}\t0\t+\n")
    pysam.tabix_compress(plain, plain + ".gz", force=True)
    pysam.tabix_index(plain + ".gz", preset="bed", force=True)
    return plain + ".gz"


def _expected_gc_byte(seq, start, stop):
    gc = sum(c in "GC" for c in seq[start:stop]) / (stop - start)
    return int(round(round(gc, 5) * 254))


# Chunks of 1 kbp: fragments sit in chunks 1, 4 and 8. Chunks 0, 2, 3, 5, 6, 7, 9 are empty.
OCCUPIED_ONLY = [(1100, 1300), (1500, 1700), (4100, 4300), (8100, 8300)]
# Adds a fragment that starts in chunk 4 and ends in chunk 5, past the chunk boundary.
WITH_BOUNDARY_SPANNER = OCCUPIED_ONLY + [(4950, 5100)]


def test_empty_chunks_do_not_load_fasta(
    synthetic_fasta, monkeypatch, cumsum_calls
):
    monkeypatch.setattr(fragments_h5_mod, "GENOMIC_CHUNK_SIZE", CHUNK)
    with tempfile.TemporaryDirectory() as tmpdir:
        bed = _write_bed(tmpdir, OCCUPIED_ONLY)
        build_fragments_h5(
            bed, os.path.join(tmpdir, "out.h5"),
            fasta_filename=synthetic_fasta, num_processes=1,
        )
    assert len(cumsum_calls) == 3, (
        f"expected one FASTA load per occupied chunk (1, 4, 8), got {len(cumsum_calls)} "
        f"of {CONTIG_LEN // CHUNK} chunks"
    )


def test_output_correct_across_skipped_chunks(
    synthetic_fasta, synthetic_seq, monkeypatch
):
    """Read-back and region queries are correct when most chunks are empty.

    The lazy load cannot change this output: empty chunks were already dropped
    before the merge. This guards a future change that skips chunks before
    dispatch, which could drop fragments or break queries across a gap.
    """
    monkeypatch.setattr(fragments_h5_mod, "GENOMIC_CHUNK_SIZE", CHUNK)
    intervals = sorted(WITH_BOUNDARY_SPANNER)
    with tempfile.TemporaryDirectory() as tmpdir:
        bed = _write_bed(tmpdir, intervals)
        out = os.path.join(tmpdir, "out.h5")
        build_fragments_h5(
            bed, out, fasta_filename=synthetic_fasta, num_processes=1,
        )

        with h5py.File(out, "r") as f:
            grp = f["data"][CONTIG]
            assert len(grp["starts"]) == len(intervals)
            np.testing.assert_array_equal(grp["starts"][:], [s for s, _ in intervals])
            np.testing.assert_array_equal(
                grp["lengths"][:], [e - s for s, e in intervals]
            )
            np.testing.assert_array_equal(
                grp["gc"][:],
                [_expected_gc_byte(synthetic_seq, s, e) for s, e in intervals],
            )

        with FragmentsH5(out) as fh5:
            # Regions chosen to straddle skipped chunks 2-3 and 5-7, and to sit wholly inside them.
            regions = [
                (None, None), (0, CONTIG_LEN), (1600, 4150), (2000, 4000),
                (4200, 8150), (5200, 8000), (9000, 9900), (1250, 1260),
            ]
            for r_start, r_stop in regions:
                starts, stops, supp = fh5.fetch_array(
                    CONTIG, r_start, r_stop, return_gc=True
                )
                lo = 0 if r_start is None else r_start
                hi = CONTIG_LEN if r_stop is None else r_stop
                expected = [(s, e) for s, e in intervals if s < hi and e > lo]
                assert list(zip(starts.tolist(), stops.tolist())) == expected, (r_start, r_stop)
                np.testing.assert_allclose(
                    supp["gc"],
                    [_expected_gc_byte(synthetic_seq, s, e) / 254 for s, e in expected],
                    atol=1e-6,
                )


def test_tsv_empty_region_does_not_load_fasta(synthetic_fasta, cumsum_calls):
    with tempfile.TemporaryDirectory() as tmpdir:
        bed = _write_bed(tmpdir, OCCUPIED_ONLY)
        empty = list(tsv_to_fragments(
            bed, CONTIG, start=2000, stop=3000, fasta_file=synthetic_fasta,
            fasta_region_start=2000, fasta_region_stop=3000 + MAX_FRAG_LENGTH,
            max_tlen=MAX_FRAG_LENGTH,
        ))
        assert empty == []
        assert len(cumsum_calls) == 0

        # Two fragments in one region: the cumsum is built once, not per fragment.
        occupied = list(tsv_to_fragments(
            bed, CONTIG, start=1000, stop=2000, fasta_file=synthetic_fasta,
            fasta_region_start=1000, fasta_region_stop=2000 + MAX_FRAG_LENGTH,
            max_tlen=MAX_FRAG_LENGTH,
        ))
        assert [(f.start, f.stop) for f in occupied] == [(1100, 1300), (1500, 1700)]
        assert all(f.gc is not None for f in occupied)
        assert len(cumsum_calls) == 1


def test_missing_fasta_contig_loads_once(synthetic_fasta, cumsum_calls):
    """get_g_or_c_cumsum returns (None, 0) for a contig absent from the FASTA.

    gc_offset 0 must still count as loaded, or every fragment re-queries the FASTA.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        bed = _write_bed(tmpdir, OCCUPIED_ONLY)
        frags = list(tsv_to_fragments(
            bed, CONTIG, fasta_file=synthetic_fasta, fasta_chrom="not_a_contig",
            max_tlen=MAX_FRAG_LENGTH,
        ))
    assert len(frags) == len(OCCUPIED_ONLY)
    assert all(f.gc is None for f in frags)
    assert len(cumsum_calls) == 1


@pytest.mark.parametrize("to_fragments", [bam_to_fragments, single_end_bam_to_fragments])
def test_bam_empty_region_does_not_load_fasta(to_fragments, cumsum_calls):
    # small.chr6.bam has no reads in chr6:19.0-20.0 Mb and 70 reads in 99.0-100.0 Mb.
    kwargs = dict(
        chrom="chr6", fasta_file=CHR6_FASTA,
        fasta_region_start=19_000_000, fasta_region_stop=20_000_000 + MAX_FRAG_LENGTH,
    )
    if to_fragments is single_end_bam_to_fragments:
        kwargs["se_max_fragment_length"] = MAX_FRAG_LENGTH
    empty = list(to_fragments(BAM, start=19_000_000, stop=20_000_000, **kwargs))
    assert empty == []
    assert len(cumsum_calls) == 0

    kwargs.update(fasta_region_start=99_000_000, fasta_region_stop=100_000_000 + MAX_FRAG_LENGTH)
    occupied = list(to_fragments(BAM, start=99_000_000, stop=100_000_000, **kwargs))
    assert len(occupied) > 0
    assert len(cumsum_calls) == 1
