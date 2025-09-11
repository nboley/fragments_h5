import os
import pytest
import subprocess
import tempfile

from fragments_h5.fragments_h5 import build_fragments_h5, FragmentsH5, bam_to_fragments

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "./data/")

# importing from datasets was giving me a circular import that i couldn't resolve, so copied this list here
GREATEST_HITS = [
    ("chr6", 99119615, 99119634),  # CTCF Constitutive
]


@pytest.fixture(scope="module")
def bam_path():
    return os.path.join(DATA_DIR, "./small.chr6.bam")


@pytest.fixture(scope="module")
def target_bam_path():
    return os.path.join(DATA_DIR, "./scATAC_breast_v1_chr6_99118615_99121634.hg38.bam")


@pytest.fixture(scope="module")
def fasta_file_path():
    return os.path.join(DATA_DIR, "./GRCh38.p12.genome.chr6_99110000_99130000.fa.gz")


@pytest.fixture(scope="module")
def small_h5_path(bam_path, fasta_file_path):
    with tempfile.TemporaryDirectory() as dirname:
        ofname = os.path.join(dirname, os.path.basename(bam_path) + ".frag.h5")
        build_fragments_h5(
            bam_path, ofname, "hg38", fasta_file=fasta_file_path
        )
        yield ofname


@pytest.fixture(scope="module")
def target_h5_path(target_bam_path, fasta_file_path):
    with tempfile.TemporaryDirectory() as dirname:
        ofname = os.path.join(
            dirname, "scATAC_breast_v1.chr6_99118615_99121634.fragments.h5"
        )

        cmd = f"""
        build-fragments-h5 \
            {target_bam_path} \
            {ofname} \
            --reference hg38 \
            --contigs chr6 \
            --fasta {fasta_file_path} \
            --verbose \
        """.strip()
        subprocess.run(cmd, shell=True, check=True)
        yield ofname


def test_build_small_h5(small_h5_path):
    # just test that the fixture builds the fragment h5
    # this also helps to make the timings more interpretable
    pass


def test_build_target_h5(target_h5_path):
    # just test that the fixture builds the fragment h5
    # this also helps to make the timings more interpretable
    pass


def fragments_eq(f1, f2):
    # we can't use the standard fragment __eq__ b/c the gc
    # storage is lossy, and so we can only check that they're close
    return (
        f1.chrom == f2.chrom
        and f1.start == f2.start
        and f1.stop == f2.stop
        and f1.mapq1 == f2.mapq1
        and f1.mapq2 == f2.mapq2
        and (
            (f1.gc is None and f2.gc is None)
            or (f1.gc - f2.gc < 1e-2)
        )
    )


def assert_fragments_identical(fs1, fs2):
    assert all(fragments_eq(f1, f2) for f1, f2 in zip(fs1, fs2))


def test_read_all(small_h5_path, bam_path, fasta_file_path):
    fmh5 = FragmentsH5(small_h5_path, "r")
    h5_fragments = list(fmh5.fetch(return_gc=True))
    bam_fragments = list(
        bam_to_fragments(
            bam_path, max_tlen=fmh5.max_fragment_length, fasta_file=fasta_file_path
        )
    )
    assert_fragments_identical(h5_fragments, bam_fragments)


def test_fetch(small_h5_path, bam_path, fasta_file_path):
    fmh5 = FragmentsH5(small_h5_path, "r")
    for contig, start, stop in GREATEST_HITS:
        h5_fragments = list(fmh5.fetch(contig, start, stop, return_gc=True))
        bam_fragments = list(
            bam_to_fragments(
                bam_path,
                chrom=contig,
                start=start,
                stop=stop,
                max_tlen=fmh5.max_fragment_length,
                fasta_file=fasta_file_path,
            )
        )
        assert_fragments_identical(h5_fragments, bam_fragments)


def test_fetches_over_target_region(target_h5_path, target_bam_path, fasta_file_path):
    fmh5 = FragmentsH5(target_h5_path, "r")

    contig, region_start, region_stop = "chr6", 99119615, 99119634
    for offset in range(0, 501, 100):
        h5_fragments = list(
            fmh5.fetch(contig, region_start - offset, region_stop + offset, return_gc=True)
        )
        bam_fragments = []
        for fragment in bam_to_fragments(
            target_bam_path,
            chrom=contig,
            start=region_start - offset,
            stop=region_stop + offset,
            max_tlen=fmh5.max_fragment_length,
            fasta_file=fasta_file_path,
        ):
            bam_fragments.append(fragment)
        assert_fragments_identical(h5_fragments, bam_fragments)
        print(offset, len(h5_fragments), len(bam_fragments))


def test_fetch_counts(target_h5_path):
    frag_h5 = FragmentsH5(target_h5_path)
    counts = frag_h5.fetch_counts("chr6", 99118615, 99121634)
    assert counts == 184


def test_fetch_cache_v_no_cache(small_h5_path):
    assert list(
        FragmentsH5(small_h5_path, "r", cache_pointers=True).fetch(*GREATEST_HITS[0])
    ) == list(
        FragmentsH5(small_h5_path, "r", cache_pointers=False).fetch(*GREATEST_HITS[0])
    )
