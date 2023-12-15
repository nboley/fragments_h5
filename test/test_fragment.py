from copy import copy
import subprocess
import numpy
import pysam
from fragments_h5.fragment import (
    Fragment,
    get_methyl_idxs_for_pair_of_alignments,
    bam_to_align,
    alignment_to_fragment_start_stop,
)


def bash_safe(cmd, print_cmd=True):
    """
    Runs a bash command as a subprocess, and check that it completes successfully.  Stdout/Stderr are forwarded
    to the user's terminal as a normal subprocess.run is.  However, this forwarding behavior is broken in jupyter
    because they change the default sys.stdout/sys.stderr, so beware.

    This function is called bash_safe because it prepends the command with set -eo pipefail.
    This causes the process to exit immediately if any commands return non-zero (set -e), or any commands in a
    pipeline fail (set -o pipefail)

    :param cmd: the bash command to run
    :param print_cmd: print the command to stdout
    :return:
    """
    if print_cmd:
        print("+ ", cmd)
    subprocess.run(
        "set -eo pipefail && " + cmd, shell=True, executable="/bin/bash", check=True
    )


def test_fragment_eq():
    f1 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    f2 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    assert f1 == f2

    f1 = Fragment(chrom="chr1", start=629257, stop=629421)
    f2 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    assert f1 != f2

    f1 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    f2 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13)
    assert f1 == f2

    f1 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    f2 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=0.5)
    assert f1 != f2


def _make_aligns_from_frag(
    alignment_file,
    chrom,
    start,
    stop,
    read_length=75,
    is_duplicate=False,
    is_supplementary=False,
    n_soft_clips=0,
):
    tid = alignment_file.get_tid(chrom)
    if read_length > stop - start:
        raise ValueError(
            f"read length is larger than fragment length.  Although this exists in real data, this "
            f"scenario is unsupported by this function"
        )

    def make_align(reference_start, is_read1, tlen):
        align = pysam.AlignedSegment()
        align.tid = tid
        align.is_read1 = is_read1
        # align.mate_is_reverse = is_read1
        # align.is_reverse = not is_read1
        align.reference_start = reference_start
        align.query_name = "read1"
        align.mapping_quality = 30
        align.query_sequence = "A" * read_length
        align.query_qualities = [30] * read_length
        cigarstring = f"{read_length - n_soft_clips}M"
        if n_soft_clips > 0:
            cigarstring = f"{n_soft_clips}S{cigarstring}"
        align.cigarstring = cigarstring
        align.tlen = tlen
        align.is_duplicate = is_duplicate
        align.is_supplementary = is_supplementary
        # MQ = mate quality
        align.set_tag("MQ", 10)
        align.set_tag("MC", "30M")
        return align

    a1 = make_align(start, is_read1=True, tlen=stop - start)
    a2 = make_align(stop - read_length, is_read1=False, tlen=-(stop - start))
    return a1, a2


def test_methyl_bam_to_frag():
    return True
    bam_file = "s3://rbio-pipe/v5/SEQRUN48/output/fastq_to_bam/BO_14_enzymatic_methyl/align.bam"
    chrom, start, stop = "chr4", 101791300, 101791400

    fasta_file = pysam.FastaFile("./GRCh38.p12.genome.fa.gz")

    align_iter = bam_to_align(
        bam_file,
        chrom=chrom,
        start=start,
        stop=stop,
        fasta_file=fasta_file,
        include_neg_tlen=True,
    )

    alignments = [align for align in align_iter]

    for alt_pair_idxs, tlen, expected_gc_pos, expected_meth in zip(
        [[115, 116], [111, 139], [226, 264]],
        [
            135,  # test case 1: perfect overlap between reads
            285,  # test case 2: some overlap between reads
            315,  # test case 3: no overlap between reads
        ],
        [
            [3, 48, 74, 131],
            [
                11,
                56,
                82,
                139,
                151,
                167,
                180,
                188,
                206,
                215,
                217,
                228,
                243,
                266,
                268,
                283,
            ],
            [2, 57, 62, 73, 75, 94, 203, 205, 219, 238, 256, 269, 274, 288],
        ],
        [[False] * 4, [False] * 16, [False] * 4 + [True] + [False] * 9],
    ):

        align_1 = alignments[alt_pair_idxs[0]]
        align_2 = alignments[alt_pair_idxs[1]]
        phred_1 = copy(align_1.qual)
        phred_2 = copy(align_2.qual)
        assert align_1.tlen == tlen

        gc_pos, is_meth = get_methyl_idxs_for_pair_of_alignments(
            align_1, align_2, fasta_file
        )
        assert (gc_pos == numpy.array(expected_gc_pos)).all()
        assert (is_meth == numpy.array(expected_meth)).all()

        frag_start, frag_stop = alignment_to_fragment_start_stop(align_1)
        frag_seq = fasta_file.fetch(align_1.reference_name, frag_start, frag_stop)

        # Check that the sequence has GCs at the expected positions
        for pos in expected_gc_pos:
            assert frag_seq[pos : pos + 2] == "CG"

        # Check that align_1 has GC at the expected position
        for gg, _gc_pos in enumerate(gc_pos):
            # Falls in the first fragment
            if _gc_pos < len(align_1.seq):
                assert (align_1.seq[_gc_pos : _gc_pos + 2] == "CG") == is_meth[gg]
            # Falls in the second fragment (not exclusive of falling in first fragment)
            if _gc_pos > len(frag_seq) - len(align_2.seq):
                pos_offset = len(frag_seq) - len(align_2.seq)
                assert (
                    align_2.seq[_gc_pos - pos_offset : _gc_pos - pos_offset + 2] == "CG"
                ) == is_meth[gg]

        # Change the unmethylated "CA" to methylated "CG" or vice-versa
        for gg, _gc_pos in enumerate(gc_pos):
            if _gc_pos < len(align_1.seq):
                # If it's methylated "CG", change to "CA"
                # We want to overwrite align_2 if we change align_1 for testing
                orig_seq = copy(align_1.seq)
                # From pysam docs:
                # One issue to look out for is that the sequence should always be set before the quality scores.
                # Setting the sequence will also erase any quality scores that were set previously
                if is_meth[gg]:
                    align_1.seq = (
                        align_1.seq[:_gc_pos] + "CA" + align_1.seq[_gc_pos + 2 :]
                    )
                else:
                    assert not is_meth[gg]
                    align_1.seq = (
                        align_1.seq[:_gc_pos] + "CG" + align_1.seq[_gc_pos + 2 :]
                    )
                align_1.qual = phred_1[: _gc_pos + 1] + chr(76) + phred_1[_gc_pos + 2 :]

                modified_gc_pos, modified_is_meth = get_methyl_idxs_for_pair_of_alignments(
                    align_1, align_2, fasta_file
                )
                assert (
                    modified_is_meth[gg] != is_meth[gg]
                ), f"{orig_seq} modified to {align_1.seq} at position {_gc_pos}, with reference sequence {frag_seq}"

            if _gc_pos > len(frag_seq) - len(align_2.seq):
                pos_offset = len(frag_seq) - len(align_2.seq)
                orig_seq = copy(align_2.seq)
                if is_meth[gg]:
                    align_2.seq = (
                        align_2.seq[: _gc_pos - pos_offset]
                        + "CA"
                        + align_2.seq[_gc_pos - pos_offset + 2 :]
                    )
                else:
                    assert not is_meth[gg]
                    align_2.seq = (
                        align_2.seq[: _gc_pos - pos_offset]
                        + "CG"
                        + align_2.seq[_gc_pos - pos_offset + 2 :]
                    )
                align_2.qual = (
                    phred_2[: _gc_pos - pos_offset + 1]
                    + chr(76)
                    + phred_2[_gc_pos - pos_offset + 2 :]
                )
                modified_gc_pos, modified_is_meth = get_methyl_idxs_for_pair_of_alignments(
                    align_1, align_2, fasta_file
                )
                assert (
                    modified_is_meth[gg] != is_meth[gg]
                ), f"{orig_seq} modified to {align_2.seq} at position {_gc_pos}, with reference sequence {frag_seq}"

            # Test the case where PHREDs are identical, and methylations disagree between the reads
            if len(align_1.seq) == len(frag_seq) == len(align_2.seq):
                align_1.seq = align_1.seq[:_gc_pos] + "CA" + align_1.seq[_gc_pos + 2 :]
                align_2.seq = align_2.seq[:_gc_pos] + "CG" + align_2.seq[_gc_pos + 2 :]
                align_2.qual = phred_2
                align_1.qual = align_2.qual
                pos_is_meth = []
                numpy.random.seed(0)
                for ii in range(100):
                    disagree_gc_pos, disagree_is_meth = get_methyl_idxs_for_pair_of_alignments(
                        align_1, align_2, fasta_file, methylation_prior=0.2
                    )
                    pos_is_meth.append(disagree_is_meth[gg])
                assert 10 <= numpy.array(pos_is_meth).sum() <= 30, (
                    gg,
                    numpy.array(pos_is_meth).sum(),
                )
