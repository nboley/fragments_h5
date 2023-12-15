from collections import Counter

import funcsigs
import logging
import numpy
import re
from dataclasses import dataclass

import pysam


from typing import Union, Optional

DEFAULT_MIN_MAPQ = 0


logging.basicConfig(format="%(levelname)s\t%(asctime)-15s\t%(message)s")
log = logging.getLogger(__name__)


def intervals_intersect(x_start, x_stop, y_start, y_stop):
    """Test whether two closed-open intervals intersect.

    >>> intervals_intersect(0,1, 0,1)
    True
    >>> intervals_intersect(0,1, 1,2)
    False
    >>> intervals_intersect(0,1, 0,0)
    False
    >>> intervals_intersect(1,2, 0,1)
    False
    >>> intervals_intersect(1,2, 0,2)
    True
    """
    return x_start < y_stop and y_start < x_stop


def intervals_intersect_with_none(x_start, x_stop, y_start, y_stop):
    """
    >>> intervals_intersect_with_none(0, 1, None, 2)
    True
    >>> intervals_intersect_with_none(0, 1, None, None)
    True
    >>> intervals_intersect_with_none(None, 2, 2, 4)
    False
    >>> intervals_intersect_with_none(None, None, 2, 4)
    True
    >>> intervals_intersect_with_none(2, 4, 4, None)
    False
    """
    if x_start is None:
        x_start = 0
    if y_start is None:
        y_start = 0
    if x_stop is None:
        x_stop = float("inf")
    if y_stop is None:
        y_stop = float("inf")

    return intervals_intersect(x_start, x_stop, y_start, y_stop)


def get_percent_gc(seq):
    """
    >>> get_percent_gc('ACTG')
    0.5
    >>> get_percent_gc('AA')
    0.0
    >>> get_percent_gc('GC')
    1.0
    >>> get_percent_gc('N')
    nan
    >>> get_percent_gc('CH')
    Traceback (most recent call last):
    ...
    AssertionError: Counter({'H': 1}) contains unexpected characters...
    """
    c = Counter(seq)
    try:
        gc = c.pop("G", 0) + c.pop("C", 0)
        at = c.pop("A", 0) + c.pop("T", 0)
        c.pop("N", None)
        assert len(c) == 0, f"{c} contains unexpected characters"
        return gc / (at + gc)
    except ZeroDivisionError:
        return numpy.nan


def unparse_frag_info(x):
    return "" if x is None else str(x)


class ClassPropertyDescriptor(object):
    """
    used by classproperty() to make @classproperty decorator
    """

    def __init__(self, fget, fset=None):
        self.fget = fget
        self.fset = fset

    def __get__(self, obj, klass=None):
        if klass is None:
            klass = type(obj)
        return self.fget.__get__(obj, klass)()

    def __set__(self, obj, value):
        if not self.fset:
            raise AttributeError("can't set attribute")
        type_ = type(obj)
        return self.fset.__get__(obj, type_)(value)

    def setter(self, func):
        if not isinstance(func, (classmethod, staticmethod)):
            func = classmethod(func)
        self.fset = func
        return self


def classproperty(func):
    """
    @classproperty decorator
    """
    if not isinstance(func, (classmethod, staticmethod)):
        func = classmethod(func)

    return ClassPropertyDescriptor(func)


@dataclass
class Fragment:
    __slots__ = [
        "chrom",
        "start",
        "stop",
        "mapq1",
        "mapq2",
        "gc",
        "strand",
        "cell_barcode",
        "num_cpgs",
        "num_meth_cpgs",
    ]

    def __init__(
        self,
        chrom: str,
        start: int,
        stop: int,
        mapq1: Optional[int] = None,
        mapq2: Optional[int] = None,
        gc: Optional[float] = None,
        strand: Optional[str] = None,
        cell_barcode: Optional[str] = None,
        num_cpgs: Optional[int] = None,
        num_meth_cpgs: Optional[int] = None,
    ):
        assert isinstance(start, (int, numpy.int64, numpy.int32))
        if strand in [b"+", b"-"]:
            strand = strand.decode()
        assert strand in [None, "+", "-"]
        for k, v in locals().items():
            # @py_assert2 check is to avoid a pytest error introduced by it changing state
            if k != "self" and not k.startswith("@py_assert"):
                setattr(self, k, v)

    def mapq_gte(self, threshold):
        """
        >>> Fragment('chr1', 0, 1, mapq1=9).mapq_gte(10)
        False
        >>> Fragment('chr1', 0, 1, mapq2=10).mapq_gte(10)
        True

        None mapqs always pass
        >>> Fragment('chr1', 0, 1).mapq_gte(90)
        True
        """
        return self.mapq12_min is None or self.mapq12_min >= threshold

    @property
    def mapq12_min(self):
        """Returns the minimum of self.mapq1 and self.mapq2 ignoring None values
        >>> Fragment('chr1', 0, 1, mapq1=1, mapq2=0).mapq12_min
        0
        >>> Fragment('chr1', 0, 1, mapq1=0, mapq2=0).mapq12_min
        0
        >>> Fragment('chr1', 0, 1, mapq1=None, mapq2=2).mapq12_min
        2
        >>> Fragment('chr1', 0, 1, mapq1=1, mapq2=None).mapq12_min
        1
        >>> Fragment('chr1', 0, 1, mapq1=None, mapq2=None).mapq12_min is None
        True
        """

        def none_to_inf(x):
            return float("inf") if x is None else x

        if self.mapq1 is None and self.mapq2 is None:
            return None
        elif self.mapq1 is None and self.mapq2 is not None:
            return self.mapq2
        elif self.mapq1 is not None and self.mapq2 is None:
            return self.mapq1
        else:
            return min(self.mapq1, self.mapq2)

    def replace(self, **kwargs):
        d = {field_name: getattr(self, field_name) for field_name in self.field_names}
        d.update(**kwargs)
        return self.__class__(**d)

    @classproperty
    def field_names(cls):
        return cls.__slots__

    @classproperty
    def field_types(cls):
        sig = funcsigs.signature(cls.__init__)
        return [sig.parameters[field].annotation for field in cls.field_names]

    @property
    def tlen(self):
        return self.stop - self.start

    @property
    def length(self):
        return self.tlen

    @property
    def midpoint(self):
        return self.start + (self.stop - self.start) // 2

    def __repr__(self):
        items = ",".join(
            "%s=%r" % (k, getattr(self, k))
            for k in self.field_names
            if getattr(self, k) is not None
        )
        return f"Fragment({items})"

    def __eq__(self, other):
        return [getattr(self, k) for k in self.field_names] == [
            getattr(other, k) for k in self.field_names
        ]


def cigar_fragment_end_matches(cigarstring, is_reverse):
    """
    Returns True of the fragment represented by an alignment's cigarstring matches at the fragment start or end.
    Since SAM alignments are always represented on the forward strand, we must look at the end of the

    >>> cigar_fragment_end_matches('1M', is_reverse=True)
    True
    >>> cigar_fragment_end_matches('10M', is_reverse=False)
    True
    >>> cigar_fragment_end_matches('10S10M', is_reverse=True)
    True
    >>> cigar_fragment_end_matches('10S10M', is_reverse=False)
    False
    """
    if is_reverse:
        return cigarstring.endswith("M")
    else:
        return bool(re.search(r"^\d+M", cigarstring))


def get_methyl_idxs_for_pair_of_alignments(
    align_1, align_2, fasta_file, methylation_prior=0.5
):
    # assumes align_1 has positive tlen
    assert align_1.query_name == align_2.query_name
    assert alignment_to_fragment_start_stop(
        align_1
    ) == alignment_to_fragment_start_stop(align_2)
    assert align_1.tlen > 0
    assert align_1.tlen == -align_2.tlen

    def get_random_methylation(methylation_prior=0.5):
        # In the case of a tie, and equal PHRED scores, we randomly choose methylation with some probability
        return numpy.random.random() < methylation_prior

    frag_start, frag_stop = alignment_to_fragment_start_stop(align_1)
    ref_seq = fasta_file.fetch(align_1.reference_name, frag_start, frag_stop)

    len_1 = len(align_1.seq)
    len_2 = len(align_2.seq)
    quals_1 = align_1.qual
    quals_2 = align_2.qual

    ref_cpg_idxs = numpy.array(
        [match.start() for match in re.finditer("CG", ref_seq)], dtype=int
    )

    # Case 1: perfect overlap between reads (tlen == read length)
    if len_1 == len_2 == align_1.tlen:

        methyl_idxs_1 = numpy.array(
            [align_1.seq[idx : idx + 2] == "CG" for idx in ref_cpg_idxs], dtype=bool
        )
        methyl_idxs_2 = numpy.array(
            [align_2.seq[idx : idx + 2] == "CG" for idx in ref_cpg_idxs], dtype=bool
        )

        # First, assume calls are identical
        methyl_idxs = methyl_idxs_1

        # Calls are different, different phred scores
        mismatches = methyl_idxs_1 != methyl_idxs_2
        if mismatches.any():
            for ii in list(numpy.where(mismatches)[0]):
                # We want to check if we see a G or not
                pos_idx = ref_cpg_idxs[ii] + 1

                # We skip the case where the phred score of read 1 is higher, as that's what the scores start as
                if ord(quals_1[pos_idx]) > ord(quals_2[pos_idx]):
                    assert methyl_idxs[ii] == methyl_idxs_1[ii]
                # If the phred score of read 2 is higher, we set the methylation to read 2's
                if ord(quals_2[pos_idx]) > ord(quals_1[pos_idx]):
                    methyl_idxs[ii] = methyl_idxs_2[ii]
                # If the phred scores are the same, randomly choose a methylation status with a prior
                elif ord(quals_1[pos_idx]) == ord(quals_2[pos_idx]):
                    # Randomly choose methylated or unmethylated.
                    methyl_idxs[ii] = get_random_methylation(methylation_prior)

    # Case 2: some overlap between reads, but not completely disjoint
    elif len_1 + len_2 > align_1.tlen:
        # read 1 |------>|
        # read 2      |<------|
        # CpGs   |  x   x   x |
        # idxs   0 2 4 6 8 0 2
        # ref_cpg_idxs = [3, 7, 11]
        # region_len = 14
        # read_len = 9

        # The positions of CpGs in read 1 (left read)
        # We select the cpg indexes that are less than the length of read 1
        # [3, 7, 11][[3, 7, 11] < 9] = [3, 7, 11][T, T, F] = [3, 7]
        # We subtract the length by 1 to account for the edge case of a CG being split off at the end of read 1,
        # i.e., read 1 ending with a C, which WOULD be followed by a G if it continued
        cpg_idxs_1 = ref_cpg_idxs[ref_cpg_idxs < len_1 - 1]

        # The positions of CpGs in read 2, relative to left side of read 2 (right read)
        # We select the cpg indexes that are in the last `len_2` bases of the region
        # [3, 7, 11][[3, 7, 11] > 14 - 9] - (14 - 9)
        #  = [3, 7, 11][[3, 7, 11] > 5] - (5)
        #  = [3, 7, 11][F, T, T] - (5)
        #  = [7, 11] - 5
        #  = [2, 6]
        start_idx_of_overlap = len(ref_seq) - len_2
        cpg_idxs_2 = (
            ref_cpg_idxs[ref_cpg_idxs >= start_idx_of_overlap] - start_idx_of_overlap
        )

        methyl_idxs_1 = numpy.array(
            [align_1.seq[idx : idx + 2] == "CG" for idx in cpg_idxs_1], dtype=bool
        )
        methyl_idxs_2 = numpy.array(
            [align_2.seq[idx : idx + 2] == "CG" for idx in cpg_idxs_2], dtype=bool
        )
        num_cpgs_overlapping = len(cpg_idxs_1) + len(cpg_idxs_2) - len(ref_cpg_idxs)
        assert num_cpgs_overlapping >= 0

        # Create the idxs where we see methylation
        methyl_idxs = numpy.concatenate(
            [methyl_idxs_1, methyl_idxs_2[num_cpgs_overlapping:]]
        )
        if num_cpgs_overlapping > 0:

            # Correct errors wherever we see a mismatch in methylation status
            for ii in numpy.where(
                methyl_idxs_1[-num_cpgs_overlapping:]
                != methyl_idxs_2[:num_cpgs_overlapping]
            )[0]:
                # We skip the case where the phred score of read 1 is higher
                # We add a 1 because we want the quals of the 'G' in 'CG'
                # ref_cpg_idxs is the indexes of the reference sequence at which we see a CG
                # ref_cpg_idxs_overlap_start is the index at which ref_cpg_idxs starts the overlap
                ref_cpg_idxs_overlap_start = len(cpg_idxs_1) - num_cpgs_overlapping
                pos_idx_1 = ref_cpg_idxs[ref_cpg_idxs_overlap_start + ii] + 1
                pos_idx_2 = cpg_idxs_2[ii] + 1

                # If the quality of read 1 base call is higher, we don't need to do anything, because by default,
                #  we use start with methylation status of read 1 wherever we overlap
                assert (
                    methyl_idxs[ref_cpg_idxs_overlap_start + ii]
                    == methyl_idxs_1[ref_cpg_idxs_overlap_start + ii]
                )

                # If the quality of read 2 base call is higher, we set the methylation status to read 2's status
                if ord(quals_2[pos_idx_2]) > ord(quals_1[pos_idx_1]):
                    methyl_idxs[ref_cpg_idxs_overlap_start + ii] = methyl_idxs_2[ii]

                # If the qualities of the reads is identical, we randomly choose a methylation status
                elif ord(quals_1[pos_idx_1]) == ord(quals_2[pos_idx_2]):
                    methyl_idxs[
                        ref_cpg_idxs_overlap_start + ii
                    ] = get_random_methylation(methylation_prior)

    # Case 3: no overlap (may be just touching). We remove all mentions of a CpG in the non-sequenced insert
    else:
        assert len_1 + len_2 <= align_1.tlen
        # The positions of CpGs in read 1 (left read)
        cpg_idxs_1 = ref_cpg_idxs[ref_cpg_idxs < len_1]
        # The positions of CpGs in read 2, relative to read 2 start (right read)
        cpg_idxs_2 = ref_cpg_idxs[ref_cpg_idxs >= (len(ref_seq) - len_2)] - (
            len(ref_seq) - len_2
        )
        ref_cpg_idxs = numpy.concatenate(
            [
                ref_cpg_idxs[: len(cpg_idxs_1)],
                ref_cpg_idxs[len(ref_cpg_idxs) - len(cpg_idxs_2) :],
            ]
        )

        methyl_idxs_1 = numpy.array(
            [align_1.seq[idx : idx + 2] == "CG" for idx in cpg_idxs_1], dtype=bool
        )
        methyl_idxs_2 = numpy.array(
            [align_2.seq[idx : idx + 2] == "CG" for idx in cpg_idxs_2], dtype=bool
        )

        methyl_idxs = numpy.concatenate([methyl_idxs_1, methyl_idxs_2])

    return ref_cpg_idxs, methyl_idxs


def alignment_to_fragment_start_stop(align):
    # re-enable the assert after we've fixed the is_proper_pair flag
    # assert align.is_proper_pair
    if align.tlen < 0:
        frag_start = align.aend + align.tlen
        frag_stop = align.aend
    elif align.tlen > 0:
        frag_start = align.pos
        frag_stop = align.pos + align.tlen
    else:
        raise AssertionError(f"properly paired reads should not have a tlen of 0")

    return frag_start, frag_stop


def bam_to_align(
    alignment_file: Union[str, pysam.AlignmentFile],
    chrom=None,
    start=None,
    stop=None,
    fasta_file: pysam.FastaFile = None,
    max_tlen=1000,
    include_neg_tlen=False,
    min_mapq=DEFAULT_MIN_MAPQ,
):
    """
    Simplified version of bam_to_fragments which relies on tlen being set properly by the aligner.


    code to confirm that tlen represents where on the reference two reads align to

    align = pysam.AlignmentFile('/ssd/rbio-pipe/SEQRUN2/output/fastq_to_bam/W_1/align.bam')
    aligns = list(itertools.islice(align, 0, 100000))

    for qname, aligns in itertools.groupby(sorted(aligns, key=lambda a: (a.qname, -a.tlen)), lambda a:a.qname):
        aligns = list(aligns)
        left, right = aligns
        all_coords = sorted((left.pos, left.reference_end, right.pos, right.reference_end))
        start, stop = all_coords[0], all_coords[-1]
        aligns[0].tlen
        if stop - start != aligns[0].tlen:
            raise Exception('tlen does not match reference_end')
    """
    if isinstance(alignment_file, str):
        alignment_file = pysam.AlignmentFile(alignment_file)
        close_alignment_file = True
    else:
        close_alignment_file = False

    if isinstance(fasta_file, str):
        fasta_file = pysam.FastaFile(fasta_file)
        close_fasta_file = True
    else:
        close_fasta_file = False

    contig_lengths = {dict_["SN"]: dict_["LN"] for dict_ in alignment_file.header["SQ"]}

    if stop is None and chrom is not None:
        stop = contig_lengths[chrom]

    if start is None and chrom is not None:
        start = 0

    if chrom is None:
        align_iter = iter(alignment_file)
    else:
        align_iter = alignment_file.fetch(
            chrom,
            max(start - max_tlen - 1, 0),
            min(stop + max_tlen + 1, contig_lengths[chrom]),
        )

    def align_is_invalid(align_, min_mapq_):

        is_invalid = (
            align_.tlen == 0
            or abs(align_.tlen) > max_tlen
            or align_.is_qcfail
            or align_.is_supplementary
            or align_.is_duplicate
            or align_.is_unmapped
            or align_.mate_is_unmapped
        )
        if min_mapq_ is not None:
            is_invalid |= align_.mapq < min_mapq_ or (
                align_.has_tag("MQ") and align_.get_tag("MQ") < min_mapq_
            )

        if align_.has_tag("MC"):
            # Check that fragment start/end has a cigar MATCH operator.
            # Otherwise, fragment length will be off.
            # We only do this if MC is set, because otherwise we would get different results depending on whether
            # we used align1 or align2 to do the filtering.
            # It should always be set in internal BAMs, but may not be for external ones.
            is_invalid |= not cigar_fragment_end_matches(
                align_.cigarstring, align_.is_reverse
            ) or not cigar_fragment_end_matches(
                align_.get_tag("MC"), align_.mate_is_reverse
            )

        return is_invalid

    for align in align_iter:
        # only look at read 1 alignments, and use tlen to infer the fragment length
        if include_neg_tlen or 0 < align.tlen <= max_tlen:
            if align_is_invalid(align, min_mapq_=min_mapq):
                continue
            yield align

    if close_alignment_file:
        alignment_file.close()
    if close_fasta_file:
        fasta_file.close()


def bam_to_fragments(
    alignment_file: Union[str, pysam.AlignmentFile],
    chrom=None,
    start=None,
    stop=None,
    fasta_file: pysam.FastaFile = None,
    max_tlen=1000,
    min_mapq=DEFAULT_MIN_MAPQ,
):
    if isinstance(alignment_file, str):
        alignment_file = pysam.AlignmentFile(alignment_file)
        close_alignment_file = True
    else:
        close_alignment_file = False

    if isinstance(fasta_file, str):
        fasta_file = pysam.FastaFile(fasta_file)
        close_fasta_file = True
    else:
        close_fasta_file = False

    align_iter = bam_to_align(
        alignment_file=alignment_file,
        chrom=chrom,
        start=start,
        stop=stop,
        fasta_file=fasta_file,
        max_tlen=max_tlen,
        min_mapq=min_mapq,
    )
    for align in align_iter:
        frag_start, frag_stop = alignment_to_fragment_start_stop(align)

        # chrom is None implies that we aren't doing any filtering. e.g. start and stop are also none
        if chrom is not None and not intervals_intersect_with_none(
            frag_start, frag_stop, start, stop
        ):
            continue

        if fasta_file is not None:
            gc = round(
                get_percent_gc(
                    fasta_file.fetch(align.reference_name, frag_start, frag_stop)
                ),
                5,
            )
        else:
            gc = None

        if align.is_read1:
            strand = "+"
        elif align.is_read2:
            strand = "-"
        else:
            strand = None

        if not align.has_tag("CB"):
            cell_barcode = None
        else:
            cell_barcode = align.get_tag("CB")

        frag = Fragment(
            align.reference_name,
            frag_start,
            frag_stop,
            mapq1=align.mapq,
            mapq2=align.get_tag("MQ") if align.has_tag("MQ") else None,
            gc=gc,
            strand=strand,
            cell_barcode=cell_barcode,
        )

        yield frag

    if close_alignment_file:
        alignment_file.close()
    if close_fasta_file:
        fasta_file.close()
