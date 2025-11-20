from collections import Counter

import funcsigs
import logging
import numpy
import re
from dataclasses import dataclass

import pysam

from typing import Union, Optional

from sequence import one_hot_encode_sequences

DEFAULT_MIN_MAPQ = 0


logging.basicConfig(format="%(levelname)s\t%(asctime)-15s\t%(message)s")
log = logging.getLogger(__name__)


@dataclass(slots=True, frozen=True)
class MethylCounts:
    unconverted_cytosines: int
    converted_cytosines: int
    unconverted_cpgs: int
    converted_cpgs: int

    _machine = re.compile("unconverted_cytosines:(\d+);converted_cytosines:(\d+);unconverted_cpgs:(\d+);converted_cpgs:(\d+)")

    def items(self):
        for key in self.__slots__:
            yield key, getattr(self, key)

    @classmethod
    def init_from_yn_tag(cls, tag_str):
        return cls(*map(int, cls._machine.fullmatch(tag_str).groups()))


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


def get_percent_gc_slow(seq):
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
    c = Counter(seq.upper())
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
        "methyl_counts",
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
        methyl_counts: Optional[dict] = None,
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
        return cigarstring.endswith("M") or cigarstring.endswith("=")
    else:
        return bool(re.search(r"^\d+(M|=)", cigarstring))


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

    # if a fasta file has been provided, ensure that any matching contig lengths are the same
    if fasta_file is not None:
        fasta_contig_lens = dict(zip(fasta_file.references, fasta_file.lengths))
        for contig, contig_len in contig_lengths.items():
            if contig in fasta_contig_lens and fasta_contig_lens[contig] != contig_len:
                raise ValueError(
                    f"The bam has contig '{contig}' with a length of '{contig_len}' but the fasta file has a contig length of '{fasta_contig_lens[contig]}' for '{contig}'"
                )

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
    chrom,
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

    if fasta_file is not None:
        contig_seq = fasta_file.fetch(chrom)
        # contig_seq = one_hot_encode_sequences([contig_seq.encode().ravel()])
        # we add an 'a' to the beginning of the seqeunce because the fragment interval is open closed
        # for example, image that we have a contig with sequence ccccaaaaaa (1 c's and 4 a's and 5 c's) and a fragment
        # with a start of 0 and an end of 1 (a one basepair fragment with sequence 'c'). The cumsum is
        #    x =  array([1., 1., 1., 1., 1., 2., 3., 4., 5., 6.], dtype=float32)
        # so x[1] - x[0] == 1 - 1 == 0
        # padding the beginning with zero fixes this edge condition
        seq = one_hot_encode_sequences([('a' + contig_seq).encode()])[0]
        g_or_c_cumsum = seq[:, (1,2)].sum(axis=1).cumsum()

    else:
        contig_seq = None

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
            if frag_stop == frag_start:
                gc = None
            else:
                gc = round(float(g_or_c_cumsum[frag_stop] - g_or_c_cumsum[frag_start])/float(frag_stop - frag_start), 5)
        else:
            gc = None

        if align.is_read1:
            if align.is_forward:
                strand = "+"
            elif align.is_reverse:
                strand = "-"
            else:
                assert False
        elif align.is_read2:
            if align.is_forward:
                strand = "-"
            elif align.is_reverse:
                strand = "+"
            else:
                assert False
        else:
            strand = None

        if align.has_tag("CB"):
            cell_barcode = align.get_tag("CB")
        else:
            cell_barcode = None

        if align.has_tag("YM"):
            methyl_counts = MethylCounts.init_from_yn_tag(align.get_tag("YM"))
        else:
            methyl_counts = None

        frag = Fragment(
            align.reference_name,
            frag_start,
            frag_stop,
            mapq1=align.mapq,
            mapq2=align.get_tag("MQ") if align.has_tag("MQ") else None,
            gc=gc,
            strand=strand,
            cell_barcode=cell_barcode,
            methyl_counts=methyl_counts,
        )

        yield frag

    if close_alignment_file:
        alignment_file.close()
    if close_fasta_file:
        fasta_file.close()

def single_end_bam_to_fragments(
        alignment_file: Union[str, pysam.AlignmentFile],
        chrom=None,
        start=None,
        stop=None,
        fasta_file: pysam.FastaFile = None,
        max_tlen=1000,
        min_mapq=DEFAULT_MIN_MAPQ,
    ):
    assert fasta_file is None

    if chrom is None:
        align_iter = pysam.AlignmentFile(alignment_file)
    else:
        align_iter = pysam.AlignmentFile(alignment_file).fetch(
            chrom,
        )

    for align in align_iter:
        if (
            align.is_qcfail
            or align.is_supplementary
            or align.is_duplicate
            or align.is_unmapped
            or align.mapq < min_mapq
        ): continue

        if align.is_forward:
            strand = "+"
        elif align.is_reverse:
            strand = "-"
        else:
            assert False

        frag = Fragment(
            align.reference_name,
            align.pos,
            align.aend,
            mapq1=align.mapq,
            mapq2=None,
            gc=None,
            strand=strand,
            cell_barcode=None,
            methyl_counts=None,
        )
        yield frag
    return
