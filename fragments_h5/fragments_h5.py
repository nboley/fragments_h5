"""Code to store fragment info in an h5 file.

This is a fast replacement for fragment beds. It gives retrieval times 10x that of a fragment bed.

During development we did some careful profiling of this code path. I think that it would be possible
to squeeze another 10% out of this code, but the majority of the time is spent accessing the h5 arrays.
h5 isn't particularly fast and mem-mapped arrays could probably give another 2x speedup or so, but it's
probably not worth the maintenance cost.


PROFILE NOTES:

I'm profiling using chr10 from two bams:
BO_4_1.capture.bam
BO_4_1.wgs.bam

the benchmark is a little un-realistic because a lot of the time is loading the data from
disk, and I re-run these tests so that the data is cached. However, it still provides something
so long as we are aware that, everything else considered, small files can be significantly better.

Test code:
    regions = [args for args in GREATEST_HITS if args[0] == 'chr10']
    frag_h5 = FragmentsH5(sys.argv[1])
    for i in range(2500):
        for region in regions:
            _ = frag_h5.fetch_array(*region, max_frag_len=frag_len)

Interleaved starts and stops:
int64 start & stop size
CHUNK_SIZE = 50000
INDEX_BLOCK_SIZE = 10000
MIN_MAPQ = 20
MAX_FRAG_LENGTH = 65535

Filesizes:
wgs: 121M
cap: 14M

max frag len is default (65k):
WGS: Execution time: 3.1765377521514893
Capture: Execution time: 7.31060004234314

max frag len is 511:
WGS: Execution time: 1.7233939170837402
Capture: Execution time: 7.120009183883667

Interleaved starts and stops:
int32 start & stop size
CHUNK_SIZE = 50000
INDEX_BLOCK_SIZE = 10000
MIN_MAPQ = 20
MAX_FRAG_LENGTH = 65535

Filesizes:
wgs: 61M
cap: 7.3M

max frag len is default (65k):
WGS: Execution time: 3.1670446395874023
Capture: Execution time: 6.935131549835205

max frag len is 511:
WGS: Execution time: 1.7471604347229004
Capture: Execution time: 6.849348068237305

# Decreasing the block and index size didn't produce a significant
# performance change

Separate starts and stops:
int32 start & unit16 length
CHUNK_SIZE = 50000
INDEX_BLOCK_SIZE = 10000
MIN_MAPQ = 20
MAX_FRAG_LENGTH = 65535

Filesizes:
wgs: 46M
cap: 5.7M

max frag len is default (65k):
WGS: Execution time: 3.6290781497955322
Capture: Execution time: 4.345688343048096

max frag len is 511:
WGS: Execution time: 3.2602548599243164
Capture: Execution time: 4.274357318878174


Separate starts and stops and mapqs and gc:
int32 start & uint16 length & uint8 mapqs and gcs

Also optimized the fetches, which would have improved the prev benchmark as well

CHUNK_SIZE = 50000
INDEX_BLOCK_SIZE = 10000
MIN_MAPQ = 20
MAX_FRAG_LENGTH = 65535

Filesizes:
wgs: 70M
cap: 8.3M

max frag len is default (65k):
WGS: Execution time: 2.004547119140625
Capture: Execution time: 2.5716748237609863

max frag len is 511:
WGS: Execution time: 1.7754967212677002
Capture: Execution time: 2.5291435718536377

max frag len is 511 and return mapqs and gc:
WGS: Execution time: 4.104072332382202
Capture: Execution time: 6.2343056201934814


## After converting to use read_direct

max frag len is 511:
WGS: Execution time: 1.4718513488769531
Capture: Execution time: 2.2286031246185303

max frag len is 511 and return mapqs and gc:
WGS: Execution time: 3.755711555480957
Capture: Execution time: 5.995445013046265

## Index block size change to 5000, and using contiguous arrays
Filesizes:
wgs: 68M
cap: 9.8M

max frag len is 511:
WGS: Execution time: 1.45387601852417
Capture: Execution time: 1.7130119800567627

max frag len is 511 and return mapqs and gc:
WGS: Execution time: 3.130281925201416
Capture: Execution time: 5.176678895950317

"""
import os

import h5py
import numpy
import pandas
import s3fs

import pysam

import fragments_h5._logging as logging
from fragments_h5.fragment import bam_to_fragments, Fragment


logger = logging.getLogger(__name__)


class FragmentsH5:
    """This data structure is designed for quickly finding all fragments that overlap an interval.

    The fragments are stored in two arrays: a sorted array of starts, and then the fragment lengths associated
    with each start. The stops can then be calculated as starts + frag_lens.

    BISECTION ALGORITHM DESCRIPTION:

    The starts array (f['data/{contig}/starts/]) contains a sorted list of all fragment starts. We wish to find
    all fragments that overlap an interval [i_start, i_stop). This corresponds to finding all fragments [f_start, f_stop)
    such that f_stop > i_start and f_start < i_stop. If we know that the maximum fragment length is <= max_frag_len,
    then the search criteria become f_start > i_start - max_frag_len and f_start < i_stop (of course we'll need to do
    a check at the end to filter out any fragments that don't actually overlap).

    This reduces to two queries:
    (q1) We want to find the smallest fragment start >= i_start - max_frag_len - 1. Since the fragment starts are
         sorted, we can use a bisection algorithm to find this start. With numpy, we can run:
            numpy.searchsorted(starts, i_start - max_frag_len - 1, side='left')

    (q2) We want to find the largest fragment start <= i_stop - 1. Since the fragment starts are sorted, we
         can use a bisection algorithm to find this start. With numpy, we can run:
            numpy.searchsorted(starts, i_stop - 1, side='right')

    INDEX DESCRIPTION:

    The above algorithm is algorithmically fast, but can be slow in practice because it requires that we access the
    entire starts array, which essentially means loading it all from disk.  We can improve on this by storing a mapping
    from genomic positions to index locations. For each genomic location that is evenly divisible by the index block size
    we store the index to the first fragment that is greater than or equal to each genomic location.

    For example, if we set the block size to 10000, INDEX[5] will return the index of the first fragment with genomic
    location >= 5*10000. We can modify the above queries as follows:
    (q1)
            search_pos = i_start - max_frag_len - 1
    (q2)
            search_pos = i_stop - 1
    (shared)
            index_lb = INDEX[search_pos//10000]
            index_ub = INDEX[1 + search_pos//10000]
            res = numpy.searchsorted(starts[index_lb:index_ub], search_pos, side=('left' for q1, 'right' for q2))
            # correct for the slice offset
            res += index_lb

    ADDITIONAL PERFORMANCE IMPROVEMENTS:

    The above algorithm performs two index searches, but it actually appears to be faster to perform an to-memory
    copy of the index chosen regions and then do all of the region tuning in the numpy array. That's what we actually
    implement in "fetch" below.
    """

    def __enter__(self) -> "FragmentsH5":
        return self

    def __exit__(self, *args, **kwargs):
        self.close()

    def close(self):
        self._f.close()
        if self._s3fs is not None:
            self._s3fs.close()

    @property
    def name(self):
        """Return the underlying filename's file path"""
        return os.path.abspath(self._f.filename)

    @property
    def filename(self):
        return self.name


    def cache_pointers(self):
        # load the index into memory
        self.index = {}
        for key in self._f["index"]:
            self.index[key] = self._f["index"][key][:]

        # cache all of the keys into python dictionaries. Accessing this info from the h5 is
        # actually pretty expensive, so we cache it all to speed up subsequent lookups
        self.data = {}
        for key in self._f["data"]:
            self.data[key] = {}
            for i_key in self._f["data"][key]:
                self.data[key][i_key] = self._f["data"][key][i_key]

    def __init__(
        self,
        fname,
        mode="r",
        cache_pointers=True,
        use_s3_fs=False,
        s3_fs_read_timeout=60,
    ):
        if use_s3_fs:
            assert fname.startswith("s3://"), fname

            self._s3fs = h5_file = s3fs.S3FileSystem(
                config_kwargs=dict(read_timeout=s3_fs_read_timeout)
            ).open(fname)
        else:
            self._s3fs = None
            h5_file = fname

        self._f = h5py.File(h5_file, mode)

        # set the metadata as attributes
        self.contig_lengths = eval(self._f.attrs["_contig_lengths_str"])
        self.index_block_size = self._f.attrs["index_block_size"]
        self.ref = self._f.attrs["ref"]
        self.max_fragment_length = self._f.attrs["max_fragment_length"]
        self.sample_id = self._f.attrs["sample_id"]
        if "fragment_length_counts" in self._f:
            self.fragment_length_counts = self._f["fragment_length_counts"][:]

        self.index = self._f["index"]
        self.data = self._f["data"]

        if cache_pointers:
            self.cache_pointers()

    @property
    def has_methyl(self):
        return any(
            "num_cpgs" in self.data[contig] for contig in self.contig_lengths.keys()
        )

    @property
    def has_strand(self):
        # A long time ago, Artur (@beyondtheproof) messed up and made fragment h5s for the "small frag" assay
        # with strand taking up two bits. The thought was that one bit would signify "+", the other "-", and if
        # neither were set, it would be a None strand. This was changed to a single bit, that was None if the data
        # didn't exist. Artur never remade the fragment h5s and probably never will. He will forever live in shame.
        return any(
            ("strand" in self.data[contig])
            and (len(self.data[contig]["strand"].shape) == 1)
            for contig in self.contig_lengths.keys()
        )

    def fetch_array(
        self,
        contig,
        region_start=None,
        region_stop=None,
        max_frag_len=None,
        return_mapqs=False,
        return_gc=False,
        return_strand=False,
        return_methyl=False,
        filter_to_midpoint_frags=False,
    ):
        """Fetch all fragments that overlap contig:[region_start, region_stop)

        Returns:
        starts -> numpy.int32 array containing fragment starts
        stops -> numpy.int32 array containing fragment stops
        supp_data -> {
            (only set if return_mapqs is set to True)
            mapq -> numpy.int32 array with two columns containing mapq1, mapq2
                     unknown mapqs are set to -1

            (only set if return_gc is set to True)
            gc    -> numpy.float32 array containing fraction of fragment that is a g or c
                     unknown GC's are set to NaN
                     note that GC is stored in a char internally, so there are only two
                     significant digits

            strand -> numpy.char1 array containing '+' or '-' for the fragment strand
        }

        See the class docstring for details about the datastructure and search algorithm.
        """
        if region_start is None:
            region_start = 0

        if region_stop is None:
            region_stop = self.contig_lengths[contig]

        def empty():
            rv = [numpy.empty(0, "int32"), numpy.empty(0, "int32"), {}]
            if return_mapqs:
                rv[2]["mapq"] = numpy.empty((0, 2), "int32")
            if return_gc:
                rv[2]["gc"] = numpy.empty(0, "float32")
            if return_strand:
                # "+" denotes positive strand
                # "-" denotes negative strand
                rv[2]["strand"] = numpy.empty(0, "|S1")
            if return_methyl:
                rv[2]["num_cpgs"] = numpy.empty(0, "uint8")
                rv[2]["num_meth_cpgs"] = numpy.empty(0, "uint8")
            return rv

        if max_frag_len is None:
            max_frag_len = self.max_fragment_length

        if contig not in self.index:
            return empty()

        # find the lower bound position that needs to be included, accounting for the
        # maximum fragment length. The -1 accounts for the fact that the fragments
        # positions are [closed, open)
        lower_bound_inclusive = max(0, region_start - max_frag_len - 1)
        # find a lower bound for the index of the fragment, using the index
        fragment_lower_bound_index = self.index[contig][
            lower_bound_inclusive // self.index_block_size
        ]
        # find the upper bound position that needs to be included. The -1 accounts
        # for the fact that the region bounds are [closed, open)
        upper_bound_inclusive = max(0, region_stop - 1)
        # find a upper bound for the index of the fragment, using the index
        fragment_upper_bound_index = self.index[contig][
            upper_bound_inclusive // self.index_block_size + 1
        ]

        # if there are no matching fragments
        if fragment_upper_bound_index == fragment_lower_bound_index:
            return empty()

        # find the subset of the array that is known from the index to contain all of the fragments
        sub_array = numpy.empty(
            fragment_upper_bound_index - fragment_lower_bound_index, dtype="int32"
        )
        self.data[contig]["starts"].read_direct(
            sub_array,
            source_sel=numpy.s_[fragment_lower_bound_index:fragment_upper_bound_index],
        )

        start_index = numpy.searchsorted(sub_array, lower_bound_inclusive, side="left")
        stop_index = numpy.searchsorted(sub_array, upper_bound_inclusive, side="right")

        # if there are no matching fragments
        if stop_index == start_index:
            return empty()

        # find all potentially overlapping fragments
        # because we searched over the max_fragment_length, it's likely that we
        # included fragments that don't actually overlap.
        indices_slice = numpy.s_[
            fragment_lower_bound_index
            + start_index : fragment_lower_bound_index
            + stop_index
        ]
        slice_len = stop_index - start_index
        # special case the starts array so that we don't need to do another fetch
        starts = sub_array[start_index:stop_index]

        # use numpy empty and read_direct b/c it's about 10% faster
        lengths = numpy.empty(slice_len, dtype="uint16")
        self.data[contig]["lengths"].read_direct(lengths, source_sel=indices_slice)

        # # there may be fragments that don't overlap in this group b/c we searched on start - max_fragment_length.
        # # filter these out. Note that a fragment doesn't overlap if the fragment stop is <= the region start
        # # or the fragment start is >= than the region stop (b/c regions are (closed, open). This implies that a region
        # # overlaps if (fragment_start < region_stop and fragment_stop > region_start)
        # first build the stops array, whihc we need for the mask
        stops = starts + lengths
        mask = (
            (starts < region_stop) & (stops > region_start) & (lengths <= max_frag_len)
        )

        if filter_to_midpoint_frags:
            lengths = stops - starts
            midpoints = starts - region_start + lengths // 2
            mask &= (midpoints >= 0) & (midpoints < region_stop - region_start)

        starts = starts[mask]
        stops = stops[mask]

        supp_data = {}

        if return_mapqs:
            supp_data["mapq"] = numpy.empty((slice_len, 2))
            self.data[contig]["mapq"].read_direct(
                supp_data["mapq"], source_sel=indices_slice
            )
            # filter out non-overlapping fragments, and cast to an int32
            supp_data["mapq"] = supp_data["mapq"][mask].astype("int32")
            # set the empty mapqs to -1
            supp_data["mapq"][supp_data["mapq"] == 255] = -1

        if return_gc:
            _tmp_gc = numpy.empty(slice_len, dtype="uint8")
            self.data[contig]["gc"].read_direct(_tmp_gc, indices_slice)
            _tmp_gc = _tmp_gc[mask]
            # find all of the 'None' gc's. We need to do this before the float
            # case to avoid any rounding errors
            gc_none_mask = _tmp_gc == 255
            # cast to a float, and re-scale
            _tmp_gc = _tmp_gc.astype("float32", copy=False) / 254.0
            # set all of the "None's" to NaN
            _tmp_gc[gc_none_mask] = numpy.NaN
            supp_data["gc"] = _tmp_gc

        if return_strand:
            if "strand" not in self.data[contig]:
                raise ValueError("The referenced h5 file does not contain strand info.")
            supp_data["strand"] = numpy.empty(slice_len, dtype="|S1")
            self.data[contig]["strand"].read_direct(
                supp_data["strand"], source_sel=indices_slice
            )
            supp_data["strand"] = supp_data["strand"][mask]

        if return_methyl:
            if ("num_cpgs" not in self.data[contig]) or (
                "num_meth_cpgs" not in self.data[contig]
            ):
                raise ValueError(
                    "The referenced h5 file does not contain methylation info."
                )

            supp_data["num_cpgs"] = numpy.empty(slice_len, dtype="uint8")
            self.data[contig]["num_cpgs"].read_direct(
                supp_data["num_cpgs"], source_sel=indices_slice
            )
            supp_data["num_cpgs"] = supp_data["num_cpgs"][mask]

            supp_data["num_meth_cpgs"] = numpy.empty(slice_len, dtype="uint8")
            self.data[contig]["num_meth_cpgs"].read_direct(
                supp_data["num_meth_cpgs"], source_sel=indices_slice
            )
            supp_data["num_meth_cpgs"] = supp_data["num_meth_cpgs"][mask]

        return starts, stops, supp_data

    def fetch(
        self,
        contig=None,
        start=None,
        stop=None,
        max_frag_len=None,
        return_mapqs=True,
        return_gc=True,
        return_strand=True,
        return_methyl=False,
    ):
        """Iterate over all fragments that overlap contig:start-stop.

        If contig is None, then iterate over all contigs. If start or stop are None then
        iterate from chrom start and/or chrom end.

        Note: max_frag_len is set to 65000 by default. Decreasing that value can lead to
        a significant performance improvement.

        .. warning:: this is MUCH slower than using .fetch_array()
        """
        if contig is None:
            assert start is None and stop is None
            contigs = list(self.contig_lengths.keys())
        else:
            contigs = [contig]
        del contig  # make sure that we don't re-use this

        def _negative_1_to_none(x):
            # used for cleaning up mapq scores
            # we store "None" as -1 b/c ints don't support Nones
            if x == -1:
                return None
            return x

        def _nan_to_none(x):
            # used for cleaning up GC scores
            # we store "None" as NaN b/c floats don't support Nones
            if numpy.isnan(x):
                return None
            return x

        def _optional_binary_to_strand(x):
            if x == b"-":
                return x
            if x == b"+":
                return x
            if x is True:
                return b"-"
            elif x is False:
                return b"+"
            elif x is None:
                return None
            else:
                raise ValueError(f"Unrecognized stand value '{x}'")

        for contig in contigs:
            contig_stop = self.contig_lengths[contig] if stop is None else stop
            contig_start = 0 if start is None else start

            frag_starts, frag_stops, supp_data = self.fetch_array(
                contig,
                contig_start,
                contig_stop,
                max_frag_len=max_frag_len,
                return_mapqs=return_mapqs,
                return_gc=return_gc,
                return_strand=return_strand,
                return_methyl=return_methyl,
            )

            for frag_start, frag_stop, mapq, gc, strand, num_cpgs, num_meth_cpgs in zip(
                frag_starts,
                frag_stops,
                supp_data["mapq"],
                supp_data["gc"],
                supp_data.get("strand", [None] * len(frag_starts)),
                supp_data.get("num_cpgs", [None] * len(frag_starts)),
                supp_data.get("num_meth_cpgs", [None] * len(frag_starts)),
            ):
                yield Fragment(
                    chrom=contig,
                    start=frag_start,
                    stop=frag_stop,
                    mapq1=_negative_1_to_none(mapq[0]),
                    mapq2=_negative_1_to_none(mapq[1]),
                    gc=_nan_to_none(gc),
                    strand=_optional_binary_to_strand(strand),
                    num_cpgs=num_cpgs,
                    num_meth_cpgs=num_meth_cpgs,
                )

    def fetch_counts(self, *args, **kwargs):
        """
        Returns the total number of fragments that overlap contig:[region_start, region_stop)
        :param args: Same input as fetch_array
        :param kwargs: Same input as fetch_array
        :return: The number of counts in a region
        """
        # FIXME: This returns the number of fragments that overlap the region, but not necessarily
        #  have a midpoint in the region. Some fragments may be considered outside the region.
        starts, _, _ = self.fetch_array(*args, **kwargs)
        return len(starts)

    def _add_fragment_length_counts(self):
        """Calculate the fragment length distribution and save it in the h5 file.

        """
        logger.info("Adding fragment length counts")
        fragment_lengths = numpy.zeros(self.max_fragment_length + 1)
        for contig in self._f["data"]:
            fls, counts = numpy.unique(
                self._f["data"][contig]["lengths"], return_counts=True
            )
            for fl, count in zip(fls, counts):
                fragment_lengths[fl] += count
        self._f["fragment_length_counts"] = fragment_lengths


def build_fragments_h5(
    input_fname,
    ofname,
    sample_id,
    reference,
    fasta_file=None,
    allowed_contigs=None,
    set_mapq_255_to_none=False,
    # file containing information about the cell barcodes
    read_strand=True,
    read_methyl=False,
):
    """Write a fragments h5 from the fragments in input_fname to ofname'
    input_fname can be either a bam file or a fragments tsv / bed

    set_mapq_255_to_none:
    By the SAM spec mapq 255 is to be reserved for unknown MAPQs, but not all of
    aligners follow this so we raise an error by default when we encounter a MAPQ of
    255. This setting overrides that.
    """
    INDEX_BLOCK_SIZE = 5000

    # the fragment lengths are stored as a uint16, so the max fragment length is
    # 2**16-1 == 65536-1 == 65535
    MAX_FRAG_LENGTH = 65535
    # this parameter is for dynamically resizing the array as it grows
    CHUNK_SIZE = 1000000

    if fasta_file is not None:
        if isinstance(fasta_file, (str, bytes)):
            fasta_file = pysam.FastaFile(fasta_file)
        elif isinstance(fasta_file, pysam.FastaFile):
            pass
        else:
            raise TypeError(
                "'in_fasta_file' type must be a pysam.FastaFile, str, bytes, or None if GC calcs aren't needed."
            )
        read_gc = True
    else:
        read_gc = False

    f = h5py.File(ofname, "x")

    # Add the attributes
    f.attrs["index_block_size"] = INDEX_BLOCK_SIZE
    f.attrs["ref"] = reference
    f.attrs["max_fragment_length"] = MAX_FRAG_LENGTH
    f.attrs["sample_id"] = sample_id

    # save metadata about the bam into the h5.
    # In particular, find the contigs and lengths, and save them into the h5.
    # h5s dont support dictionaries natively, so
    # we save the dictionary as a string and then use python's eval to load it.
    assert input_fname.endswith(".bam")

    logger.info("Extracting contig lengths from bam")
    with pysam.AlignmentFile(input_fname) as bam_fp:
        f.attrs["_bam_header"] = str(bam_fp.header)
        if allowed_contigs is None:
            allowed_contigs = bam_fp.references
        f.attrs["_contig_lengths_str"] = str(dict(zip(allowed_contigs, bam_fp.lengths)))
        num_mapped = {
            x.contig: x.mapped // 2
            for x in bam_fp.get_index_statistics()
            if allowed_contigs is None or x.contig in allowed_contigs
        }
        num_mapped_cnt = sum(num_mapped.values())
    input_to_fragments = bam_to_fragments
    input_type = "bam"

    contig_lengths = eval(f.attrs["_contig_lengths_str"])
    logger.debug(f"Processing contigs: '{contig_lengths.keys()}'")

    count = 0
    logger.info("Loading fragments for insertion into h5")
    for contig_i, (contig, contig_length) in enumerate(contig_lengths.items()):
        logger.info(f"Converting {contig} ({contig_i+1}/{len(contig_lengths)})")
        # initialize all of the storage arrays
        # we do this in numpy arrays that we copy into the h5 for performance reasons
        starts_arr = numpy.zeros(0, dtype="int32")
        lengths_arr = numpy.zeros(0, dtype="uint16")
        mapq_arr = numpy.zeros((0, 2), dtype="uint8")
        gc_arr = numpy.zeros(0, dtype="uint8")
        strand_arr = numpy.zeros(0, dtype="|S1")
        # We record the number of CpGs not in the whole fragment, but the overlap of the pair of sequenced reads
        # Each read has a maximum length of 150, for a total of 300. CpGs take 2 bp, so max number is 150 < 256
        num_cpgs_arr = numpy.zeros(0, dtype="uint8")
        num_meth_cpgs_arr = numpy.zeros(0, dtype="uint8")

        # if there aren't any fragments then ff is never set. This fixes that edge case.
        ff = 0
        for fragment in input_to_fragments(
            input_fname, chrom=contig, max_tlen=MAX_FRAG_LENGTH, fasta_file=fasta_file
        ):
            # if we have filled this chunk then resize the array
            if ff % CHUNK_SIZE == 0:
                count += CHUNK_SIZE
                logger.debug(
                    f"Finished processing {sample_id} through position "
                    f"{count}/{num_mapped_cnt} "
                    f"({round(100*count/num_mapped_cnt, 2)})"
                )
                if input_type == "bam":
                    logger.debug(
                        f"Finished processing {sample_id} through position "
                        f"{count}/{sum(num_mapped.values())} "
                        f"({round(100*count/sum(num_mapped.values()), 2)})"
                    )
                if input_type == "bed":
                    logger.debug(
                        f"Finished processing {sample_id} through position {count}"
                    )
                starts_arr.resize(starts_arr.shape[0] + CHUNK_SIZE)
                lengths_arr.resize(lengths_arr.shape[0] + CHUNK_SIZE)
                mapq_arr.resize((mapq_arr.shape[0] + CHUNK_SIZE, 2))
                if read_gc:
                    gc_arr.resize(gc_arr.shape[0] + CHUNK_SIZE)
                if read_strand:
                    strand_arr.resize(strand_arr.shape[0] + CHUNK_SIZE)
                if read_methyl:
                    num_cpgs_arr.resize(num_cpgs_arr.shape[0] + CHUNK_SIZE)
                    num_meth_cpgs_arr.resize(num_meth_cpgs_arr.shape[0] + CHUNK_SIZE)

            # set all of the data
            starts_arr[ff] = fragment.start
            lengths_arr[ff] = fragment.length

            if fragment.mapq1 == 255:
                if set_mapq_255_to_none:
                    fragment.mapq1 = None
                else:
                    assert False, (
                        f"{fragment} has MAPQ1 of 255. Make sure this is as intended, "
                        f"and retry with set_mapq_255_to_none=True"
                    )
            if fragment.mapq2 == 255:
                if set_mapq_255_to_none:
                    fragment.mapq2 = None
                else:
                    assert False, (
                        f"{fragment} has MAPQ2 of 255. Make sure this is as intended,"
                        f" and retry with set_mapq_255_to_none=True"
                    )

            mapq_arr[ff, 0] = 255 if fragment.mapq1 is None else fragment.mapq1
            mapq_arr[ff, 1] = 255 if fragment.mapq2 is None else fragment.mapq2

            if read_gc:
                if not (pandas.isnull(fragment.gc) or 0 <= fragment.gc <= 1):
                    raise ValueError(f"{fragment.gc} is invalid for {fragment}")
                gc_arr[ff] = (
                    255 if pandas.isnull(fragment.gc) else int(round(fragment.gc * 254))
                )

            if read_strand:
                strand_arr[ff] = fragment.strand

            if read_methyl:
                num_cpgs_arr[ff] = fragment.num_cpgs
                num_meth_cpgs_arr[ff] = fragment.num_meth_cpgs

            # Increment number of fragments once the fragment is processed
            ff += 1

        # move the data into the h5 file
        f.create_dataset(
            f"data/{contig}/starts", data=starts_arr[: ff + 1], dtype="int32"
        )
        assert MAX_FRAG_LENGTH <= 2 ** 16 - 1
        f.create_dataset(
            f"data/{contig}/lengths", data=lengths_arr[: ff + 1], dtype="uint16"
        )
        f.create_dataset(f"data/{contig}/mapq", data=mapq_arr[: ff + 1], dtype="uint8")
        if read_gc:
            f.create_dataset(f"data/{contig}/gc", data=gc_arr[: ff + 1], dtype="uint8")
        if read_strand:
            f.create_dataset(
                f"data/{contig}/strand", data=strand_arr[: ff + 1], dtype="|S1"
            )
        if read_methyl:
            f.create_dataset(
                f"data/{contig}/num_cpgs", data=num_cpgs_arr[: ff + 1], dtype="uint8"
            )
            f.create_dataset(
                f"data/{contig}/num_meth_cpgs",
                data=num_meth_cpgs_arr[: ff + 1],
                dtype="uint8",
            )

    logger.info("Creating index")
    contig_lengths = eval(f.attrs["_contig_lengths_str"])
    # Build the index
    # See the class notes for a description of the index
    for contig, contig_length in contig_lengths.items():
        block_indices = numpy.array(list(range(0, contig_length, INDEX_BLOCK_SIZE)))
        index_poss = numpy.searchsorted(
            f[f"data/{contig}/starts"][:], block_indices, side="left"
        )
        # add a final entry so that the index lookup still works even if we're at the end of the contig.
        # (e.g. index_ub = INDEX[1 + search_pos//10000] could be out of bounds)
        index_poss = numpy.append(index_poss, len(f[f"data/{contig}/starts"]))
        f[f"index/{contig}"] = index_poss

    f.close()

    # open an h5 using the class interface in r/w mode so that we can
    # add the fragment length information
    fm_h5 = FragmentsH5(ofname, "r+")
    fm_h5._add_fragment_length_counts()
    fm_h5.close()
