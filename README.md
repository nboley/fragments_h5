# Fragments H5

Fragments h5 is a python library that implements a fast and memory efficient method for storing DNA sequencing fragments. It was developed for the analysis of cell-free-DNA fragments using fragmentomics approaches. See this paper for a brief overview of the field: https://www.nature.com/articles/s41416-021-01635-z

## Quick Start

First, convert your bam to a fragment h5 using the build-fragments-h5 command. 
```
build-fragments-h5 test/data/small.chr6.bam ./small.chr6.fragments.h5 --sample-id test_data --reference hg38 
```

We can use h5ls to inspect the h5 file.

```
> user@machine: ~/src/fragments_h5/$ h5ls -r ./small.chr6.fragments.h5
/                        Group
/data                    Group
/data/chr6               Group
/data/chr6/lengths       Dataset {677}
/data/chr6/mapq          Dataset {677, 2}
/data/chr6/starts        Dataset {677}
/data/chr6/strand        Dataset {677}
/fragment_length_counts  Dataset {65536}
/index                   Group
/index/chr6              Dataset {34163}
```

Next, open up your favorite python interface to load the fragments h5.
```
Python 3.7.1 | packaged by conda-forge | (default, Feb 18 2019, 01:42:00) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.1.1 -- An enhanced Interactive Python. Type '?' for help.

In [1]: from fragments_h5 import FragmentsH5                                                                                                                                                                                                                                              
In [2]: fh5 = FragmentsH5("./small.chr6.fragments.h5")                                                                                                                                                                                                                                    
```

FragmentsH5 implements a fetch method which returns an iterator over Fragment objects. This is a very slow interface, but can be useful for testing and inspection. Here's an example of grabbing a single fragment:
```
In [3]: next(fh5.fetch())                                                                                                                                                                                                                                                  
Out[3]: Fragment(chrom='chr6',start=634343,stop=634501,mapq1=60,mapq2=60,strand='+')
```

The better way to load fragments is through the fetch_array method. This returns a tuple containing 3 elements: a numpy array of starts, a numpy array of stops, and a dictionary containing supplemental data.
```
In [19]: x = fh5.fetch_array('chr6')                                                                                                                                                                                                                                                      

In [20]: x[0].shape                                                                                                                                                                                                                                                                       
Out[20]: (676,)

In [21]: x[1].shape                                                                                                                                                                                                                                                                       
Out[21]: (676,)

In [22]: x[2]                                                                                                                                                                                                                                                                             
Out[22]: {}
```
The returned starts and ends can then be used to construct an array of midpoints, look at raw starts and ends, or similar. 


## Interface

### FragmentsH5 Class

#### Attributes

```
filename               : return an absolute path to the h5 filename
name                   : alias for filename
ref                    : the name of the reference genome
sample_id              : identifier for the sample that this data originates from (set at creation time)
has_methyl             : whether or not the fragment h5 contains cpg and converted cpg counts
has_strand             : where or not the fragment h5 conmtains strand information
max_fragment_length    : the maximum fragment length stored
fragment_length_counts : an array of fragment counts for each fragment length from 0 to max_fragment_length
```

#### FragmentsH5.init(...)

```
Args:
    fname (str)             : path to fragment h5 file.
    mode (str)              : mode to load the h5 file in. (defaults to 'r')
    cache_pointers (bool)   : Whether or not to load the index into memory. (defaults to False)
        This is a memory-intensive operation but doubles the fetch speed. This is usually worth
        it when you're analyzing a single fragment h5 and want to look across lots of regions,
        but can quickly lead to memory issue.
    use_s3_fs (bool)        : If set then allow `s3:// ... ` paths to be loaded. Requires s3fs to be installed.
    s3_fs_read_timeout (int): Tiemout time for s3fs. (Defaults to 60)

Returns:
    FragmentsH5 instance
```

#### FragmentsH5.fetch_array(...)

```
Fetch all fragments that overlap contig:[region_start, region_stop)

Args:
    contig (str): contig or chromosome to fetch fragments from
    region_start (int, optional) : inclusive start of region to fetch fragments from (defaults to whole contig)
    region_stop (int, optional): exclusive stop of region to fetch fragments from (defaults to whole contig)
    max_frag_len (int, optional): max length of fragments to fetch (defaults to no filtering)
    return_mapqs (bool, optional): return fragments' mapq scores (defaults to False)
    return_gc (bool, optional): return fragments' gc content (defaults to False)
    return_strand (bool, optional): return fragments' strand (defaults to False)
    return_methyl (bool, optional): return fragments' cpg and converted cpg counts (defaults to False)
    filter_to_midpoint_frags (bool, optional): only return fragments whose midpoints are contained in
                                               the filter region (defaults to returning all overlapping
                                               fragments)


Returns:
    starts -> numpy.int32 array containing fragment starts
    stops -> numpy.int32 array containing fragment stops
    supp_data -> {
        (only set if return_mapqs is set to True)
        mapq          -> numpy.int32 array with two columns containing mapq1, mapq2
                         unknown mapqs are set to -1

        (only set if return_gc is set to True)
        gc            -> numpy.float32 array containing fraction of fragment that is a g or c
                         unknown GC's are set to NaN
                         Note: GC is stored in a uint8 so there are only two significant digits

        (only set if return_strand is set to True)
        strand        -> numpy.char1 array containing '+' or '-' for the fragment strand

        (only set if return_gc is set to True)
        num_cpgs      -> numpy.uint8 array containing number of cpgs in the fragment
        num_meth_cpgs -> numpy.uint8 array containing number of converted cpgs in the fragmnet
    }
```

### build-fragments-h5
```
build-fragments-h5 --help
usage: build-fragments-h5 [-h] [--quiet | --verbose | --debug]
                          [--log-format LOG_FORMAT]
                          [--log-filename LOG_FILENAME]
                          [--log-file-verbosity-level {CRITICAL,ERROR,WARNING,INFO,DEBUG,NOTSET}]
                          [--log-file-format LOG_FILE_FORMAT] --reference
                          {hg16,hg17,hg18,hg19,hg38} --sample-id SAMPLE_ID
                          [--fasta FASTA] [--contigs CONTIGS [CONTIGS ...]]
                          [--set-mapq-255-to-none] [--exclude-strand]
                          [--read-methyl]
                          input_bam_or_bed output_frags_h5

positional arguments:
  input_bam_or_bed      bam or bed file to read fragments from
  output_frags_h5       where to write the new fragments h5

optional arguments:
  -h, --help            show this help message and exit
  --quiet, -q           Only output error log messages (and above) to the
                        output stream.
  --verbose, -v         Output info level log messages (and above) to the
                        output stream.
  --debug               Output debug level log messages (and above) to the
                        output stream.
  --reference {hg16,hg17,hg18,hg19,hg38}
                        The reference genome of input_bam (hg19 or hg38).
  --sample-id SAMPLE_ID
                        The sample_id of the bam (should correspond to an
                        entry in SampleDataFrame.
  --fasta FASTA         Path to a fasta file containing the reference genome.
  --contigs CONTIGS [CONTIGS ...]
                        Restrict building the fragment h5 over these contigs.
  --set-mapq-255-to-none
                        set mapqs of 255 to None
  --exclude-strand      Exclude strand info
  --read-methyl         Read in methylation frag beds

logging:
  --log-format LOG_FORMAT
                        Format string to use for log messages.
  --log-filename LOG_FILENAME
                        Write log messages to both the default handler and
                        --log-filename. (default: do not write messages to a
                        log file)
  --log-file-verbosity-level {CRITICAL,ERROR,WARNING,INFO,DEBUG,NOTSET}
                        Logging level to use for the log file handler.
                        (default: log all messages)
  --log-file-format LOG_FILE_FORMAT
                        Format string to use for log messages written to file
                        (see LogRecord in the logging module docs).
```

## Algorithm Description:
This data structure is designed for quickly finding all fragments that overlap an interval.

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

The above is algorithmically fast, but can be slow in practice because it requires that we access the
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

ADDITIONAL PERFORMANCE IMPROVEMENT FROM THE ABOVE DESCRIPTION:

The above algorithm performs two index searches, but it actually appears to be faster to perform an to-memory
copy of the index chosen regions and then do all of the region tuning in the numpy array. That's what we actually
implement in the the class.


## Notes on Profiling:

During development we did some careful profiling of this code path. I think that it would be possible to squeeze another 10% out of this code, but the majority of the time is spent accessing the h5 arrays. h5 isn't particularly fast and mem-mapped arrays could probably give another 2x speedup or so, but it's probably not worth the maintenance cost.

I'm profiling using chr10 from two bams:
BO_4_1.capture.bam
BO_4_1.wgs.bam

This benchmark is a little un-realistic because a lot of the time is loading the data from disk, and I re-run these tests so that the data is cached. However, it still provides something so long as we are aware that, everything else considered, small files can be significantly better.

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

Decreasing the block and index size didn't produce a significant performance change

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


After converting to use read_direct

max frag len is 511:
WGS: Execution time: 1.4718513488769531
Capture: Execution time: 2.2286031246185303

max frag len is 511 and return mapqs and gc:
WGS: Execution time: 3.755711555480957
Capture: Execution time: 5.995445013046265

Index block size change to 5000, and using contiguous arrays
Filesizes:
wgs: 68M
cap: 9.8M

max frag len is 511:
WGS: Execution time: 1.45387601852417
Capture: Execution time: 1.7130119800567627

max frag len is 511 and return mapqs and gc:
WGS: Execution time: 3.130281925201416
Capture: Execution time: 5.176678895950317
