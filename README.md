# Fragments H5

Fragments h5 is a python library that implements a fast and memory efficient method for storing DNA sequencing fragments. It was developed for the analysis of cell-free-DNA fragments using fragmentomics approaches. See this paper for a brief overview of the field: https://www.nature.com/articles/s41416-021-01635-z

## Quick Start

First, convert your bam to a fragmnet h5 using the build-fragments-h5 command.


### Command Line Interface

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

### PROFILE NOTES:

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
